import ee
import geemap

import logging
import multiprocessing
import requests
import shutil
from retry import retry
import os
import csv
import sys

import retry

from .utils import createGrid

class extractor:

    def __init__(self, covariates, aoi, scale, dd, target = None, spcvGridSize = None):
        """
        Prepares explanatoy variables and response variables for data exatraction.
        A grid for spatial cross validation (spcv) is also added as a band.
        The directory to be used to save downloaded data is set.
        
        Args:
            covariates (ee.List, ee.Image, ee.ImageCollection): List should contain images (ee.Image). If imagecollection,
                each image should be a single covariate.

            aoi (ee.Feature): The area of interest that designates the study area.

            scale (int): The spatial resolution for exported data.

            dd (string): The destination to save exported/downloaded data.

            target (ee.FeatureCollection, ee.Image): The response variable. Optionally, this can be included in the covariates list.
                If ee.FeatureCollection, covariates will be sampled at all target pixels.

            spcvGridSize (int): The size of blocks to use for spcv. All samples within a block are 
                assigned the same unique id.
        
        Returns:
            ee objects required for data extraction.

            covariates (ee.Image or ee.ImageCollection): ee.Image if single band image or ee.ImageCollection
                if multiband image or multiple covariates.
            
            target (ee.Image): The target
            
        """
        self.covariates = covariates
        self.aoi = aoi
        self.scale =scale
        self.dd = dd
        self.target = target
        self.spcvGridSize = spcvGridSize
        
        # grid for spcv
        if self.spcvGridSize is not None:
            #Generate Grid for spcv
            self.spcvGrid, self.ids = createGrid(self.spcvGridSize, aoi = self.aoi, vect = False)

            #Compile covariates
            if self.covariates.name() in ['Image','ImageCollection']:
                bandnames = self.covariates.bandNames().getInfo()
                finalCovariates = self.spcvGrid.addBands(ee.ImageCollection(self.covariates).toBands()).float().rename(['id'] + bandnames)
            elif self.covariates.name() == 'List':
                finalCovariates = self.spcvGrid.addBands(ee.ImageCollection.fromImages(self.covariates).toBands().float())
            else:
                raise TypeError("Only takes single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
        else:
            #Compile covariates
            if self.covariates.name() in ['Image','ImageCollection']:
                finalCovariates = ee.ImageCollection(self.covariates).toBands().float()
            elif self.covariates.name() == 'List':
                finalCovariates = ee.ImageCollection.fromImages(self.covariates).toBands().float()
            else:
                raise TypeError("Only takes single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
        
        # single or multi band image
        self.covariates = finalCovariates

        #Accomodate target
        if self.target is not None and self.target.name() == 'ee.Image':
            self.target = self.target.rename('target').float().clip(self.aoi)

    def extractAoi(self):

        """
        Extract image patches of covariates within an aoi(specified in prepareForExtraction).
                    
        Returns:
            Default: Data (image patches as GeoTiff) exported to download directory (dd).
            
        """
        
        # Download patch as tif ()
        if not os.path.exists(self.dd + '/X/'):
                os.makedirs(self.dd + '/X/')
        os.chdir(self.dd + '/X/')
        geemap.download_ee_image(self.covariates, crs= 'EPSG:4326', filename= os.path.join(self.dd, f"X/X.tif"),\
                                scale = self.scale, region= self.aoi.geometry())

        if self.target is not None:
            #Set directory
            if not os.path.exists(self.dd + '/Y/'):
                os.makedirs(self.dd + '/Y/')
            os.chdir(self.dd + '/Y/')
            # Download patch as tif
            geemap.download_ee_image(self.target, crs= 'EPSG:4326', filename= os.path.join(self.dd, f"Y/Y.tif"),\
                                scale = self.scale, region= self.aoi.geometry())

        
    def extractPoints(self, batchSize, prefix = 'task_'):
        """
        Extract covariate data at points when points is set to true in extractor.
        
       Args:
           batchSize (int): The number of batches to split job into
           prefix (str): Defaults to 'task_'. Used as file names.
       Returns:
           If points = True, data (csv) exported to download directory (dd).
        """
        #Set working directory
        if not os.path.exists(self.dd):
                os.makedirs(self.dd)
        os.chdir(self.dd)
        
        # Sample all pixels at points
        if self.target.name() == 'Image':
            projection = self.target.projection()
            # Extract data
            points = self.target.sample(**{'dropNulls': True, 'factor': None,\
                                'numPixels': None,  'region': self.aoi,\
                            'scale': self.scale,'projection': projection,'geometries':True})
        else:
            points = self.target.filterBounds(self.aoi)
        
        size = points.size().getInfo()
        pointsList = points.toList(size)
        
        batchSize = batchSize
        
        for i, batch in enumerate(range(0, size, batchSize)):
            fc = ee.FeatureCollection(pointsList.slice(i, i+batchSize))
            task_no = str(i+1)
            print(f'Extracting data for task {i+1} with {fc.size().getInfo()} point(s)')
            tsk_name = prefix + task_no

            data = self.covariates.reduceRegions(fc, ee.Reducer.first(), self.scale)
            properties = data.first().propertyNames()

            output = data.map(lambda ft: ft.set('output', properties.map(lambda prop: ft.get(prop))))
            result = output.aggregate_array('output').getInfo()

            # Write the results to a file.
            filename = f'{tsk_name}.csv'
            with open(filename, 'w') as f:
                writer = csv.writer(f)
                # write the header
                writer.writerow(properties.getInfo())
                # write multiple rows
                writer.writerows(result)
            
    def extractRegions(self, reduce = True):
        """
        Extract summary statistics of covariates for regions.
        
            
        Returns:
            Data (csv) exported to download directory (dd).
            
        """