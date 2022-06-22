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
class extractor:

    def __init__(self, covariates, grid, aoi, scale, dd, workers, target = None, spcvGridSize = None):
        """
        Prepares explanatoy variables and response variables for data exatraction.
        A grid for spatial cross validation (spcv) is also added as a band.
        The directory to be used to save downloaded data is set.
        
        Args:
            covariates (ee.List, ee.Image, ee.ImageCollection): List should contain images (ee.Image). If imagecollection,
                each image should be a single covariate.

            grid (ee.Image or ee.FeatureCollection): A grid created using the createGrid function since the unique patch
                id's need to match grid values.

            aoi (ee.Feature): The area of interest that designates the study area.

            scale (int): The spatial resolution for exported data.

            dd (string): The destination to save exported/downloaded data.

            target (ee.Image): The response variable. Optionally, this can be included in the covariates list.

            spcvGridSize (int): The size of blocks to use for spcv. All samples within a block are 
                assigned the same unique id.
            
        
        Returns:
            ee objects required for data extraction.

            covariates (ee.Image or ee.ImageCollection): ee.Image if single band image or ee.ImageCollection
                if multiband image or multiple covariates.
            
            target (ee.Image): The target

            grid (ee.Image or ee.FeatureCollection): used as workers during data extraction

            aoi (ee.Feature): The area of interest

            scale (int): The scale (GSD) of the output (covariate and/or target) in metres.
        """
        self.covariates = covariates
        self.grid = grid
        self.aoi = aoi
        self.scale =scale
        self.dd = dd
        self.target = target
        self.spcvGridSize = spcvGridSize
        self.workers = workers
        
        if self.spcvGridSize is not None:
            #Generate Grid for spcv
            self.spcvGrid, _ = createGrid(self.spcvGridSize, aoi = self.aoi, vect = False)

            #Compile covariates
            if self.covariates.name() in ['Image','ImageCollection']:
                finalCovariates = self.spcvGrid.addBands(ee.ImageCollection(self.covariates).toBands()).float()
            elif self.covariates.name() == 'List':
                for covariate in self.covariates:
                    self.spcvGrid = self.spcvGrid.addBands(ee.Image(covariate))
                finalCovariates = self.spcvGrid.float()
            else:
                raise TypeError("Only takes single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
        else:
            #Compile covariates
            if self.covariates.name() in ['Image','ImageCollection']:
                finalCovariates = ee.ImageCollection(self.covariates).toBands().float()
            elif self.covariates.name() == 'List':
                for self.covariate in covariates:
                    finalCovariates = ee.Image(covariate).float()
            else:
                raise TypeError("Only takes single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
        self.covariates = finalCovariates

        #Accomodate target
        if self.target is not None:
            self.target = self.target.rename('target').float()

    def extractAoi(self, index, fid):

        """
        Extract image patches of covariates within an aoi(specified in prepareForExtraction).
        
        Args:
            index (int): index returned by enumerate
            fid (int): unique id of grid. Each grid cell corresponds to a worker.
            
        Returns:
            Data (image patches as GeoTiff) exported to download directory (dd) specified in prepareForExtraction function.
        """
        assert ('linux' in sys.platform), "This code runs on Linux only."
        # Select patch
        if self.grid.name() == 'Image':
            patchSize = self.grid.projection().nominalScale()
            cell = self.grid.updateMask(self.grid.eq(fid).selfMask())
            vCell = cell.addBands(cell).reduceToVectors(ee.Reducer.mean(), self.aoi, patchSize, maxPixels=1e13)
        elif self.grid.name() == 'FeatureCollection':
            vCell = ee.Feature(self.grid.filter(ee.Filter.eq('id', fid)).first())
      
        # Download patch as tif ()
        if self.target is not None:
            #Set directory
            if os.path.exists(self.dd):
                os.chdir(self.dd+'X/')
            else:
                os.makedirs(self.dd+'X/')
            geemap.ee_export_image(self.covariates, crs= 'EPSG:4326', filename= os.path.join(dd, f"X_{fid}.tif"),\
                                scale = self.scale, region= ee.Feature(vCell.first()).geometry(), file_per_band=False)

            #Set directory
            if os.path.exists(self.dd):
                os.chdir(self.dd+'Y/')
            else:
                os.makedirs(self.dd+'Y/')
            # Download patch as tif ()
            geemap.ee_export_image(self.target, crs= 'EPSG:4326', filename= os.path.join(dd, f"Y_{fid}.tif"),\
                                scale = self.scale, region= ee.Feature(vCell.first()).geometry(), file_per_band=False)
        else:
            geemap.ee_export_image(self.covariates, crs= 'EPSG:4326', filename= os.path.join(dd, f"X_{fid}.tif"),\
                                scale = self.scale, region= ee.Feature(vCell.first()).geometry(), file_per_band=False)
    #extract data at points
    def extractPoints(index, fid):
        """
        Extract covariate data at points.
        
        Args:
            index (int): index returned by enumerate
            fid (int): unique id of grid. Each grid cell corresponds to a worker.
            
        Returns:
            Data (csv) exported to download directory (dd) specified in prepareForExraction function.
        """
        assert ('linux' in sys.platform), "This code runs on Linux only."
        # Select patch
        if grid.name() == 'Image':
            patchSize = grid.projection().nominalScale()
            cell = grid.updateMask(grid.eq(fid).selfMask())
            vCell = cell.addBands(cell).reduceToVectors(ee.Reducer.mean(), aoi, patchSize)
        elif grid.name() == 'FeatureCollection':
            vCell = ee.Feature(grid.filter(ee.Filter.eq('id', fid)).first())
        
        #Filter points to patch
        locations = points.filterBounds(vCell)
        
        #Extract data
        size = locations.size().getInfo()
        if  size == 0:
            print(f'Skipping block {index} with zero features')
        else:
            print(f'Extracting block {index} data for {size} gedi point(s)')
            data = covariates.reduceRegions(locations, ee.Reducer.first(), scale)
            properties = covariates.propertyNames()

            output = data.map(lambda ft: ft.set('output', properties.map(lambda prop: ft.get(prop))))
            result = output.aggregate_array('output').getInfo()

            # Write the results to a file.
            filename = f'block_{fid}_results_{index}.csv'
            with open(filename, 'w') as f:
                writer = csv.writer(f)
                # write the header
                writer.writerow(properties.getInfo())
                # write multiple rows
                writer.writerows(result)

    #Extract data around points
    def extractSparsePatches(index, fid):
        """
        Extract summary statistics for neighbourhoods around points.
        
        Args:
            index (int): index returned by enumerate
            fid (int): unique id of grid. Each grid cell corresponds to a worker.
            
        Returns:
            Data (GeoTif) exported to download directory (dd) specified in prepareForExraction function.
            Band 1 cooresponds to the spcvGrid and the last band corresponds to the target variable.
            
        """
        assert ('linux' in sys.platform), "This code runs on Linux only."
        # Select patch
        if grid.name() == 'Image':
            #patchSize = grid.projection().nominalScale()
            cell = grid.updateMask(grid.eq(fid).selfMask())
            vCell = cell.addBands(cell).reduceToVectors(**{'reducer':ee.Reducer.first()\
                                    , 'geometry':aoi, 'scale': scale, 'maxPixels':1e13})
        elif grid.name() == 'FeatureCollection':
            vCell = ee.Feature(grid.filter(ee.Filter.eq('id', fid)).first())
        
            target = finalCovariates.select('target').clip(vCell)
            vTarget = target.int().addBands(target).reduceToVectors(ee.Reducer.last(), vCell, 25)#fix scale value
            size = vTarget.size().getInfo()
        
        if size==0:
            print(f'Skipping patch {fid} with zero features')
        else:
            print(f'Extracting patch {fid} with {size} features')
        # Download patch as tif ()
        geemap.ee_export_image(finalCovariates, filename= os.path.join(dd, f"covs_{fid}.tif"),\
                crs= 'EPSG:4326', scale= scale, region= vTarget.geometry(), file_per_band=False)

    #extract data at points
    def extractPointsNeighbourhood(index, fid):
        """
        Extract covariate data at points.
        
        Args:
            index (int): index returned by enumerate
            fid (int): unique id of grid. Each grid cell corresponds to a worker.
            
        Returns:
            Data (csv) exported to download directory (dd) specified in prepareForExraction function.
        """
        assert ('linux' in sys.platform), "This code runs on Linux only."
        # Select patch
        if grid.name() == 'Image':
            patchSize = grid.projection().nominalScale()
            cell = grid.updateMask(grid.eq(fid).selfMask())
            vCell = cell.addBands(cell).reduceToVectors(**{'reducer':ee.Reducer.first()\
                                    , 'geometry':aoi, 'scale': patchSize, 'maxPixels':1e13})
        elif grid.name() == 'FeatureCollection':
            vCell = ee.Feature(grid.filter(ee.Filter.eq('id', fid)).first())
        
        #Create neighbourhood around points
        target = finalCovariates.select('target')
        points = target.sample(region = vCell, scale = patchSize)\
                .map(lambda x: x.buffer(12.5).bounds())
        
        #Extract data
        size = points.size().getInfo()
        if  size == 0:
            print(f'Skipping block {index} with zero features')
        else:
            print(f'Extracting block {index} data with {size} features')
            # Download patch as tif ()
            geemap.ee_export_image(finalCovariates.updateMask(target), filename= os.path.join(dd, f"covs_{fid}.tif"),\
                crs= 'EPSG:4326', scale= scale, region= vCell.geometry(), file_per_band=False)

    @retry(tries=10, delay=1, backoff=2)
    def pExtract(extractFunction, workers, nProcesses = 25):
        """
        Extract data from gee using parrallel processing
        
        Args:
            extractFunction : One of extractAoi, extractPoints, extractSparsePoints or extractPointNeighbourhood
            
            nProcesses (int): The Number of parrallel processes
            
            workers (list): A list of unique ids corresponding to unique grid id's
            
        Returns:
            Downloads files (csv or Geotiff) to Google Drive
        """

        if __name__ == '__main__':
            assert ('linux' in sys.platform), "This code runs on Linux only."
            logging.basicConfig()

            pool = multiprocessing.Pool(nProcesses)
            pool.starmap(extractFunction, enumerate(workers))

            pool.close()
            pool.join()