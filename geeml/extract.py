import ee
import geedim as gd

import logging
import requests
import shutil
import os
import sys

from .utils import createGrid

import urllib.request
import os
import pathlib
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from itertools import product
import csv
from csv import writer

from tqdm import TqdmWarning
from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

class extractor:

    def __init__(self, covariates, aoi, scale, dd, crs = 'EPSG:4326', target = None, spcvGridSize = None, num_threads = 5):
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

            crs (str): The coordinate reference system for exported data.

            target (ee.FeatureCollection, ee.Image): The response variable. Optionally, this can be included in the covariates list.
                If ee.FeatureCollection, covariates will be sampled at all target pixels.

            spcvGridSize (int): The size of blocks to use for spcv. All samples within a block are 
                assigned the same unique id.
                
            num_threads (int): The number of threads to use simulataneously to extract data.
        
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
        self.crs = crs
        self.target = target
        self.spcvGridSize = spcvGridSize
        self.num_threads = num_threads
        
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
        gd.download.BaseImage(self.covariates).download(crs= self.crs, filename= os.path.join(self.dd, f"X/X.tif"),\
                                scale = self.scale, region= self.aoi.geometry(), num_threads= self.num_threads, dtype= 'float64')

        if self.target is not None:
            #Set directory
            if not os.path.exists(self.dd + '/Y/'):
                os.makedirs(self.dd + '/Y/')
            os.chdir(self.dd + '/Y/')
            # Download patch as tif
            gd.download.BaseImage(self.target).download(crs= self.crs, filename= os.path.join(self.dd, f"Y/Y.tif"),\
                                scale = self.scale, region= self.aoi.geometry(), num_threads= self.num_threads, dtype= 'float64')
            

    def geomFeatures(self, grid, item):
        """
        If target Sparse Image, Gets points within a single grid geometry that correspond to item label.
        If target is FeatureCollection, gets features within a single grid geometry that correspond to item label.
        Uses FilterBounds

        
        Args:
            grid (ee.FeatureCollection): The grid to parralelise over
            item (int): corresponds to an id within the grid (from createGrid object).
        
        Returns:
            points that are within the specified item
        
        """
        geom = ee.Feature(grid.filter(ee.Filter.eq('label', item)).first()).geometry()

        # Sample all pixels at points
        if self.target.name() == 'Image':
            projection = self.target.projection()
            # Extract data
            points = self.target.clip(geom).sample(**{'dropNulls': True, 'factor': None,\
                                'numPixels': None,  'region': geom,\
                            'scale': self.scale,'projection': projection,'geometries':True})
        else:
            points = self.target.filterBounds(geom)
        
        return points

    def geomGridCells(self, grid, item):
        """
        Convert target grid cells that intersect item grid cell to featurecollection
        
        
        Args:
            grid (ee.FeatureCollection): The grid to parralelise over
            item (int): corresponds to an id within the grid (from createGrid object).
        
        Returns:
            gridcells that are within the specified item
        
        """
        geom = ee.Feature(grid.filter(ee.Filter.eq('label', item)).first()).geometry()

        reduction = self.target.clip(geom).reduceRegion(ee.Reducer.frequencyHistogram(), self.aoi, maxPixels=1e13)
        values = ee.Dictionary(reduction.get(self.target.rename('id')))\
                    .keys()\
                    .map(lambda x: ee.Number.parse(x))

        
        cellsOfInterest = ee.ImageCollection(values.map(lambda x: grid.eq(ee.Number(x)))).max().selfMask().rename('id')

        # Convert cells within item to vector grid cells
        gridCells = self.target.updateMask(cellsOfInterest).reduceToVectors(**{'geometry': geom, 'scale': self.scale, 'geometryType': 'polygon',\
                                'eightConnected': False, 'labelProperty': 'id', 'reducer': ee.Reducer.mode()})
        
        return gridCells

        
    def extractPoints(self, gridSize = 50000, batchSize = None, filename = 'output.csv'):
        """
        Extract covariate data at points.

       Args:
           gridSize (int): The tile size used to filter features. Runs in parralel.
           batchSize (int): The number of batches to split job into. If large gridSize results in Out of Memory errors
                specify a batchSize smaller than the number of samples. Runs in sequence.
           filename (str): The output file name.
           
       Returns:
           Data (csv) exported to download directory (dd).
        """
        
        logger = logging.getLogger(__name__)

        max_threads = self.num_threads or min(32, (os.cpu_count() or 1) + 4)
        
        #Set working directory
        if not os.path.exists(self.dd):
                os.makedirs(self.dd)
        os.chdir(self.dd)
        
        self.properties = self.covariates.bandNames().getInfo()
        self._properties = self.covariates.bandNames()
        
        # add target band name to properties
        if self.target.name == 'ee.Image':
            self.properties = ee.List(self._properties).add(self.target.bandNames()).getInfo()

        # Create grid
        grid, items = createGrid(gridSize, ee.Feature(self.aoi))
        
        self.batchSize = batchSize
        
        desc = 'Points'#str(item)
        bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
        bar = tqdm(total = grid.size().getInfo(), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

        warnings.filterwarnings('ignore', category=TqdmWarning)
        redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

        with redir_tqdm, bar:
            def downloadPoints(item):
                
                points = self.geomFeatures(grid, item)
                
                size = points.size().getInfo()
                if size>0:
                    pointsList = points.toList(size)

                    for batch in range(0, size+1, self.batchSize):
                        fc = ee.FeatureCollection(pointsList.slice(batch, batch+batchSize))

                        data = self.covariates.reduceRegions(fc, ee.Reducer.first(), self.scale)

                        output = data.map(lambda ft: ft.set('output', self._properties.map(lambda prop: ft.get(prop))))
                        result = output.aggregate_array('output').getInfo()

                        file_exists = os.path.isfile(filename)
                        # Write the results to a file.
                        csv_writer_lock = threading.Lock()
                        with csv_writer_lock:
                            with open(filename, 'a', newline='') as f:
                                writer = csv.writer(f)
                                if not file_exists:
                                    # write the header
                                    writer.writerow(self.properties)
                                    # write multiple rows
                                    writer.writerows(result)
                                    f.flush()
                                    f.close()
                                else:
                                    writer.writerows(result)
                                    f.flush()
                                    f.close()
                            
            with ThreadPoolExecutor(max_workers = max_threads) as executor:
                            # Run the tile downloads in a thread pool
                            futures = [executor.submit(downloadPoints, tile) for tile in items]
                            try:
                                for future in as_completed(futures):
                                    future.result()
                                    bar.update(1)
                                    
                            except Exception as ex:
                                logger.info('Cancelling...')
                                executor.shutdown(wait=False, cancel_futures=True)
                                raise ex        
            
    def extractPolygons(self, reduce = True, reducer = None, gridSize = 50000, batchSize = None, filename = 'output.csv'):
        """
        Extract summary statistics of covariates for regions.
        
        Args:
           reduce (bool): default True. if False, each pixel within a polygon is downloaded.
           reducer (ee.Reducer): The reducer(s) to use to summarise data. If multiple reducers need to be applied, use combined reducers.
           sparse (bool): default True. Only used when target is ee.Image. If True, samples covaraites at sparse valid target pixels (for eg., GEDI).
           If False, raster grid is converted to polygons that are used by reducRegions.
           gridSize (int): The tile size used to filter features. Runs in parralel.
           batchSize (int): The number of batches to split job into. If large gridSize results in Out of Memory errors
                specify a batchSize smaller than the number of samples. Runs in sequence.
           filename (str): The output file name.
            
        Returns:
            Data (csv) exported to download directory (dd).
            
        """
        
        logger = logging.getLogger(__name__)

        max_threads = self.num_threads or min(32, (os.cpu_count() or 1) + 4)
        
        #Set working directory
        if not os.path.exists(self.dd):
                os.makedirs(self.dd)
        os.chdir(self.dd)
        
        self._properties = self.covariates.bandNames()
        self.properties = self._properties.getInfo()
        
        # add target band name to properties
        if self.target.name == 'ee.Image':
            self.properties = self._properties.add(self.target.bandNames()).getInfo()

        # Create grid
        grid, items = createGrid(gridSize, ee.Feature(self.aoi))
        
        self.batchSize = batchSize
        
        desc = filename
        bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
        bar = tqdm(total = grid.size().getInfo(), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

        warnings.filterwarnings('ignore', category=TqdmWarning)
        redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

        with redir_tqdm, bar:
            def downloadPolygons(item):
                
                if sparse:
                    features = self.geomFeatures(grid, item)
                else:
                    features = self.geomGridCells(self.target, grid, item)
                
                size = features.size().getInfo()
                if size>0:
                    featuresList = features.toList(size)

                    for batch in range(0, size+1, self.batchSize):
                        fc = ee.FeatureCollection(featuresList.slice(batch, batch+batchSize))
                        
                        if reduce:
                            data = self.covariates.reduceRegions(fc, reducer, self.scale)
                        else:
                            data =  self.covariates.sampleRegions(fc, scale = self.scale)
                        
                        self._properties = data.first().propertyNames()
                        self.properties = self._properties.getInfo()

                        output = data.map(lambda ft: ft.set('output', self._properties.map(lambda prop: ft.get(prop))))
                        result = output.aggregate_array('output').getInfo()

                        file_exists = os.path.isfile(filename)
                        # Write the results to a file.
                        csv_writer_lock = threading.Lock()
                        with csv_writer_lock:
                            with open(filename, 'a', newline='') as f:
                                writer = csv.writer(f)
                                if not file_exists:
                                    # write the header
                                    writer.writerow(self.properties)
                                    # write multiple rows
                                    writer.writerows(result)
                                    f.flush()
                                    f.close()
                                else:
                                    writer.writerows(result)
                                    f.flush()
                                    f.close()
                            
            with ThreadPoolExecutor(max_workers = max_threads) as executor:
                            # Run the tile downloads in a thread pool
                            futures = [executor.submit(downloadPolygons, tile) for tile in items]
                            try:
                                for future in as_completed(futures):
                                    future.result()
                                    bar.update(1)
                                    
                            except Exception as ex:
                                logger.info('Cancelling...')
                                executor.shutdown(wait=False, cancel_futures=True)
                                raise ex

    def extractRasterGrid(self, reduce = True, reducer = None, gridSize = 50000, batchSize = None, filename = 'output.csv'):
        """
        Extract summary statistics of covariates for regions.
        
        Args:
           reduce (bool): default True. if False, each pixel within a polygon is downloaded.
           reducer (ee.Reducer): The reducer(s) to use to summarise data. If multiple reducers need to be applied, use combined reducers.
           gridSize (int): The tile size used to filter features. Runs in parralel.
           batchSize (int): The number of batches to split job into. If large gridSize results in Out of Memory errors
                specify a batchSize smaller than the number of samples. Runs in sequence.
           filename (str): The output file name.
            
        Returns:
            Data (csv) exported to download directory (dd).
            
        """
        
        logger = logging.getLogger(__name__)

        max_threads = self.num_threads or min(32, (os.cpu_count() or 1) + 4)
        
        #Set working directory
        if not os.path.exists(self.dd):
                os.makedirs(self.dd)
        os.chdir(self.dd)
        
        self._properties = self.covariates.bandNames()
        self.properties = self._properties.getInfo()
        
        # add target band name to properties
        if self.target.name == 'ee.Image':
            self.properties = self._properties.add(self.target.bandNames()).getInfo()

        # Create grid
        grid, items = createGrid(gridSize, ee.Feature(self.aoi))
        
        self.batchSize = batchSize
        
        desc = filename
        bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
        bar = tqdm(total = grid.size().getInfo(), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

        warnings.filterwarnings('ignore', category=TqdmWarning)
        redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

        with redir_tqdm, bar:
            def downloadPolygons(item):
                
                features = self.geomGridCells(grid, item)
                
                size = features.size().getInfo()
                if size>0:
                    featuresList = features.toList(size)

                    for batch in range(0, size+1, self.batchSize):
                        fc = ee.FeatureCollection(featuresList.slice(batch, batch+batchSize))
                        
                        if reduce:
                            data = self.covariates.reduceRegions(fc, reducer, self.scale, crs= self.crs)
                        else:
                            data =  self.covariates.sampleRegions(fc, scale = self.scale)
                        
                        self._properties = data.first().propertyNames()
                        self.properties = self._properties.getInfo()

                        output = data.map(lambda ft: ft.set('output', self._properties.map(lambda prop: ft.get(prop))))
                        result = output.aggregate_array('output').getInfo()

                        file_exists = os.path.isfile(filename)
                        # Write the results to a file.
                        csv_writer_lock = threading.Lock()
                        with csv_writer_lock:
                            with open(filename, 'a', newline='') as f:
                                writer = csv.writer(f)
                                if not file_exists:
                                    # write the header
                                    writer.writerow(self.properties)
                                    # write multiple rows
                                    writer.writerows(result)
                                    f.flush()
                                    f.close()
                                else:
                                    writer.writerows(result)
                                    f.flush()
                                    f.close()
                            
            with ThreadPoolExecutor(max_workers = max_threads) as executor:
                            # Run the tile downloads in a thread pool
                            futures = [executor.submit(downloadPolygons, tile) for tile in items]
                            try:
                                for future in as_completed(futures):
                                    future.result()
                                    bar.update(1)
                                    
                            except Exception as ex:
                                logger.info('Cancelling...')
                                executor.shutdown(wait=False, cancel_futures=True)
                                raise ex
        