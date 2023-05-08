import logging
import os

import ee
import geedim as gd
from .utils import createGrid, eeprint

import threading
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import pandas as pd

from tqdm import TqdmWarning
from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

class extractor:

    def __init__(self, covariates, aoi, scale, dd, crs = 'EPSG:4326', target = None, spcvGridSize = None, num_threads = 5):
        """
        Prepares explanatoy variables and response variables for data exatraction .
        A grid for spatial cross validation (spcv) is also added as a band if spcvGridSize is specified.
        The directory where data will be saved is set.
        
        Args:
            covariates (ee.List, ee.Image, ee.ImageCollection): if ee.List, list should contain images (ee.Image).
                If imagecollection, each image should be a single covariate.

            aoi (ee.Feature): The area of interest that designates the download extent.

            scale (int): The output spatial resolution for exported data.

            dd (string): The destination to save exported/downloaded data.

            crs (str): The coordinate reference system for exported data.

            target (ee.FeatureCollection, ee.Image): The response variable. Include as a band (ee.Image) in ee.Image
                if it also needs to be extracted when using extractByGrid, extractPoints or extractFeatures. If
                ee.FeatureCollection, covariates will be sampled at all target pixels.

            spcvGridSize (int): The size of blocks to use for spatial cross-validation (spcv). All samples within a block are 
                assigned the same id.
                
            num_threads (int): The number of threads to use when extracting data.
            
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
            bandnames = self.covariates.bandNames().getInfo()
            if self.covariates.name() == 'ImageCollection':
                finalCovariates = self.spcvGrid.addBands(ee.ImageCollection(self.covariates).toBands()).float().rename(['id'] + bandnames)
            elif self.covariates.name() == 'List':
                finalCovariates = self.spcvGrid.addBands(ee.ImageCollection.fromImages(self.covariates).toBands().float())
            elif self.covariates.name() == 'Image':
                finalCovariates = self.spcvGrid.addBands(self.covariates.float()).rename(['id'] + bandnames)
            else:
                raise TypeError("Only takes a single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
        else:
            #Compile covariates
            if self.covariates.name() == 'ImageCollection':
                finalCovariates = ee.ImageCollection(self.covariates).toBands().float()
            elif self.covariates.name() == 'List':
                finalCovariates = ee.ImageCollection.fromImages(self.covariates).toBands().float()
            elif self.covariates.name() == 'Image':
                finalCovariates = self.covariates
            else:
                raise TypeError("Only takes a single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
        
        # single or multi band image
        self.covariates = finalCovariates

        #Format target
        if self.target is not None and self.target.name() == 'ee.Image':
            self.target = self.target.rename('target').float()


    def extractAoi(self):

        """
        Download an ee.Image (single or multiband image). If ee.Imagecollection consider using the mosaic() 
        or the toBands() functions to convert to an ee.Image.
                    
        Returns:
            Image (GeoTiff) exported to download directory (dd). Within the specified download directory (dd)
            Explanatory variables/features are saved a folder called 'X' and response variables/targets are saved
            to a folder called 'Y'.

        Example:
            from geeml.extract import extractor
            import geedim as gd
            import ee

            ee.Authenticate()
            ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com') # initialise earth engine

            # Load data
            covariates = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate('2018-01-01', '2018-12-31').mosaic()
            target = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('occurrence')
            aoi = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017').filter(ee.Filter.eq('country_na', 'Kenya'))
            
            # Create extractor object
            ex = extractor(covariates, aoi, 30, 'C:/Users/username/Desktop/kenya', target = target)

            # Extract data
            ex.extractAoi()
            
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
            

    def _geomFeatures(self, grid, item):
        """
        If target is an ee.Image, for a grid (ee.FeatureCollection), filter the grid to a grid cell and return points for
        all the pixels within the selected grid cell. If the target is an ee.FeatureCollection, returns the features that
        overlap the selected grid cell.
        
        Args:
            grid (ee.FeatureCollection): The grid to parralelise over
            item (int): corresponds to an id within the grid (from createGrid function).
        
        Returns:
            features (point or polygons) that overlap the specified grid item

        Example:
            from geeml.extract import extractor
            import ee

            ee.Authenticate()
            ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com') # initialise earth engine

            # Create grid
            grid, ids = createGrid(50000, aoi = aoi, vect = False, list = True)

            # Filter grid to a single grid cell
            item = ids[0]
            features = _geomFeatures(grid, item)
                    
        """
        geom = ee.Feature(grid.get(item)).geometry()

        # Sample all pixels at points
        if self.target.name() == 'Image':
            projection = self.target.projection()
            # Extract data
            features = self.target.sample(**{'dropNulls': True, 'factor': None,\
                                'numPixels': None,  'region': geom,\
                            'scale': self.scale,'projection': projection,'geometries':True})
        else:
            features = self.target.filterBounds(geom)
        
        return features

    def _geomGridCells(self, grid, item):
        """
        If target is an ee.Image, for a grid (ee.FeatureCollection), filter the grid to a grid cell and return points for
        all the pixels within the selected grid cell. If the target is an ee.FeatureCollection, returns the features that
        overlap the selected grid cell. Raster equivalent of _geomFeatures. It is recommended to use _geomFeatures as it
        is more efficient.
        
        Args:
            grid (ee.Image): The grid to parralelise over
            item (int): corresponds to an id within the grid (from createGrid function).
        
        Returns:
            features (point or polygons) that overlap the specified grid item

        Example:
            from geeml.extract import extractor
            import ee

            ee.Authenticate()
            ee.Initialize() # initialise earth engine

            # Create grid
            grid, ids = createGrid(50000, aoi = aoi, vect = False, list = True)

            # Filter grid to a single grid cell
            ids = ids.getInfo()
            item = ids[0]

            features = _geomFeatures(grid, item)
                            OR
            features = _geomFeatures(grid, 0)
        
        """
        geom = ee.Feature(grid.filter(ee.Filter.eq('label', item)).first()).geometry()

        reduction = self.target.clip(geom).rename('id').reduceRegion(ee.Reducer.frequencyHistogram(), self.aoi, maxPixels=1e13)
        values = ee.Dictionary(reduction.get('id'))\
                    .keys()\
                    .map(lambda x: ee.Number.parse(x))

        
        cellsOfInterest = ee.ImageCollection(values.map(lambda x: self.target.eq(ee.Number(x)))).max().selfMask().rename('id')

        # Convert cells within item to vector grid cells
        gridCellsImg = self.target.updateMask(cellsOfInterest)
        gridCells = gridCellsImg.addBands(gridCellsImg).reduceToVectors(**{'geometry': self.aoi, 'scale': self.scale, 'geometryType': 'polygon',\
                                'eightConnected': False, 'labelProperty': 'id', 'reducer': ee.Reducer.first()})
        
        return gridCells

        
    def extractPoints(self, gridSize = 50000, batchSize = 2500, filename = 'output.csv'):
        """
        Extract covariate data at points.

        Args:
           gridSize (int): The grid size used to filter features. Runs in parralel.
           batchSize (int): The number of batches to split job into. If large gridSize results in Out of Memory errors
                either specify a smaller batchSize or a smaller gridSize.
           filename (str): The output csv file name.
           
        Returns:
           Data (csv) exported to download directory (dd).

        Example:
            from geeml.extract import extractor
            import ee

            ee.Authenticate()
            ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com') # initialise earth engine

            # Create extractor object
            covariates = ee.ImageCollection('COPERNICUS/S2_SR').filterDate('2019-01-01', '2019-12-31')

            # Define 100 random points
            points = ee.FeatureCollection.randomPoints(aoi, 100, 123)

            aoi = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
            aoi = aoi.filter(ee.Filter.eq('country_na', 'South Africa'))
            ex = extractor(covariates, aoi, 30, 'C:/Users/Desktop/kenya', target = points)

            # Extract data
            ex.extractPoints()

        """
        
        logger = logging.getLogger(__name__)

        #Set download directory if it exists. If not creates directory and sets as download directory
        if not os.path.exists(self.dd):
                os.makedirs(self.dd)
        os.chdir(self.dd)

        max_threads = self.num_threads or min(32, (os.cpu_count() or 1) + 4)

        self._properties = self.covariates.bandNames()
        self.properties = self._properties.getInfo()
        self.batchSize = batchSize

        # add target band name to properties
        if self.target.name == 'ee.Image':
            self.properties = self._properties.add(self.target.bandNames()).getInfo()

        # Create grid
        grid, items = createGrid(gridSize, ee.Feature(self.aoi), crs = self.crs, list = True)
        items = items.getInfo()
        
        desc = filename
        bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
        bar = tqdm(total = len(items), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

        warnings.filterwarnings('ignore', category=TqdmWarning)
        redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

        with redir_tqdm, bar:
            def downloadPoints(item):
                
                points = self._geomFeatures(grid, item)
                
                size = points.size().getInfo()
                if size>0:
                    pointsList = points.toList(size)

                    for batch in range(0, size+1, self.batchSize):
                        fc = ee.FeatureCollection(pointsList.slice(batch, batch+self.batchSize))

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
            
    def extractPolygons(self, reduce = True, reducer = None, gridSize = 50000, batchSize = 1000, filename = 'output.csv'):
        """
        Extract summary statistics of covariates for irregular polygons.
        
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
        grid, items = createGrid(gridSize, ee.Feature(self.aoi), crs = self.crs, list = True)
        items = items.getInfo()
        
        self.batchSize = batchSize
        
        desc = filename
        bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
        bar = tqdm(total = len(items), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

        warnings.filterwarnings('ignore', category=TqdmWarning)
        redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

        with redir_tqdm, bar:
            def downloadPolygons(item):
                
                features = self._geomFeatures(grid, item)
                
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

    # def byBands(self, reduce = True, reducer = None, batchSize = None, filename = 'output.csv'):
    #     """
    #     Extract or reduce pixel values by points/polygons/grid band-wise. Bands can be added to a pre-existing csv file.
    #     If download fails, the function will attempt to download the failed bands again with different settings (batchSize) 
    #     to ensure download completes.
        
    #     Args:
    #         reduce (bool): default True. if False, each pixel within a polygon is downloaded.
    #         reducer (ee.Reducer): The reducer(s) to use to summarise data. If multiple reducers need to be applied, use combined reducers.
    #         batchSize (int): The number of batches to split job into. If large gridSize results in Out of Memory errors
    #             specify a batchSize smaller than the number of samples. Runs in sequence.
    #         filename (str): The output file name.
            
    #     Returns:
    #         Data (csv) exported to download directory (dd).
            
    #     """
    #     logger = logging.getLogger(__name__)

    #     max_threads = self.num_threads or min(32, (os.cpu_count() or 1) + 4)
        
    #     #Set working directory
    #     if not os.path.exists(self.dd):
    #             os.makedirs(self.dd)
    #     os.chdir(self.dd)
        
    #     self._properties = self.covariates.bandNames()
    #     self.properties = self._properties.getInfo()
        
    #     # add target band name to properties
    #     if self.target.name == 'ee.Image':
    #         self.properties = self._properties.add(self.target.bandNames()).getInfo()
        
    #     size = self.target.size().getInfo()
        
    #     if batchSize is None:
    #         self.batchSize = size+1
    #     else:
    #         self.batchSize = batchSize

    #     featuresList = self.target.toList(size)
        
    #     desc = filename
    #     bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
    #     bar = tqdm(total = len(self.properties), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

    #     warnings.filterwarnings('ignore', category=TqdmWarning)
    #     redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

    #     with redir_tqdm, bar:
    #         def downloadPolygons(band):
    #             file_exists = os.path.isfile(filename)
    #             if file_exists:
    #                 # check how many variables have been downloaded
    #                 df = pd.read_csv(filename)
    #                 header = df.columns
    #                 # remove bands from properties
    #                 self.properties = [x for x in self.properties if x not in header]
    #             else:
    #                 df = pd.DataFrame()

    #             for band in self.properties:
    #                 result = pd.DataFrame()
    #                 for batch in range(0, size+1, self.batchSize):
    #                     fc = ee.FeatureCollection(featuresList.slice(batch, batch + self.batchSize))
                                
    #                     if reduce:
    #                         data = self.covariates.select(band).reduceRegions(fc, reducer, self.scale).map(lambda ft: ft.set('id', ft.id()))
    #                     else:
    #                         data =  self.covariates.select(band).sampleRegions(fc, scale = self.scale)

    #                     propNames = data.first().propertyNames().slice(0,2)
    #                     output = data.map(lambda ft: ft.set('output', propNames.map(lambda prop: ft.get(prop))))
    #                     result = output.aggregate_array('output').getInfo()
    #                     result = pd.DataFrame(result, columns = ['id', band])

    #                     if df.shape[0] == 0:
    #                         df = result
    #                     else:
    #                         df = df.merge(result, on ='id', how ='outer')
    #                     writer_lock = threading.Lock()
    #                     with writer_lock:
    #                         df.to_csv(filename)

    #         with ThreadPoolExecutor(max_workers = max_threads) as executor:
    #                         # Run the tile downloads in a thread pool
    #                         futures = [executor.submit(downloadPolygons, band) for band in self.properties]
    #                         try:
    #                             for future in as_completed(futures):
    #                                 future.result()
    #                                 bar.update(1)
                                    
    #                         except Exception as ex:
    #                             logger.info('Cancelling...')
    #                             executor.shutdown(wait=False, cancel_futures=True)
    #                             raise ex


    def extractByGrid(self, reduce = True, reducer = None, gridSize = 50000, batchSize = 250, filename = 'output.csv'):
        """
        Extract pixel values (if reducer = True) or summary statistics (reducer = True) of covariates by target grid.
        
        Args:
           reduce (bool): default True. if False, each pixels' values within a grid cell is extracted.
           reducer (ee.Reducer): The reducer(s) to use to summarise data. If multiple reducers need to be applied,
                                 use combined reducers (recommended for multi-reducers).
           gridSize (int): The tile size used to filter features. Runs in parralel.
           batchSize (int): The number of batches to split a job into. If large gridSize results in Out of Memory errors
                specify a  smaller batchSize.
           filename (str): The output file name.
            
        Returns:
            Data (csv) exported to download directory (dd).

        Example:
            from geeml.extract import extractor
            from geeml.utils import createGrid
            import ee

            ee.Authenticate()
            ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')

            # Define an area of interest of South Africa
            aoi = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017').filter(ee.Filter.eq('country_na', 'South Africa'))
            # Define SRTM to extract
            srtm = ee.Image('USGS/SRTMGL1_003')
            # Create a Grid of the aoi
            grid, ids = createGrid(10000, ee.Feature(aoi), vect = False)

            # Define a extractor object
            ex = extractor(srtm, aoi, 1000, 'C:/Users/username/Desktop/kenya', target = table)

            # Extract each SRTM value in each 10x10 km grid cell at a 1000m resolution
            ex.extractByGrid(reduce = False, gridSize = 10000, filename = 'srtm.csv')

            # Extract the mean SRTM value per grid cell
            ex.extractByGrid(reduce = True, reducer = ee.Reducer.mean(), scale = 30, gridSize = 10000, batchSize = 2500)

            # Extract the mean and standard deviation SRTM value per grid cell at a 1km resolution for South Africa
            combinedReducer = ee.Reducer.mean().combine(reducer2 = ee.Reducer.stdDev(), sharedInputs = True)
            ex.extractByGrid(reduce = True, reducer = combinedReducer, scale = 1000, gridSize = 10000, batchSize = 2500)

            
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
        if self.target is not None and self.target.name == 'ee.Image':
            self.properties = self._properties.add(self.target.bandNames()).getInfo()

        # Create grid
        grid, items = createGrid(gridSize, ee.Feature(self.aoi), crs = self.crs, list=True)
        items = items.getInfo()
        
        self.batchSize = batchSize
        
        desc = filename
        bar_format = ('{desc}: |{bar}| [{percentage:5.1f}%] in {elapsed:>5s} (eta: {remaining:>5s})')
        bar = tqdm(total = len(items), desc=desc, bar_format=bar_format, dynamic_ncols=True, unit_scale=True, unit='B')

        warnings.filterwarnings('ignore', category=TqdmWarning)
        redir_tqdm = logging_redirect_tqdm([logging.getLogger(__package__)])  # redirect logging through tqdm

        with redir_tqdm, bar:
            def downloadPolygons(item):
                features = self._geomFeatures(grid, item)
                
                size = features.size().getInfo()
                if size>0:
                    featuresList = features.toList(size)

                    for batch in range(0, size+1, self.batchSize):
                        fc = ee.FeatureCollection(featuresList.slice(batch, batch+self.batchSize))
                        
                        if reduce:
                            data = self.covariates.reduceRegions(fc, reducer, self.scale, crs= self.crs)
                        else:
                            data =  self.covariates.sampleRegions(fc, scale = self.scale)

                        if data.size().getInfo()>0:
                            self._properties = data.first().propertyNames()
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
                            futures = [executor.submit(downloadPolygons, item) for item in items]
                            try:
                                for future in as_completed(futures):
                                    future.result()
                                    bar.update(1)
                                    
                            except Exception as ex:
                                logger.info('Cancelling...')
                                executor.shutdown(wait=False, cancel_futures=True)
                                raise ex        