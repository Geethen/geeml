import ee
import geemap

import logging
import multiprocessing
import requests
import shutil
from retry import retry
import os
import csv

#extract data at points
@retry(tries=10, delay=1, backoff=2)
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

#Useful during inference
@retry(tries=10, delay=1, backoff=2)
def extractAoi(index, fid):
    """
    Extract image patches of covariates within an aoi(specified in prepareForExtraction).
    
    Args:
        index (int): index returned by enumerate
        fid (int): unique id of grid. Each grid cell corresponds to a worker.
        
    Returns:
        Data (image patches as GeoTif) exported to download directory (dd) specified in prepareForExtraction function.
    """
    assert ('linux' in sys.platform), "This code runs on Linux only."
    # Select patch
    if grid.name() == 'Image':
        patchSize = grid.projection().nominalScale()
        cell = grid.updateMask(grid.eq(fid).selfMask())
        vCell = cell.addBands(cell).reduceToVectors(ee.Reducer.mean(), aoi, patchSize, maxPixels=1e13)
    elif grid.name() == 'FeatureCollection':
        vCell = ee.Feature(grid.filter(ee.Filter.eq('id', fid)).first())
    
    # Download patch as tif ()
    geemap.ee_export_image(covariates, crs= 'EPSG:4326', filename= os.path.join(dd, f"X/X_{fid}.tif"),\
                           scale= scale, region= ee.Feature(vCell.first()).geometry(), file_per_band=False)
    
    # Download patch as tif ()
    geemap.ee_export_image(target, crs= 'EPSG:4326', filename= os.path.join(dd, f"Y/Y_{fid}.tif"),\
                           scale= scale, region= ee.Feature(vCell.first()).geometry(), file_per_band=False)

#Extract data around points
@retry(tries=10, delay=1, backoff=2)
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
@retry(tries=10, delay=1, backoff=2)
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

def pExtract(extractFunction, nProcesses = 25, workers = items):
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