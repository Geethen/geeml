import ee
import os

def getCountry(point, simple= True):
    """
    Returns country adminstartive boundary that point falls within
    
    Args:
        point (ee.Geometry.Point): A point within the country of interest
        simple (Bool): if True uses simplified geometry LSIB adminstrative boundaries.
                        False uses detailed LSIB adminstartive boundaries.
        
    Returns:
        Country adminstrative boundary (ee.Feature)
    """
    if simple:
        countries = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
    else:
        countries = ee.FeatureCollection('USDOS/LSIB/2017')
    
    country = countries.filterBounds(point)
    return country

def createGrid(patchSize, aoi, units='distance', scale=None, vect = True):
    """
    Generate a grid with a specified spacing.
    
    Args:
        patchSize (Int): The spatial resolution or spacing between grid cells (metres).
        aoi (geometry obj): The extent to return the grid. if None, a global grid is generated.
        vect (Bool): if True returns an ee.Featurecollection containing the grid.
                    if False, returns a grid as an ee.Image.
        units (string): Either distance (metres) or pixels. if pixels are specified, also specify scale.
        scale (int): Factor to multiply pixels. Only specify if units equals 'pixels'.
    
    Returns: 
        ee.FeatureCollection for vect=True or ee.Image for vect=False and list of
        unique ids
        
    """
    seed = 123
    if units == 'pixels':
      cellSize = patchSize*scale
    else:
      cellSize = patchSize
    proj = ee.Projection("EPSG:5070").atScale(cellSize)
    grid = ee.Image.random(seed).multiply(1e6).int().reproject(proj).rename('id').clip(aoi)
    
    if vect:
        grid = grid.reduceToVectors(**{ 'geometry': aoi.geometry().buffer(cellSize)})
        values = grid.aggregate_array('label').getInfo()
        return grid, values
    else:
        reduction = grid.reduceRegion(ee.Reducer.frequencyHistogram(), aoi, maxPixels=1e13);
        values = ee.Dictionary(reduction.get('id'))\
                    .keys()\
                    .map(lambda x: ee.Number.parse(x)).getInfo() 
        return grid, values
    return grid, values

def prepareForExtraction (covariates, grid, aoi, scale, dd, target = None, spcvGridSize = None):

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
    #Set directory
    if os.path.exists(dd):
        os.chdir(dd)
    else:
        os.makedirs(dd)
    
    if spcvGridSize is not None:
        #Generate Grid for spcv
        spcvGrid, _ = createGrid(spcvGridSize, aoi=aoi, vect=False)

        #Compile covariates
        if covariates.name() in ['Image','ImageCollection']:
            finalCovariates = spcvGrid.addBands(ee.ImageCollection(covariates).toBands()).float()
        elif covariates.name() == 'List':
            for covariate in covariates:
                spcvGrid = spcvGrid.addBands(ee.Image(covariate))
            finalCovariates = spcvGrid.float()
        else:
            raise TypeError("Only takes single image (ee.Image), list of image id's or ee.Imagecollection as covariates")
    else:
        #Compile covariates
        if covariates.name() in ['Image','ImageCollection']:
            finalCovariates = ee.ImageCollection(covariates).toBands().float()
        elif covariates.name() == 'List':
            for covariate in covariates:
                finalCovariates = ee.Image(covariate).float()
        else:
            raise TypeError("Only takes single image (ee.Image), list of image id's or ee.Imagecollection as covariates")

    #Accomodate target
    if target is not None:
        target = target.rename('target').float()
    
    return covariates, target, grid, aoi, scale