import ee
import uuid
from IPython.display import display, HTML

import json
import yaml

def eeprint(obj):
    """
    A function to pretty print earth engine outputs. Works on both jupyter locally and on colab.
    You do not need to use .getInfo()

    Args:
        obj (GEE object): Takes any Google earth engine object.

    Returns:
        A collapsible JSON containing attributes/metadata.

    Examples:
        # Load NASADEM image
        image = ee.Image("NASA/NASADEM_HGT/001")
        # Print image metadata
        eeprint(image)
        # Print image bandnames
        eeprint(image.bandNames())
    """
    id = str(uuid.uuid4())
    if isinstance(obj, ee.computedobject.ComputedObject):
        obj = obj.getInfo()
    try:
        obj = yaml.load(str(obj), Loader=yaml.SafeLoader)
    except:
        pass
    if isinstance(obj, dict):
        json_str = json.dumps(obj)   
    elif isinstance(obj, str):
        return display(HTML('<div id="{}" style="height: auto; width:100%;">{}</div>'.format(id, obj)))
    else:
        json_str = str(obj)
    display(HTML('<div id="{}" style="height: auto; width:100%;"></div>'.format(id)))
    display(HTML("""
        <script src="/static/components/requirejs/require.js"></script> <!-- Needed in Colab -->
        <script>
            require(["https://rawgit.com/caldwell/renderjson/master/renderjson.js"], function() {
              renderjson.set_show_to_level(1)
              document.getElementById('%s').appendChild(renderjson(%s))
            });
        </script>
    """ % (id, json_str)))

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

def createGrid(patchSize, aoi, vect = True, list = False, crs = 'EPSG:4326'):
    """
    Generate a grid with a specified spacing.
    
    Args:
        patchSize (Int): The spatial resolution or spacing between grid cells (metres).
        aoi (ee.Geometry): The extent to return the grid.
        vect (Bool): Deafault True. if True returns a grid as an ee.Featurecollection.
                    if False, returns a grid as an ee.Image.
        list (Bool): Default False. If True returns vector grid as a ee.list of ee.Geometry.Polygons.
        crs (String): The coordinate reference system of the grid. Default is 'EPSG:4326'.
    
    Returns: 
        ee.FeatureCollection for vect=True or ee.Image for vect=False and unique grid ids (ee.List)

    Examples:
        # Create a 10km vector grid
        grid, ids = createGrid(10000, aoi)

        # Create a 10km raster grid
        grid, ids = createGrid(10000, aoi, vect = False)

        # Create a 10km vector grid as a list of ee.Geometry.Polygons
        grid, ids = createGrid(10000, aoi, list = True)

    """
    seed = 123
    
    if vect:
        grid = ee.FeatureCollection(aoi.geometry().coveringGrid(crs, patchSize))
        gridLen = ee.Number(grid.size())
        values = ee.List.sequence(0, gridLen.subtract(1))

        if list == True:
            grid = grid.toList(gridLen)
        return grid, values
    else:
        proj = ee.Projection(crs).atScale(patchSize)
        grid = ee.Image.random(seed).multiply(1e6).int().reproject(proj).rename('gid').clip(aoi)
        reduction = grid.reduceRegion(ee.Reducer.frequencyHistogram(), aoi, maxPixels=1e13, scale = patchSize)
        values = ee.Dictionary(reduction.get('gid'))\
                    .keys()\
                    
        return grid, values