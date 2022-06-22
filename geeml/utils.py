import ee
import uuid
from IPython.core.display import display, HTML

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