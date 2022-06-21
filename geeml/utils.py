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
    

