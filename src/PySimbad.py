import urllib2
import urllib, sys

class SimbadError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def simbad(psr):
    """Search simbad database for a pulsar's information, especially the coordinates."""
    data = {}
    data['output.format'] = 'ASCII'
    data['Ident'] = psr
    url_values = urllib.urlencode(data)
    url = 'http://simbad.u-strasbg.fr/simbad/sim-id'
    full_url = url + '?' + url_values
    data = urllib2.urlopen(full_url)
    results = data.read()
    info = {}
    failed = False
    for lines in results.split('\n')[7:]:
        try:
            key,value = lines.split(':')
        except ValueError:
            try:
                value+=lines
            except UnboundLocalError:
                key='Search Result'
                value = lines
                failed = True
        except:
            pass
        info.update({key:value})
    for key in info.keys():
        value = info[key]
        value.replace('\n','')
        if not value.find('~ ~') == -1:value=None
        info[key] = value
    if failed:
        for key in info.keys():
            print '%s:\n%s\n' % (key, info[key])
        raise SimbadError('Failed to identify the object %s' % psr)
    else:
        return info


def SimbadCoord(psr):
    info = simbad(psr)
    #print info
    try: 
        coords = info['Coordinates(ICRS,ep=2000,eq=2000)']
    except KeyError:
        try:
            coords = info['Coordinates(FK5,ep=J2000,eq=2000)']
        except KeyError:
            coords = info['Coordinates(FK4,ep=2000,eq=2000)']
    return coords
