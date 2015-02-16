import urllib2
import urllib

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

def simcoo(coordinates, dec=None):
    """Search simbad database for a identifier based on coordinates."""
    data = {}
    data['output.format'] = 'ASCII'
    if dec is None:
        data['Coord'] = str(coordinates)
    else:
        data['Coord'] = '%f %f' % (coordinates, dec)
    data['Radius'] = 3
    data['Radius.unit'] = 'arcsec'
    url_values = urllib.urlencode(data)
    url = 'http://simbad.u-strasbg.fr/simbad/sim-coo'
    full_url = url + '?' + url_values
    data = urllib2.urlopen(full_url)
    results = data.read()
    for lines in results.split('\n'):
        if lines.startswith('Object'):
            ls = lines.split('  ---  ')
            return ls[0][7:]
    print results


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

if __name__ == '__main__':
    #print simcoo('123.332697 -5.558109')
    #print simcoo('123.570997 -5.900061')
    #print simcoo('123.448826 -5.623558')
    #print simcoo('123.325151 -5.640358')
    #print simcoo('123.325125 -5.64067')
    m48 = simbad('M48')
    for k in m48.keys():
        print k,'\t',m48[k]
    exit()
    print simcoo('123.564614 -5.720953')
    print simcoo(123.369415, -5.804038)
    print simcoo(123.541259,-5.883464)
    print simcoo(123.678062,-5.710619)