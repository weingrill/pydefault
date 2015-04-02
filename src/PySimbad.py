import urllib2
import urllib

class SimbadError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SimbadObject(dict):
    def __init__(self, identifier):
        data = {}
        data['output.format'] = 'ASCII'
        data['Ident'] = identifier
        url_values = urllib.urlencode(data)
        url = 'http://simbad.u-strasbg.fr/simbad/sim-id'
        full_url = url + '?' + url_values
        data = urllib2.urlopen(full_url)
        results = data.read()
        """['Proper motions', 
        'Flux V ', 
        'Redshift', 
        'Notes (0) ', 
        'Identifiers (7)', 
        'Spectral type', 
        'Coordinates(Gal,ep=J2000,eq=2000)', 
        'Coordinates(FK4,ep=B1950,eq=1950)', 
        'Morphological type', 
        'cz', 
        'Bibcodes  1850-2015 () (1613)', 
        'Parallax', 
        'Angular size', 
        'Flux B ', 
        'Coordinates(FK5,ep=J2000,eq=2000)', 
        'nullCoordinates(ICRS,ep=J2000,eq=2000)', 
        'Radial Velocity']"""
        
        def nfloat(s):
            try:
                result = float(s)
            except ValueError:
                return None
            return result

        def nint(s):
            try:
                result = int(s)
            except ValueError:
                return None
            return result

        def nstr(s):
            if s == '~': return None
            else: return s
            
        def flux(splitvalues):
            return [nfloat(splitvalues[0]),
                         nfloat(splitvalues[1].strip('[]')),
                         splitvalues[2],
                         splitvalues[3]]
        for lines in results.split('\n')[7:]:
            try:
                key,value = lines.split(':')
            except ValueError:
                try:
                    value+=lines
                except UnboundLocalError:
                    key='Search Result'
                    value = lines
            except:
                pass
            self.update({key:value})
            
        for key in self.keys():
            sv = self[key].split()    
            if key == 'Proper motions':
                value = [nfloat(sv[0]), 
                         nfloat(sv[1]),
                         nfloat(sv[2].strip('[]')),
                         nfloat(sv[3]),
                         nfloat(sv[4].strip('[]')),
                         nstr(sv[5]),
                         nstr(sv[6])]
            elif key == 'Flux U ':
                value = flux(sv)
            elif key == 'Flux V ':
                value = flux(sv)
            elif key == 'Flux B ':
                value = flux(sv)
            elif key == 'Flux H ':
                value = flux(sv)
            elif key == 'Flux J ':
                value = flux(sv)
            elif key == 'Flux K ':
                value = flux(sv)
            elif key == 'Radial Velocity':
                value = [nfloat(sv[0]),
                         nfloat(sv[1].strip('[]'))]
            elif key == 'cz':
                value = [nfloat(sv[0]),
                         nfloat(sv[1].strip('[]')),
                         nstr(sv[2]),
                         nstr(sv[3])]
            elif key == 'Redshift':
                value = [nfloat(sv[0]),
                         nfloat(sv[1].strip('[]')),
                         nstr(sv[2]),
                         nstr(sv[3])]
            elif key == 'Parallax':
                value = [nfloat(sv[0]),
                         nfloat(sv[1].strip('[]')),
                         nstr(sv[2]),
                         nstr(sv[3])]
            elif key == 'Angular size':
                value = [nfloat(sv[0]),
                         nfloat(sv[1]),
                         nint(sv[2]),
                         nstr(sv[3].strip('()')),
                         nstr(sv[4]),
                         nstr(sv[5]),
                         nstr(sv[6])]
                
            elif key == 'Spectral type':
                value = [nstr(s) for s in sv]
            elif key == 'Morphological type':
                value = [nstr(s) for s in sv]
                
            elif key[:11] == 'Identifiers':
                value = self.pop(key)
                sv = value.split('  ')
                sv = [s.strip() for s in sv if len(s)>=1]
                key = 'Identifiers'
                value = sv
            elif key[:8] == 'Bibcodes':
                value = self.pop(key)
                key = 'Bibcodes'
                value = value.split()
            elif key[:5] == 'Notes':
                value = self.pop(key).strip()
                key = 'Notes'
            elif key[:15] == 'Coordinates(Gal':
                value = self.pop(key)
                key = 'Galactic Coordinates'
                value = [float(v) for v in value.split()]
            elif key[:15] == 'Coordinates(FK4':
                value = self.pop(key).strip()
                key = 'FK4 Coordinates'
            elif key[:15] == 'Coordinates(FK5':
                value = self.pop(key).strip()
                key = 'FK5 Coordinates'
            elif key == 'nullCoordinates(ICRS,ep=J2000,eq=2000)':
                value = self.pop(key).strip()
                key = 'ICRS Coordinates'    
                coords, _, rest = value.partition('(')
                src, _, rest = rest.partition(')')
                ident, _, rest = rest.partition('[')
                err, _, ref = rest.partition(']')
                err = err.split()
                errarr = [nfloat(err[0]),nfloat(err[1]), nint(err[2])]
                value = [coords.strip(), 
                         nstr(src.strip()), 
                         nstr(ident.strip()), 
                         errarr, 
                         ref.strip()]
            else:
                print key,':'
                print sv
            self.update({key:value})
            key, value = None, None

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
    #s = SimbadObject('HD192043')
    s = SimbadObject('NGC6882')
    print s['ICRS Coordinates']
    s1 = SimbadObject('HD192043')
    print s1['ICRS Coordinates']
    
    #print simcoo('123.332697 -5.558109')
    m48 = simbad('M48')
    for k in m48.keys():
        print k,'\t',m48[k]
    exit()