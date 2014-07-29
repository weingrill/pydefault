'''
Created on Jul 24, 2014

@author: jwe
'''

class Frame(dict):
    '''
    classdocs
    '''


    def __init__(self, objid):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.objid = objid
        
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        
        
#         for k in keys:
#             i = keys.index(k)
#             v = values[i]
#             print k, v
#             self[k] = v
    
    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'frames';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        query = """SELECT * from frames where objid = '%s'""" % self.objid
        result = self.wifsip.query(query)
        values = [r for r in result[0]]
        return values
        
            
    def __setitem__(self, key, value):
        if value is None:
            query = """UPDATE frames 
            SET %s=NULL 
            WHERE objid='%s';""" % (key, self.objid)
        else:            
            query = """UPDATE frames 
            SET %s=%s 
            WHERE objid='%s';""" % (key, str(value), self.objid)
        self.wifsip.execute(query)
        
    def __getitem__(self, key):
        result = self.wifsip.query("""SELECT %s 
        FROM frames 
        WHERE objid='%s';""" % (key, self.objid))
        return result[0][0]    

               
if __name__ == '__main__':
    f = Frame('20140303A-0071-0001')
    print f.keys()
    print f.values()
    print f['hjd'],f['expt']