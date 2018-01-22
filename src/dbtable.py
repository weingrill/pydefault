#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 16, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
from table import Table

class DBTable(Table):
    '''
    Database table Class that interfaces a table as a list of dict.
    '''


    def __init__(self, datasource, tablename, columns='*', condition='True'):
        '''
        Constructor
        '''
        self.datasource = datasource
        self.tablename = tablename
        
        
        if type(columns)==str and '*' in columns:
            # get columns by property
            columns_str = columns
        elif type(columns)==str and ',' in columns:
            # nothing to to, columns is a string
            columns_str = columns
        elif type(columns)==list:
            columns_str = ', '.join(columns) 
        
        query = "SELECT {0} FROM {1} WHERE {2};".format(columns_str, tablename, condition)
        data = datasource.query(query)
        
        columns = self._columns
        
        # crude assumption: key is the first column
        #execute super Constructor; alternative use for super(): 
        #super(Table, self).__init__(key = columns[0], columns = columns) 
        Table.__init__(self, key = columns[0], columns = columns)
        self.append(data)

    @property
    def _columns(self):
        '''property to return the list of columns in the query'''
        result = [column.name for column in self.datasource.cursor.description]
        return result


class DBArray(object):
    def __init__(self, datasource, tablename, columns='*', condition='TRUE'):
        '''
        Constructor
        '''
        import numpy as np
        
        self.datasource = datasource
        self.tablename = tablename
        self.columns = columns
        
        query = "SELECT * FROM %s WHERE %s;" % (tablename, condition)
        data = self.datasource.query(query)
        #self.append(data)
    
        columns = self.datasource.columns(self.tablename)
        data_types = self.datasource.data_types(self.tablename)
        
        for c,d in zip(columns, data_types):
            print c,d
        arraydata = []
        for star in self._stars:
            arraydata.append(tuple(star))
        self.array = np.array(arraydata, dtype = zip(columns, data_types))    

def _test_DBTable():
    from datasource import DataSource
    import matplotlib.pyplot as plt
    import numpy as np

    ds = DataSource(database='stella', user='stella', host='pera.aip.de')
    t = DBTable(ds, 'm48phot', columns=['starid','ccdy','ra'], condition="starid LIKE '%F1%' ORDER BY ra")
    
    print t.columns
    print t['starid']
    print t['JCS-F1-003']
    
    x = t['ccdy']
    y = t['ra']
    p = np.polyfit(x, y, 1)
    print p
    plt.plot(x, y-np.polyval(p, x),'o')
    plt.show()

def _test_DBArray():
    from datasource import DataSource
    import matplotlib.pyplot as plt
    ds = DataSource(database='stella', user='stella', host='pera.aip.de')
    
    a = DBArray(ds, 'candidates')
    plt.plot(a['bv'],a['period'])
    plt.show()
    
if __name__ == '__main__':
    _test_DBTable()