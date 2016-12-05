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


    def __init__(self, datasource, tablename, condition='True'):
        '''
        Constructor
        '''
        self.datasource = datasource
        self.tablename = tablename
        
        columns = self._columns
        # crude assumption: key is the first column
        Table.__init__(self, key = columns[0], columns = columns)
        
        query = "SELECT * FROM %s WHERE %s;" % (tablename, condition)
        data = datasource.query(query)
        self.append(data)

    @property
    def _columns(self):
        query = """SELECT column_name
            FROM INFORMATION_SCHEMA.COLUMNS 
            WHERE table_name = '%s';""" % self.tablename
        result = self.datasource.query(query)
        return [r[0] for r in result]

if __name__ == '__main__':
    from datasource import DataSource
    ds = DataSource(database='stella', user='stella', host='pera.aip.de')
    
    t = DBTable(ds, 'm48phot', condition="starid LIKE '%F1%' ORDER BY ra")
    
    print t.columns
    print t['starid']
    print t['JCS-F1-003']
    
    import matplotlib.pyplot as plt
    #plt.plot(t[:,'bv'], t[:,'vmag'],',')
    import numpy as np
    x = t['ccdy']
    y = t['ra']
    p = np.polyfit(x, y, 1)
    print p
    plt.plot(x, y-np.polyval(p, x))
    plt.show()
    
    #print r.column('zero')        