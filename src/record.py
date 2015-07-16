#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 12, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Record(object):
    '''
    Record class to access records in a table by their key or their column
    '''

    def __init__(self, key=None, columns=None):
        '''
        key needs to be a string
        columns need to be a list of str
        '''
        if type(key)<>str:
            raise(TypeError)
        if type(columns)<>list:
            raise(TypeError)
        
        self.key = key
        self.columns = columns
        self.table = []
    
    def append(self, newrecord, columns=None):
        """
        appends a new record
        """
        
        if type(newrecord)==dict:
            self.table.append(newrecord)
            
        elif type(newrecord)==list:
            if columns is None:
                columns = self.columns
            for record in newrecord:
                recdict = {}
                for i in range(len(record)):
                    recdict[columns[i]] = record[i]
                self.table.append(recdict)
        else:
            raise(TypeError)
    
    def __getitem__(self, keyvalue):
        """
        Record['keyvalue'] returns the record where key=keyvalue
        Record[:, 'columnname'] returns the column with the columnname
        Record[1] returns the second record
        Record[1:3] returns record 1 and 2 (the second and the third)
          
        """
        if type(keyvalue) == str:
            # try to find the key with keyvalue
            for row in self.table:
                if row[self.key] == keyvalue:
                    return row
            if keyvalue in self.columns:
                return [row[keyvalue] for row in self.table]
            else:
                raise(ValueError)
                    
        elif type(keyvalue) == int or type(keyvalue) == slice:
            return self.table[keyvalue]
        
        elif type(keyvalue) == tuple:
            #print keyvalue[0]
            return [row[keyvalue[1]] for row in self.table[keyvalue[0]]]
        else:
            print 'unknown type ', type(keyvalue)
            
    def __len__(self):
        return len(self.table)
    
    def keys(self):
        return self.columns
    
    def __call__(self, **kwargs):
        key, value = kwargs.items()[0]
        childrecord = Record(key = self.key, columns = self.columns)
        for row in self.table: 
            if row[key] == value:
                childrecord.append(row)
        return childrecord

    def __str__(self):
        s = ''
        s += '\t'.join(self.table[0].keys()) + '\n'
        for row in self.table:
            s += '\t'.join([str(row[c]) for c in self.table[0].keys()]) + '\n'
        return s
        
if __name__ == '__main__':
    r1 = Record(key = 'name', columns=['name','email','salary','location'])
    r1.append({'name':'Joerg', 'email':'joerg@weingrill.net', 'salary': 2500.0, 'location': 'Potsdam'})
    r1.append({'name':'Katja', 'email':'kjanssen@aip.de', 'salary': 2800.0, 'location': 'Potsdam'})
    r1.append({'name':'Trey', 'email':'cmack@aip.de', 'salary': 2300.0, 'location': 'Berlin'})
    print r1['Joerg']
    print r1['Katja']
    print r1[:,'salary']
    print r1['salary']
    print 'Who is in Potsdam?', r1(location = 'Potsdam')
    #print r1['Susy']
    
    print r1[1:3]
    print r1[1:3,'salary']
    print r1
    exit()
    r = Record(key = 'objid', columns=['objid','path','filename', 'fwhm', 'zero'])
    import psycopg2

    
    query = """(SELECT objid, path, filename, fwhm, zero 
        FROM science
        WHERE wcs IS TRUE AND phot IS FALSE)
        UNION
        (SELECT science.objid, path, filename, fwhm, zero 
        FROM science, frames
        WHERE science.objid=frames.objid 
        ORDER BY objid DESC) 
        LIMIT 200;"""
     
    wifsip = psycopg2.connect(database='stella', user='stella', host='pera.aip.de')
    cur = wifsip.cursor()
    cur.execute(query)
    datas = cur.fetchall()
    r.append(datas)
    print r[:,'fwhm']
    #print r.column('fwhm')
    import matplotlib.pyplot as plt
    plt.plot(r[:,'fwhm'],r[:,'zero'],'o')
    plt.show()
    
    #print r.column('zero')