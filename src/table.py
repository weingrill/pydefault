#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 12, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Table(object):
    '''
    Table class to access records in a table by their key or their column
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
        Record['columnname'] returns the column with the columnname
        Record[1] returns the second record
        Record[1:3] returns record 1 and 2 (the second and the third)
        """
        from numpy import array
        if type(keyvalue) == str:
            # try to find the key with keyvalue
            for row in self.table:
                if row[self.key] == keyvalue:
                    return row
            if keyvalue in self.columns:
                return array([row[keyvalue] for row in self.table])
            else:
                raise(ValueError)
                    
        elif type(keyvalue) == int or type(keyvalue) == slice:
            return self.table[keyvalue]
        
        elif type(keyvalue) == tuple:
            #print keyvalue[0]
            return array([row[keyvalue[1]] for row in self.table[keyvalue[0]]])
        else:
            print 'unknown type ', type(keyvalue)
            
    def __len__(self):
        """
        returns the number of records
        """
        
        return len(self.table)
    
    def keys(self):
        """
        returns the name of the columns
        """

        return self.columns
    
    def __call__(self, **kwargs):
        """
        returns a filtered subtable
        """

        key, value = kwargs.items()[0]
        childrecord = Table(key = self.key, columns = self.columns)
        for row in self.table: 
            if row[key] == value:
                childrecord.append(row)
        return childrecord

    def __str__(self):
        """
        returns a tab separated output of the table
        """
        
        s = ''
        s += '\t'.join(self.table[0].keys()) + '\n'
        for row in self.table:
            s += '\t'.join([str(row[c]) for c in self.table[0].keys()]) + '\n'
        return s

        
if __name__ == '__main__':
    t = Table(key = 'name', columns=['name','email','salary','location'])
    t.append({'name':'Joerg', 'email':'joerg@weingrill.net', 'salary': 2500.0, 'location': 'Potsdam'})
    t.append({'name':'Katja', 'email':'kjanssen@aip.de', 'salary': 2800.0, 'location': 'Potsdam'})
    t.append({'name':'Trey', 'email':'cmack@aip.de', 'salary': 2300.0, 'location': 'Berlin'})
    print t['Joerg']
    print t['Katja']
    print t[:,'salary']
    print t['salary']
    print 'Who is in Potsdam?', t(location = 'Potsdam')
    #print r1['Susy']
    
    print t[1:3]
    print t[1:3,'salary']
    print t
