'''
Created on Apr 17, 2013

@author: jwe <jweingrill@aip.de>
'''
class DataSource(object):
    def __init__(self, database=None, user=None, host=None, dictcursor=False):
        """Constructor: opens database on host using username"""
        import psycopg2
        
        try:
            self.database = psycopg2.connect(database=database, user=user, host=host) 
        finally:
            if dictcursor:
                from psycopg2.extras import DictCursor
                self.cursor = self.database.cursor(cursor_factory=DictCursor)
            else:
                self.cursor = self.database.cursor()
    
    def query(self, querystring):
        """executes a query and returns result"""
        self.cursor.execute(querystring)
        return self.cursor.fetchall()

    def execute(self, query, commit=True):
        """executes a query"""
        try:
            self.cursor.execute(query)
        finally:
            if commit:
                self.database.commit()
    
    def truncate(self, table):
        """truncates a table"""
        try:
            self.cursor.execute('TRUNCATE '+ table)
        finally:
            self.database.commit()
            
    def droptable(self, table):
        """drops a table"""
        try:
            self.cursor.execute('DROP TABLE IF EXISTS '+ table + ';')
        finally:
            self.database.commit()
            
    def dropview(self, view):
        """drops a view"""
        try:
            self.cursor.execute('DROP VIEW IF EXISTS '+ view + ';')
        finally:
            self.database.commit()
    
    def commit(self):
        """commit changes to table"""
        self.database.commit()
            
    def __exit__(self):
        """proxy function to close database"""
        self.close()
       
    def close(self):
        """closes the database"""
        if self.database:
            self.database.close()  

class Table(dict):
    '''
    class that interfaces a table on wifsip database
    '''
    keyvalue = ''
    
    def __init__(self, tablename, keyname):
        self.keyname = keyname
        self.tablename = tablename
        self.wifsip = DataSource(database='wifsip', user='sro', host='pina.aip.de')
        

    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS 
        WHERE table_name = '%s';""" % self.tablename
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        query = """SELECT * 
        FROM %s 
        WHERE %s = '%s'""" % (self.tablename, self.keyname, self.keyvalue)
        result = self.wifsip.query(query)
        values = [r for r in result[0]]
        return values

    def __setitem__(self, key, value):
        if value is None:
            query = """UPDATE %s 
            SET %s=NULL 
            WHERE %s='%s';""" % (self.tablename, key, self.keyname, self.keyvalue)
        else:
            if type(value) is str:
                value = "'%s'" % value            
            query = """UPDATE %s 
            SET %s=%s 
            WHERE %s='%s';""" % (self.tablename, key, str(value), self.keyname, self.keyvalue)
        self.wifsip.execute(query)
        
    def __getitem__(self, key):
        result = self.wifsip.query("""SELECT %s 
        FROM %s 
        WHERE %s like '%s';""" % (key, self.tablename, self.keyname, self.keyvalue))
        return result[0][0]
