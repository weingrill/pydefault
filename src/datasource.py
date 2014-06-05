'''
Created on Apr 17, 2013

@author: jwe <jweingrill@aip.de>
'''
class DataSource(object):
    def __init__(self, database, user, host, dictcursor=False):
        """Constructor: opens database on host using username"""
        import psycopg2
        
        try:
            self.database = psycopg2.connect(database=database, user=user, host=host) 
        finally:
            if dictcursor:
                self.cursor = self.database.cursor(cursor_factory=psycopg2.extras.DictCursor)
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

class Table(object):
    def __init__(self, name, datasource):
        self.name = name
        self.datasource = datasource
    
    def column(self, colname):
        return self.datasource.query('SELECT '+colname+' FROM TABLE '+self.name)
        