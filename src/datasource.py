'''
Created on Apr 17, 2013

@author: jwe <jweingrill@aip.de>
'''
def format_block(block,nlspaces=0):
    '''Format the given block of text, trimming leading/trailing
    empty lines and any leading whitespace that is common to all lines.
    The purpose is to let us list a code block as a multiline,
    triple-quoted Python string, taking care of
    indentation concerns.'''

    import re

    # separate block into lines
    lines = str(block).split('\n')

    # remove leading/trailing empty lines
    while lines and not lines[0]:  del lines[0]
    while lines and not lines[-1]: del lines[-1]

    # look at first line to see how much indentation to trim
    ws = re.match(r'\s*',lines[0]).group(0)
    if ws:
        lines = map( lambda x: x.replace(ws,'',1), lines )

    # remove leading/trailing blank lines (after leading ws removal)
    # we do this again in case there were pure-whitespace lines
    while lines and not lines[0]:  del lines[0]
    while lines and not lines[-1]: del lines[-1]

    # account for user-specified leading spaces
    flines = ['%s%s' % (' '*nlspaces,line) for line in lines]

    return '\n'.join(flines)+'\n'

class DataSource(object):
    def __init__(self, database=None, user=None, host=None, dictcursor=False):
        """Constructor: opens database on host using username"""
        import psycopg2
        
        try:
            self.database = psycopg2.connect(database=database, user=user, host=host) 
        except psycopg2.OperationalError:
            raise(psycopg2.OperationalError)
            
        if dictcursor:
            from psycopg2.extras import DictCursor
            self.cursor = self.database.cursor(cursor_factory=DictCursor)
        else:
            self.cursor = self.database.cursor()
    
    def query(self, querystring, queryparameters=None):
        """executes a query and returns result"""
        if queryparameters is None:
            self.cursor.execute(format_block(querystring))
        else:
            self.cursor.execute(format_block(querystring), queryparameters)
        return self.cursor.fetchall()

    def execute(self, query, commit=True):
        """executes a query"""
        try:
            self.cursor.execute(format_block(query))
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
    
    def columns(self, tablename):
        """returns the columns of the table"""
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS 
        WHERE table_name = '%s';""" % tablename
        result = self.query(query)
        column_names = [r[0] for r in result]
        return column_names
    
    def data_types(self, tablename):
        """returns the datatypes of the table"""
        from numpy import bool_, int32, float16, float32
        from astropy.coordinates import SkyCoord  # @UnresolvedImport
        
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS 
        WHERE table_name = '%s';""" % tablename
        result = self.query(query)
        keys = [r[1] for r in result]
        lens = [r[2] for r in result]
        translate = {
                     'character varying': 'S', 
                     'real':    float16, 
                     'integer': int32, 
                     'double precision': float32,
                     'point': SkyCoord,
                     'boolean': bool_
                     }
        dtypes = [translate[k] for k in keys]
        for i,dtype in enumerate(dtypes):
            if dtype == 'S':
                dtypes[i] = 'S'+str(lens[i])
            
        return dtypes
        
    def character_maximum_length(self, tablename):
        """returns the datatypes of the table"""
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS 
        WHERE table_name = '%s';""" % tablename
        result = self.query(query)
        keys = [r[2] for r in result]
        return keys

    def description(self):
        return self.database.description()

            
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
