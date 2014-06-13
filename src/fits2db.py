'''
Created on Apr 26, 2013

@author: jwe
'''
def fits2db(filename, relname, database='corot', user='sro', host='pina'):
    from datasource import DataSource
    import pyfits
    import StringIO
    
    ds = DataSource(database, user, host)
    hdulist = pyfits.open(filename)
        #hdr = hdulist[0].header
    tab = hdulist[1].data
    hdulist.close()
    # drop table if exists
    ds.droptable(relname)
    
    # create table
    dtypes = {'>f8': 'double precision',
              '>f4': 'real',
              '>i2': 'smallint',
              '>i4': 'integer'}
    q = 'CREATE TABLE '+relname+' ('
    for n in tab.dtype.names:
        dtype = str(tab[n].dtype)
        if dtype in dtypes:
            dbtype = dtypes[dtype]
        elif dtype.startswith('|S'):
                dbtype = 'varchar('+dtype[2:]+')'
        print n,dtype, dbtype
        q += n.lower().replace('-','_') + ' ' + dbtype + ',\n'
    q = q.rstrip(',\n') + ');\n'
    ds.execute(q)
    
    # fill in data
    values = ''
    for t in tab:
        valstr = '\t'.join([str(ti) for ti in t])
        #valstr = valstr.lstrip('(')
        #valstr = valstr.rstrip(')')
        valstr = valstr.replace('nan','\N')
        values += valstr+'\n'
    f = StringIO.StringIO(values)
    ds.cursor.copy_from(f,relname)
    ds.commit()
   
if __name__ == '__main__':
    #fits2db('/work2/jwe/cat/Hauck1997.fit','hauck')
    fits2db('/home/jwe/Downloads/Mermilliod.fit','mermilliod')