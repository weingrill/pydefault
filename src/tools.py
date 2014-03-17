'''
Created on Dec 16, 2013

@author: jwe
'''
class FileCache(str):
    """This shoud become a cached file object"""
    def __init__(self, sourcefile=None, cachefile=None):
        """simply returns the cachefile name,
        but the work is done in the background:
        .) copies file from source to cache
        .) checks free space
        .) copies back the file to free up some space
        """ 
        
        # do we have enough diskspace?
        import sys,os,statvfs

        f = os.statvfs("/home")
        print "preferred block size", "=>", f[statvfs.F_BSIZE]
        print "fundamental block size", "=>", f[statvfs.F_FRSIZE]
        print "total blocks", "=>", f[statvfs.F_BLOCKS]
        print "total free blocks", "=>", f[statvfs.F_BFREE]
        print "available blocks", "=>", f[statvfs.F_BAVAIL]
        # |--> no: move a file back
        # |--> yes: copy file from source
        self.cachefile = cachefile
        return self.cachefile
    
def log(logfile, logstring='', silent=False):
    '''
    Created on Jan 14, 2013

    @author: jwe
    Prints logstring to a logfile
    '''
    from datetime import datetime

    now = datetime.now()
    try:
        f = open(logfile,'at')
        f.write('%s\t%s\n' % (now, logstring.rstrip('\n')))
        f.close()
    except IOError:
        print 'Can''t write logfile %s!' % logfile
    if not silent:
        print logstring    
