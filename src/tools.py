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
        import os,statvfs

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
    
def log(logfile, logstring='', silent=False, withtime=True):
    '''
    Created on Jan 14, 2013

    @author: jwe
    Prints logstring to a logfile
    '''
    raise(DeprecationWarning)
    from datetime import datetime

    now = datetime.now()
    try:
        f = open(logfile,'at')
        if withtime:
            f.write('%s\t%s\n' % (now, logstring.rstrip('\n')))
        else:
            f.write('%s\n' % (logstring.rstrip('\n')))
        f.close()
    except IOError:
        print 'Can''t write logfile %s!' % logfile
    if not silent:
        print logstring    


class Estimator(object):
    '''
    Estimator object to provide an estimate and a percentage of a process
    '''


    def __init__(self, firstitem, lastitem):
        '''
        Constructor
        '''
        from time import time
        
        self.firstitem = firstitem
        self.lastitem = lastitem
        self.s0 = time()
    
    def estimate(self, current, comment=''):
        from time import time
        s1 = time()
        if current < self.firstitem:
            self.firstitem = current  
        totalitems = self.lastitem - self.firstitem
        percent = 100.*(current-self.firstitem)/totalitems
        perfile = (s1-self.s0)/(current-self.firstitem+1)
        togo = (totalitems - current + self.firstitem)*perfile
        hours = int(togo // 3600)
        mins =  int((togo % 3600) // 60)
        secs = (togo % 3600) % 60
        togostr = '%02d:%02d:%04.1f' % ( hours, mins, secs)
        print "%5.2f%%, %.3fs/item, %s: %s" % (percent, perfile, togostr, comment)
