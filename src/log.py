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
    finally:
        f.close()
    if not silent:
        print logstring    
    