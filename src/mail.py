#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 17, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Mail(object):
    '''
    classdocs
    '''


    def __init__(self, recipient, subject, text):
        '''
        Constructor
        '''

        import subprocess
        cmd = ('echo "%s" | mailx -s "%s" %s' % (text, subject, recipient))
        #LOG.info(cmd)
        subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]

if __name__ == '__main__':
    m = Mail('jweingrill@aip.de','Python mail', """Dies ist eine mail
    aus Python heraus""")