#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 1, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
a = 1
b = 0
print 'first try'
try:
    c = a/b
except ZeroDivisionError, e:
    print 'exception %s' % e
else:
    print 'else'
finally:
    print 'finally'
print 'end'
b = 1
print 'second try'
try:
    c = a/b
except ZeroDivisionError, e:
    print 'exception %s' % e
else:
    print 'else'
finally:
    print 'finally'
print 'end'
b=0
print 'third try'
try:
    c = a/b
except IOError, e:
    print 'exception %s' % e
else:
    print 'else'
finally:
    print 'finally'
print 'end'
