#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 24, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

IDL to Python converter
'''
def convert(idlfile, pythonfile=None, verbose=False):
    import os
    #import parse

    if pythonfile is None:
        # remove .pro ending and append .py ending
        oldbase = os.path.splitext(idlfile)
        print oldbase
        pythonfile = oldbase[0]+'.py'
    f = open(idlfile)
    lines = f. readlines()
    f.close()
    docstring = False
    newlines = ''
    for l in lines:
        if ';+' in l:
            docstring = True
        if ';-' in l:
            docstring = False
        
        l = l.replace(';+','"""')
        l = l.replace(';-','"""')
        if docstring: 
            l = l.replace(';',' ')
        else:
            l = l.replace(';','#')
            l = l.replace('^','**')
            l = l.replace('GT','>')
            l = l.replace('GE','>=')
            l = l.replace('EQ','==')
            l = l.replace('LT','<')
            l = l.replace('LE','<=')
            l = l.replace('$','\\')
            l = l.replace('then',':')
            l = l.replace('endif','#endif')
            l = l.replace('begin','#begin')
            l = l.replace('else','else:')
            l = l.replace('endelse','#endelse')
            l = l.replace('endelse','#endelse')
            l = l.replace('endfor','#endfor')
            l = l.replace('!DPI','pi')
            l = l.replace('pro ','def ')
            l = l.replace('function ','def ')
            l = l.replace('&&','and')
            l = l.replace('||','or')
            l = l.replace('&',';')
            l = l.replace('print,','print ')
        newlines += l
        if verbose:
            print l.rstrip('\n')
    f = open(pythonfile,'wt')
    f.writelines(newlines)
    f.close()
            
        

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='IDL to Python converter')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose printing')
    parser.add_argument('idlfile', help='name of the inputfile')
    parser.add_argument('pythonfile', nargs='?', default=None, help='name of the inputfile')

    args = parser.parse_args()    
    convert(args.idlfile,args.pythonfile, verbose=args.verbose)
