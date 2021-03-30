#!/usr/bin/env python
# -*- coding: us-ascii -*-

#--------------------------------------------------------------------------------
# Copyright (C) 2015-2016 BigDFT group
# This file is distributed under the terms of the
# GNU General Public License, see
# or http://www.gnu.org/copyleft/gpl.txt .
#--------------------------------------------------------------------------------
from __future__ import print_function

#Redefine ArgumentParser to have the help message if no arguments
class UniParser():
    def __init__(self,*args,**kwargs):
        argss,whichone=self.kw_pop('method','argparse',**kwargs)
        if whichone=='argparse':
            try:
                import argparse
                self.argparse=True
            except:
                import optparse
                self.argparse=False
        elif whichone=='optparse':
            import optparse
            self.argparse=False
        else:
            raise ValueError('method not recognized')
        if not self.argparse: self.positional=[]
        if self.argparse:
            arg = self.kw_update('formatter_class',argparse.RawDescriptionHelpFormatter,**argss)
            self.parser=argparse.ArgumentParser(*args,**arg)
        else:
            self.parser=optparse.OptionParser(*args,**argss)
    def option(self,*args,**kwargs):
        arg,remainder=self.kw_pop('remainder',False,**kwargs)
        if self.argparse:
            import argparse
            if remainder: arg.update({'nargs':argparse.REMAINDER})
            self.parser.add_argument(*args,**arg)
        else:
            name=args[0]
            if '-' not in name or remainder:
                if remainder: 
                    name=args[1].lstrip('--').replace('-','_')
                    self.parser.add_option(*args,**arg)
                positional=self.kw_update('name',name,**kwargs)
                self.positional.append(positional)
                if remainder: positional.update({'remainder':True})
            else:
                self.parser.add_option(*args,**arg)
    def add_group(self,*args,**kwargs):
        arg,mutually_exclusive=self.kw_pop('mutually_exclusive',False,**kwargs)
        if self.argparse:
            if mutually_exclusive:
                grp=self.parser.add_mutually_exclusive_group(*args,**arg)
            else:
                grp=self.parser.add_argument_group(*args,**arg)
        else:
            import optparse
            if len(args) == 0: args += ('group',)
            grp=optparse.OptionGroup(self.parser,*args,**arg)
        if hasattr(self,'groups'):                                
            self.groups.append(grp)
        else:
            self.groups=[grp]
    def group_option(self,*args,**kwargs):
        if self.argparse:
            self.groups[-1].add_argument(*args,**kwargs)
        else:
            self.groups[-1].add_option(*args,**kwargs)
    def kw_update(self,key,val,**kwargs):
        arg=kwargs.copy()
        arg.update({key: val})
        return arg
    def kw_pop(self,*args,**kwargs):
        arg=kwargs.copy()
        key,default=args
        if key in arg:
            return arg,arg.pop(key)
        else:
            return arg,default
    def error(self, message):
        import sys
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        self.exit()
    def args(self):
        meth='argparse' if self.argparse else 'optparse'
        print('Parsing arguments with method '+meth+'...')
        if self.argparse:
            return self.parser.parse_args()
        else:
            (argstmp, argtmp) = self.parser.parse_args()
            for i,pos in enumerate(self.positional):
                remainder=pos.get('remainder',False)
                if i < len(argtmp):
                    if remainder:
                        val=[argtmp[i]]
                        for j in range(i+1,len(argtmp)):
                            val.append(argtmp[j])
                    else:
                        val=argtmp[i]
                else:
                    val=pos.get('default',[])
                name=pos['name']
                if hasattr(argstmp,name): 
                    val=[getattr(argstmp,name)]+val
                setattr(argstmp,name,val)
            return argstmp

def test(method):
    parser= UniParser(description='BigDFT suite Installer',
                            epilog='''
If you want more help type "%(prog)s help"
------------------------------------------------
                            For more information, visit www.bigdft.org''',method=method)
    parser.option('action',nargs='?',default='help',
               help='Action to be performed by the Installer.'
               ' (default: %(default)s)',choices=['help'])
    parser.option('package',nargs='?',default='spred',
               help='Package to be built by the installer. (default: %(default)s)')
    parser.option('-f','--file',
               help='Use an alternative configuration file instead of the default configuration '
               + 'given by the environment variable <foo>')
    parser.add_group(mutually_exclusive=True)
    parser.group_option("-v", "--verbose", action="store_true",help='Verbose output, default from a development branch')
    parser.group_option("-q", "--quiet", action="store_true",help='Verbosity disabled output, default from a development branch')
    parser.option('-d','--debug',action='store_true',
                   help='Verbose output, default from a development branch')
    parser.option('-c','--configure-line',remainder=True, #nargs=argparse.REMAINDER,
                   help='Specify the configure line to be passed (set BIGDFT_CONFIGURE_FLAGS variable)')
    args = parser.args()
    return args

if __name__ == "__main__":

    args=test('argparse')
    print('Arguments argparse')
    print(args.__dict__)

    args=test('optparse')
    print('Arguments optparse')
    print(args.__dict__)
