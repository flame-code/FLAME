#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
if sys.version_info[0] == 2:
    import __builtin__ as builtins
else:
    import builtins

#USE_CHECKOUT_SRC = True

#if USE_CHECKOUT_SRC:
#    sys.path.insert(0, '/home/caliste/local/jhbuild')
#    pkgdatadir = None
#    datadir = None
#    import jhbuild
#    srcdir = os.path.abspath(os.path.join(os.path.dirname(jhbuild.__file__), '..'))
#else:
#    pkgdatadir = "@pkgdatadir@"
#    datadir = "@datadir@"
#    srcdir = "@srcdir@"
#    if '@pythondir@' not in sys.path:
#        sys.path.insert(0, '@pythondir@')
#    try:
#        import jhbuild
#    except ImportError:
#        sys.path.insert(0, srcdir)
pkgdatadir = None
datadir = None
import jhbuild
srcdir = os.path.abspath(os.path.join(os.path.dirname(jhbuild.__file__), '../..'))

builtins.__dict__['PKGDATADIR'] = pkgdatadir
builtins.__dict__['DATADIR'] = datadir
builtins.__dict__['SRCDIR'] = srcdir

jhb_env = {}


def addpath_override(oldfunc, *args, **kwargs):
    global jhb_env
    variable = args[0]
    # value = args[1]
    # situation = jhb_env.setdefault(variable, [])
    # prepend = kwargs.get('prepend', True)
    # if prepend:
    #     situation.insert(0, value)
    # else:
    #     situation.append(value)
    import os
    oldenv = os.environ.get(variable, '').split(os.pathsep)
    jhbvalues = jhb_env.get(variable, [])
    # print oldenv, jhbvalues, variable
    oldfunc(*args, **kwargs)
    newenv = {i: val
              for i, val in enumerate(os.environ[variable].split(os.pathsep))}
    for i, val in newenv.items():
        if val in oldenv + jhbvalues:
            continue
        situation = jhb_env.setdefault(variable, [])
        situation.insert(-1 if i != 0 else 0, val)


import jhbuild.environment as env
old_addpath = env.addpath
from functools import partial
env.addpath = partial(addpath_override, old_addpath)

import jhbuild.main
jhbuild.main.main(sys.argv[1:])

print (jhb_env)
