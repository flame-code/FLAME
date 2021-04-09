# jhbuild - a tool to ease building collections of source packages
# Copyright (C) 2001-2006  James Henstridge
# Copyright (C) 2007-2008  Frederic Peters
#
#   makesys.py: using a simple Makefile and make.sys config file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

__metaclass__ = type

import os
import stat
import shutil
try:
    import hashlib
except ImportError:
    import md5 as hashlib

from jhbuild.errors import FatalError, BuildStateError, CommandError
from jhbuild.modtypes import \
     DownloadableModule, register_module_type, MakeModule
from jhbuild.versioncontrol.tarball import TarballBranch

__all__ = [ 'MakesysModule' ]

class MakesysModule(MakeModule, DownloadableModule):
    '''Base type for modules that are distributed with a single Makefile
    with configuration inside a make.sys file.'''
    type = 'makesys'

    PHASE_CHECKOUT = DownloadableModule.PHASE_CHECKOUT
    PHASE_FORCE_CHECKOUT = DownloadableModule.PHASE_FORCE_CHECKOUT
    PHASE_CLEAN          = 'clean'
    PHASE_BUILD          = 'build'
    PHASE_CHECK          = 'check'
    PHASE_DIST           = 'dist'
    PHASE_DISTCLEAN      = 'distclean'
    PHASE_INSTALL        = 'install'

    def __init__(self, name, branch=None,
                 makeargs='',
                 makeinstallargs='',
                 makefile='Makefile',
                 preset=None):
        MakeModule.__init__(self, name, branch=branch, makeargs=makeargs,
                            makeinstallargs=makeinstallargs, makefile=makefile)
        self.preset = preset
        self.supports_install_destdir = True
        self.supports_non_srcdir_builds = False

    def get_srcdir(self, buildscript):
        return self.branch.srcdir

    def get_builddir(self, buildscript):
        if buildscript.config.buildroot and self.supports_non_srcdir_builds:
            d = buildscript.config.builddir_pattern % (
                self.branch.checkoutdir or self.branch.get_module_basename())
            return os.path.join(buildscript.config.buildroot, d)
        else:
            return self.get_srcdir(buildscript)

    def do_build(self, buildscript):
        buildscript.set_action(_('Building'), self)
        if self.preset is not None:
            shutil.copy(os.path.join(self.get_srcdir(buildscript), "config", "make.sys." + self.preset), os.path.join(self.get_srcdir(buildscript)))
            makeargs = ""
        else:
            open(os.path.join(self.get_srcdir(buildscript), "make.sys"), 'w').close()
            makeargs = self.get_makeargs(buildscript)
        cmd = '%s %s' % (os.environ.get('MAKE', 'make'), makeargs)
        buildscript.execute(cmd, cwd = self.get_builddir(buildscript),
                extra_env = self.extra_env)
    do_build.depends = [PHASE_CHECKOUT]
    do_build.error_phases = [PHASE_FORCE_CHECKOUT,
            PHASE_CLEAN, PHASE_DISTCLEAN]

    def do_clean(self, buildscript):
        buildscript.set_action(_('Cleaning'), self)
        makeargs = self.get_makeargs(buildscript)
        cmd = '%s clean' % (os.environ.get('MAKE', 'make'))
        buildscript.execute(cmd, cwd = self.branch.srcdir,
                extra_env = self.extra_env)

    def do_install(self, buildscript):
        buildscript.set_action(_('Installing'), self)
        destdir = self.prepare_installroot(buildscript)
        cmd = '%s install DESTDIR=%s prefix=%s' % (os.environ.get('MAKE', 'make'),
                                                   destdir, self.config.prefix)
        buildscript.execute(cmd, cwd = self.get_builddir(buildscript),
                    extra_env = self.extra_env)
        self.process_install(buildscript, self.get_revision())
    do_install.depends = [PHASE_BUILD]

    def xml_tag_and_attrs(self):
        return ('makesys',
                [('id', 'name', None),
                 ('preset', 'preset make.sys', None)])

def collect_args(instance, node, argtype):
    if node.hasAttribute(argtype):
        args = node.getAttribute(argtype)
    else:
        args = ''

    for child in node.childNodes:
        if child.nodeType == child.ELEMENT_NODE and child.nodeName == argtype:
            if not child.hasAttribute('value'):
                raise FatalError(_("<%s/> tag must contain value=''") % argtype)
            args += ' ' + child.getAttribute('value')

    return instance.eval_args(args)

def parse_makesys(node, config, uri, repositories, default_repo):
    instance = MakesysModule.parse_from_xml(node, config, uri, repositories, default_repo)

    if node.hasAttribute("preset"):
        instance.preset = node.getAttribute('preset')
    instance.makeargs = collect_args (instance, node, 'makeargs')
    instance.makeinstallargs = collect_args (instance, node, 'makeinstallargs')

    return instance
register_module_type('makesys', parse_makesys)

