#!/usr/bin/env python
# -*- coding: us-ascii -*-

#--------------------------------------------------------------------------------
# Copyright (C) 2015-2016 BigDFT group
# This file is distributed under the terms of the
# GNU General Public License, see
# or http://www.gnu.org/copyleft/gpl.txt .
#--------------------------------------------------------------------------------
from __future__ import print_function

BIGDFT_CFG='BIGDFT_CONFIGURE_FLAGS'
CLEAN=' clean '
CLEANONE=' cleanone '
UNINSTALL=' uninstall '
LIST=' list '
BUILD=' build '
BUILDONE=' buildone '
TINDERBOX=' tinderbox -o buildlogs '
DOT=' dot '
DOTCMD=' | dot -Edir=back -Tpng > buildprocedure.png '
DIST='  dist --dist-only '#bigdft-suite '
RCFILE='buildrc'
MKFILE='Makefile'
SETUP=' setup '

def get_macros(fr, *expressions):
    """
    Given a file, find the macros in it.

    Args:
      fr (str): path of the file to search in.
      *expressions: list of expressions to match.
    Returns:
      (list): a list of macros (strings).
    """
    from re import compile, search

    regex = [compile(x.strip()) for x in expressions]
    result = []
    with open(fr) as ifile:
        for line in ifile:
            if all(search(r, line) for r in regex) \
              and "dnl" not in line and '#' not in line:
                result.append(line.rstrip('\n'))
    return result

def get_files(fd, *expressions):
    """
    Given a directory, find the files which match the given macros.

    Args:
      fd (str): path to the directory to look in.
      *expressions: a list of expressions to match.

    Returns:
      (list): a list of files which have lines that match the expressions.
    """
    from os.path import join
    from re import compile, search
    from glob import glob

    regex = [compile(x.strip()) for x in expressions]
    result = []
    for f in glob(join(fd,"*")):
        with open(f) as ifile:
            for line in ifile:
                if all(search(r, line) for r in regex) \
                  and "dnl" not in line and '#' not in line:
                    result.append(f)
                    break

    return result

from copy import deepcopy
CHECKMODULES= ['futile','atlab','chess','psolver','bigdft','PyBigDFT','spred']

MAKEMODULES= deepcopy(CHECKMODULES)
MAKEMODULES.insert(CHECKMODULES.index('bigdft'),'libABINIT')

MANUALMAKEMODULES = deepcopy(MAKEMODULES)
MANUALMAKEMODULES.remove('PyBigDFT')
#allowed actions and corresponding description
ACTIONS={'build':
         'Compile and install the code with the given configuration.',
         'make':
         'Recompile the bigdft internal branches, skip configuring step.',
         'clean':
         'Clean the branches for a fresh reinstall.',
         'startover':
         'Wipe out all the build directories and recompile the important parts',
         'autogen':
         'Perform the autogen in the modules which need that. For developers only.',
         'update':
         'Useful to update a pre-compiled branch after a merge',
         'dist':
         'Creates a tarfile for the suite tailored to reproduce the compilation options specified.',
         'check':
         'Perform check in the bigdft branches, skip external libraries.',
         'cleanone':
         'Clean a single module of the suite',
         'buildone':
         'Build a single module of the suite',
         'dry_run':
         "Visualize the list of modules that will be compiled with the provided configuration in the 'buildprocedure.png' file.",
         'link':
         'Show the linking line that have to be used to connect an external executable to the package (when applicable)' }

#actions which need rcfile to be executed
NEEDRC=['build','dist','dry_run','startover','buildone']

#actions which do not need Makefile creation
NOMAKEFILE=['autogen','dry_run']

TARGETS={
    'bigdft': ['bin','bigdft'],
    'spred': ['bin','mhgps'],
    'chess': ['lib','libCheSS-1.a'],
    'futile': ['lib','libfutile-1.a'],
    'PyBigDFT': ['bin','bigdft'], #not really a target
    'psolver': ['lib','libPSolver-1.a'],
    'atlab': ['lib','libatlab-1.a']
    }

class BigDFTInstaller():
    """Class for the installation of the ``bigdft-suite``.

    Attributes:
        action (str): Action to be performed.
        package (str): Package to install.
        rcfile (str): Configuration file.
        conditions (list(str)): List of conditions.
        verbose (bool): If ``True``, more verbose.
        quiet (bool): If ``True``, no messages.
        yes (bool): If ``True``, ask a question.
    """
    m4_re=['^AX_','CHECK_PYTHON','PKG_CHECK_MODULES'] #regular expressions to identify proprietary macros
    def __init__(self,action,package,rcfile,conditions,verbose,quiet,yes):
        import os
        from sys import executable as pythonpath
        self.action=action    #Action to be performed
        self.conditions = ['testing'] if conditions is None else conditions
        self.package=package  #Package
        self.yes=yes      #Ask a question
        #the position of the script
        self.scriptpath = __file__
        #look where we are
        self.srcdir = os.path.abspath(os.path.dirname(self.scriptpath))
        if self.srcdir == '': self.srcdir='.'
        #look the builddir
        self.builddir=os.getcwd()
        #look if we are building from a branch
        bigdftdir=os.path.join(self.srcdir,'bigdft')
        self.branch=os.path.isfile(os.path.join(bigdftdir,'branchfile'))
        self.verbose=(verbose or action=='check')  #verbose option
        if not self.verbose and not quiet: self.verbose=self.branch
        #To be done BEFORE any exit instruction in __init__ (get_rcfile)
        self.time0 = None

        if os.path.abspath(self.srcdir) == os.path.abspath(self.builddir) and self.action not in NOMAKEFILE:
            print(50*'-')
            print("ERROR: BigDFT Installer works better with a build directory different from the source directory, install from another directory")
            print("SOLUTION: Create a separate directory and invoke this script from it")
            exit(1)
        #hostname
        self.hostname=os.uname()[1]

        #rcfile
        self.get_rcfile(rcfile)

        #jhbuild script
        self.jhb=pythonpath+' '+os.path.join(self.srcdir,'bundler/jhbuild.py ')
        if self.rcfile != '': self.jhb += '-f '+self.rcfile

        #conditions to be added
        if len(self.conditions) > 0: self.jhb += ' --conditions=+'+self.conditions[0]
        for cond in self.conditions[1:]:
            self.jhb += ',+'+cond

        #date of bigdft executable if present
        self.time0=self.target_time()

        self.print_present_configuration()

        #now get the list of modules that has to be treated with the given command
        self.modulelist=self.get_output(self.jhb + LIST+self.package).split('\n')
        print(" List of modules to be treated:",self.modulelist)

        #then choose the actions to be taken
        getattr(self,action)()

    def target_time(self):
        import os
        dt=TARGETS['bigdft']
        dt=TARGETS.get(self.package)
        tgt=os.path.join(dt[0],dt[1])
        return self.filename_time(os.path.join(self.builddir,'install',tgt))

    def filename_time(self,filename):
        import os
        if os.path.isfile(filename):
            return os.path.getmtime(filename)
        else:
            return 0

    def hook_environment_modifications(self):
        """
        Intercepts the `py:func:addpath` method of jhbuild
        in order to understand which environment variables have been
        identified to compile the suite.
        """
        pass

    def get_rcfile(self, rcfile):
        """Determine the rcfile"""
        import os
        from six.moves import input
        # see if the environment variables BIGDFT_CFG is present
        self.rcfile = ''
        if rcfile is not None:
            self.rcfile = rcfile
        else:
            if not (BIGDFT_CFG in list(os.environ.keys())) or self.action in NEEDRC: self.rcfile=RCFILE
        # see if it exists where specified
        if os.path.exists(self.rcfile): return
        # otherwise search again in the rcfiles
        rcdir=os.path.join(self.srcdir,'rcfiles')
        if self.rcfile != '': self.rcfile=os.path.join(rcdir, self.rcfile)
        if os.path.exists(self.rcfile): return
        self.rcfile = ''
        if BIGDFT_CFG in list(os.environ.keys()): return
        #otherwise search for rcfiles similar to hostname and propose a choice
        rcs=[]
        for file in os.listdir(rcdir):
            testname=os.path.basename(file)
            base=os.path.splitext(testname)[0]
            if base in self.hostname or self.hostname in base or base.split('-')[0] in self.hostname: rcs.append(file)
        print("Search in the configuration directory '%s'" % rcdir)
        if len(rcs)==1:
            self.rcfile=os.path.join(rcdir,rcs[0])
        elif len(rcs) > 0 and (self.action in NEEDRC or not self.yes):
            print("No valid configuration file specified, found various that matches the hostname '%s'" % self.hostname)
            print('In the directory "'+rcdir+'"')
            print('Choose among the following options')
            for i,rc in enumerate(rcs):
                print(str(i+1)+'. '+rc)
            while True:
                choice=input('Pick your choice (q to quit) ')
                if choice == 'q': exit(0)
                try:
                    ival=int(choice)
                    if (ival <= 0): raise
                    ch=rcs[ival-1]
                    break
                except:
                    print('The choice must be a valid integer among the above')
            self.rcfile=os.path.join(rcdir,ch)
        elif self.action in NEEDRC:
            print('No valid configuration file provided and '+BIGDFT_CFG+' variable not present, exiting...')
            exit(1)

    def __dump(self,*msg, **kwargs):
        if self.verbose and kwargs.get('verbose',True):
            for m in msg:
                print (m)

    def print_present_configuration(self):
        import  os
        from six.moves import input
        indent = ' '*2
        print('Configuration chosen for the Installer:')
        print(indent + 'Hostname:',self.hostname)
        print(indent + 'Source directory:',os.path.abspath(self.srcdir))
        print(indent + 'Compiling from a branch:',self.branch)
        print(indent + 'Build directory:',os.path.abspath(self.builddir))
        print(indent + 'Action chosen:',self.action)
        print(indent + 'Verbose:',self.verbose)
        print(indent + 'Jhbuild baseline:',self.jhb)
        if self.rcfile=='' and self.action in NEEDRC:
            print(indent + 'Configuration options:')
            print(indent*2 + "Source: Environment variable '%s'" % BIGDFT_CFG)
            print(indent*2 + "Value: '%s'" % os.environ[BIGDFT_CFG])
        elif self.rcfile!='':
            print(indent + 'Configuration options:')
            print(indent*2 + "Source: Configuration file '%s'" % os.path.abspath(self.rcfile))
        while not self.yes:
            if not self.branch:
                print(indent + '#WARNING: You are compiling from a User Branch. Developments are discouraged in this case')
                print(indent + '# as compilation errors might lead to source deletion. For extensive developments the usage')
                print(indent + '# of a versioned branch is advised. Please ignore this warning if you are not a developer.')
            ok = input('Do you want to continue (Y/n)? ')
            if ok == 'n' or ok=='N':
                exit(0)
            elif ok != 'y' and ok != 'Y' and repr(ok) != repr(''):
                print('Please answer y or n')
            else:
                break

    def selected(self,l):
        return [val for val in l if val in self.modulelist]

    def jhbuildaction(self, action, target, options=''):
        "Call jhbuild to do action on target."
        import os
        command = action + target +options
        ret = os.system(self.jhb + command)
        if (ret != 0):
            raise RuntimeError('JHBuild failed to complete :' + command )

    def shellaction(self,path,modules,action,hidden=False):
        "Perform a shell action, dump also the result if verbose is True."
        import os
        import sys
        for mod in self.selected(modules):
            directory=os.path.join(path,mod)
            here = os.getcwd()
            if os.path.isdir(directory):
                #self.__dump('Treating directory '+directory)
                sys.stdout.write('Module '+mod+' ['+directory+']: '+action)
                sys.stdout.flush()
                os.chdir(directory)
                if hidden:
                    self.get_output(action)
                else:
                    ierr=os.system(action)
                    if ierr != 0: raise Exception('Error in action: "'+action+'" for package: "'+mod+'"')
                    #print 'Error in action: "'+action+'" for package: "'+mod+'"'
                    #sys.exit(1)
                os.chdir(here)
                #self.__dump('done.')
                sys.stdout.write(' (done)\n')
            else:
                sys.stdout.write('Cannot perform action "'+action+'" on module "'+mod+'" directory not present in the build.\n')
            sys.stdout.flush()

    def get_output(self,cmd,verbose=True):
        import subprocess
        self.__dump('executing:',cmd,verbose = verbose)
        # note the use of universal_newlines
        # in python3, this forces the output to plain text instead of binary.
        # this poorly named parameter was aliased as "text" in version 3.7
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True,universal_newlines=True)
        (out, err) = proc.communicate()
        self.__dump("program output:", out,verbose = verbose or len(out) > 0)
        return out

    def removefile(self,pattern,dirname,names):
        "Delete the files matching the pattern"
        import os,fnmatch
        for name in names:
            if fnmatch.fnmatch(name,pattern):
                self.__dump('removing',os.path.join(dirname,name))
                os.remove(os.path.join(dirname,name))

    def get_ac_argument(self,tgt,acmacro):
        "Retrieve the list of ac arguments for the macro acmacro from file tgt"
        import os
        m4args=set()
        if os.path.isfile(tgt):
            for dd in get_macros(tgt, acmacro):
                if len(dd) == 0: continue
                m4=dd.split('[')[1]
                m4args.add(m4.split(']')[0])
        return list(m4args)

    def get_m4_macros(self,tgt,previous_macros=[]):
        "Identify the name of the proprietary m4 macros used in configure.ac"
        import os
        macros=set()
        if os.path.isfile(tgt):
            for regexp in self.m4_re:
                for m4 in get_macros(tgt, regexp):
                    if len(m4)>0:
                        m4t=m4.split('(')[0]
                        if m4t not in previous_macros: macros.add(m4t)
            required=self.get_ac_argument(tgt,'AC_REQUIRE')
            for m in required:
                found=m in previous_macros
                if found: continue
                for regexp in self.m4_re:
                    found=regexp.lstrip('^') in m
                    if found: break
                if found: macros.add(m)
        return list(macros)

    def get_m4_files(self,macros):
        "Find the files needed for the definition of the proprietary macros"
        import os
        files=set()
        #localize then the associated file in the m4 repository
        tgt=os.path.join(self.srcdir,'m4')+os.sep
        #print 'XXXXX',macros
        for m in macros:
            ffs=get_files(tgt,m,'AC_DEFUN')[0]
            if ffs!='': files.add(ffs)
        #now for each of the files get all the macros which are required but not explicitly called
        files=list(files)
        newm4=set()
        for f in files:
            tgt=os.path.join(self.srcdir,f)
            newm=self.get_m4_macros(tgt,previous_macros=macros)
            #print 'new macros',newm,f
            for m4 in newm:
                newm4.add(m4)
        newm4 = list(newm4)
        #print "new macros which have to be added",newm4,files
        if len(newm4) > 0: files=self.get_m4_files(macros+newm4)
        return list(files)

    def get_m4_dir(self,mod):
        "Return the configure macro dir(s)"
        import os
        tgt=os.path.join(self.srcdir,mod,'configure.ac')
        return self.get_ac_argument(tgt,'AC_CONFIG_MACRO_DIR')

    def copyfiles(self,filelist,dest):
        import os,shutil
        if not os.path.isdir(dest): return
        for f in filelist:
            test=os.path.join(dest,os.path.basename(f))
            if os.path.isfile(test):
                if len(self.get_output('diff -q '+f+' '+test,verbose = False)) == 0: continue
            try:
                shutil.copy(f,dest)
            except:
                os.chmod(dest, 777) #?? still can raise exception
                shutil.copy(f,dest)
            print('Copied file '+f+' in directory '+dest)

    def autogen(self):
        "Perform the autogen action"
        import os
        #first copy the macros in the config.m4 directories of proprietary packages
        for mod in self.selected(MAKEMODULES):
            macros=self.get_m4_macros(os.path.join(self.srcdir,mod,'configure.ac'))
            #print 'initial macros',mod,macros
            files=self.get_m4_files(macros)
            #print 'related files',files
            #print 'directory to copy files',self.get_m4_dir(mod)
            #now copy the files, overwriting the previously existing ones
            for d in self.get_m4_dir(mod):
                self.copyfiles(files,os.path.join(self.srcdir,mod,d))
            self._setupone(mod)
        #self.jhbuildaction(SETUP, self.package)
        #self.shellaction(self.srcdir,self.modulelist,'autoreconf -fi')

    def check(self):
        "Perform the check action"
        self.shellaction('.',CHECKMODULES,'make check',hidden=not self.verbose)

    def make(self):
        "Perform the simple make action"
        self.shellaction('.',MANUALMAKEMODULES,'make -j6 && make install',hidden=not self.verbose)

    def dist(self):
        "Perform make dist action"
        import os
        tarfile=os.path.join(self.builddir,self.package+'-suite.tar.gz')
        disttime0=self.filename_time(tarfile)
        self.jhbuildaction(DIST, self.package+"-suite")
        disttime1=self.filename_time(tarfile)
        if not (disttime1 == disttime0):
            print('SUCCESS: distribution file "'+self.package+'-suite.tar.gz" generated correctly')
        else:
            print('WARNING: the dist file seems not have been updated or generated correctly')

    def build(self):
        "Build the bigdft module with the options provided by the rcfile"
        import os
        #in the case of a nonbranch case, like a dist build, force checkout
        #should the make would not work
        co = '' if self.branch else ' -C'
        if (self.verbose):
            self.jhbuildaction(BUILD, self.package,options = co)
        else:
            self.jhbuildaction(TINDERBOX, self.package, options = co)

    def _wipeout(self,mod):
        import shutil
        print('Wipe directory: ',mod)
        shutil.rmtree(mod, ignore_errors=True)

    def _cleanone(self,mod):
        import os
        self.jhbuildaction(UNINSTALL, mod)
        self.jhbuildaction(CLEANONE, mod)
        #here we should eliminate residual .mod files
        os.walk(mod,self.removefile,"*.mod")
        os.walk(mod,self.removefile,"*.MOD")
        if not self.branch: #this is necessary as we come from a tarfile
            self._wipeout(mod)

    def cleanone(self):
        self._cleanone(self.package)

    def buildone(self):
        self._buildone(self.package)

    def clean(self):
        "Clean files in the build directory"
        for mod in self.selected(MAKEMODULES[::-1]): #invert cleaning order for elegance
            self._cleanone(mod)

    def _buildone(self,mod):
        print('Resetting: ',mod)
        self._setupone(mod)
        print('Building: ',mod)
        self.jhbuildaction(BUILDONE,mod)

    def _setupone(self,mod):
        self.jhbuildaction(SETUP,mod,options=' -t '+mod)

    def startover(self):
        "Wipe files in the makemodules directory"
        if not self.branch:
            print('ERROR: The action "startover" is allowed only from a developer branch')
            exit(1)
        import os
        for mod in self.selected(MAKEMODULES):
            self.get_output(self.jhb+UNINSTALL+mod)
            self._wipeout(mod)
        print('Building again...')
        for mod in self.selected(MAKEMODULES):
            self._buildone(mod)
        self.build()

    def dry_run(self):
        "Do dry build"
        self.get_output(self.jhb+DOT+self.package+DOTCMD)

    def link(self):
        "Show the linking line, when applicable"
        import os
        PPATH="PKG_CONFIG_PATH"
        addpath=os.path.join(self.builddir,'install','lib','pkgconfig')
        if PPATH in os.environ:
            if addpath not in os.environ[PPATH].split(':'):
                os.environ[PPATH]+=':'+addpath
        else:
            os.environ[PPATH]=addpath
        includes=self.get_output('pkg-config --cflags '+self.package)
        libs=self.get_output('pkg-config --libs '+self.package)
        #add the external linalg at the end to avod linking problems
        linalg=self.get_output('pkg-config --variable linalglibs '+self.package)
        plugin=self.get_output('pkg-config --variable plugin '+self.package)
        print('--------- Linking line to build with package "'+self.package+'":')
        print("  "+includes+" "+libs)

    def makefile_dump(self):
        "Build the Makefile that the installation of BigDFT creates for performing the same actions on this build"
        import os
        sflist=[]
        sflist.append("""
        #This Makefile is automatically generated from the BigDFT installer to avoid the calling to
        #the installer again for future action on this build
        #Clearly such actions only _assume_ that the build is fully functional and almost nothing
        #can be done with this file if a problem might arise.
        #Otherwise stated: this is an automatic message, please do not reply.
all: build
        """)
        #this should be done only if buildrc is there
        for a in ACTIONS:
            sflist.append(a+':  ')
            sflist.append('\t'+
                          os.path.abspath(self.scriptpath)+' '+a+' '+self.package+' -y')
        mkfile=open(MKFILE,'w')
        for item in sflist:
            mkfile.write("%s\n" % item)
            #rcfile.write("\n")
        mkfile.close()


    def rcfile_from_env(self):
        "Build the rcfile information from the chosen "+BIGDFT_CFG+" environment variable"
        import os
        if os.path.isfile(self.rcfile) and not os.path.isfile(RCFILE):
            from shutil import copyfile
            copyfile(self.rcfile,RCFILE)
            print('The configuration file used has been copied in the build tree, file "'+RCFILE+'"')
            return
        if BIGDFT_CFG not in list(os.environ.keys()) or os.path.isfile(RCFILE): return
        print('The suite has been built from a single configure line.')
        rclist=[]
        rclist.append("""#This is the configuration file for the BigDFT installer""")
        rclist.append("""#This is a python script which is executed by the build suite """)
        rclist.append(" ")
        rclist.append("""#Add the condition testing to run tests and includes PyYaml""")
        for cond in self.conditions:
            rclist.append('conditions.add("'+cond+'")')
        rclist.append("""#List the module the this rcfile will build""")
        rclist.append("modules = ['"+self.package+"',]")
        sep=' """ '
        confline=sep+os.environ[BIGDFT_CFG]+sep
        rclist.append("#example of the potentialities of the python syntax in this file")
        rclist.append("def env_configuration():")
        rclist.append("    return "+confline)
        rclist.append("#the following command sets the environment variable to give these settings")
        rclist.append("#to all the modules")
        rclist.append("import os")
        rclist.append("os.environ['"+BIGDFT_CFG+"']=env_configuration()")
        rclist.append("#here follow the configuration instructions for the modules built")
        rclist.append("#we specify the configurations for the modules to customize the options if needed")
        rclist.append("module_autogenargs.update({")
        rclist.append("   ")
        for mod in self.modulelist:
            rclist.append("'"+mod+"': env_configuration(),")
            rclist.append("   ")
        rclist.append("})")
        #then write the file
        rcfile=open(RCFILE,'w')
        for item in rclist:
            rcfile.write("%s\n" % item)
            #rcfile.write("\n")
        rcfile.close()
        print("Your used configuration options have been saved in the file '%s'." % RCFILE)
        print("Such file will be used for next builds, you might also save it in the 'rcfiles/'.")
        print("Directory of the source for future use. The name might contain the hostname.")

    def __del__(self):
        import os
        print(50*'-')
        print('Thank you for using the Installer of BigDFT suite.')
        print('The action considered was:',self.action)
        try:
           if self.time0 is not None:
               if not (self.time0==self.target_time()) and self.target_time()!=0:
                   print('SUCCESS: The Installer seems to have built correctly',self.package,' bundle')
                   print('All the available executables and scripts can be found in the directory')
                   print('"'+os.path.join(os.path.abspath(self.builddir),'install','bin')+'"')
                   if self.action in NEEDRC: self.rcfile_from_env()
               elif (self.action == 'build' or self.action == 'make'):
                   print('WARNING: The Installer seems NOT have created or updated',self.package,' binaries')
                   print('        (maybe everything was already compiled?)')
                   print('ACTION: check the compiling procedure.')
                   if self.branch:
                       print('HINT: It appears you are compiling from a branch source tree. Did you perform the action "autogen"?')
                   if not self.verbose and self.action == 'build':
                      print('  HINT: Have a look at the file index.html of the buildlogs/ directory to find the reason')
                   if self.yes:
                      exit(1)
        except:
            print('Goodbye...')
        if self.action not in NOMAKEFILE: self.makefile_dump() #temporary

if __name__ == '__main__':
    #import the uniparse module
    import UniParse

    #Redefine ArgumentParser to have the help message if no arguments
    class Installer_Parser(UniParse.UniParser):
        def error(self, message):
            import sys
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            self.exit()

    parser=Installer_Parser(description='BigDFT suite Installer',
                                epilog='''
    If you want more help type "%(prog)s help"
    ------------------------------------------------
    For more information, visit www.bigdft.org''')

    parser.option('action',nargs='?',default='help',
                        help='Action to be performed by the Installer.'
                        ' (default: %(default)s)',choices=['help']+[a for a in ACTIONS])
    parser.option('package',nargs='?',default='spred',
                        help='Package to be built by the installer. (default: %(default)s)',
                        choices=CHECKMODULES)
    parser.option('-f','--file',
                       help='Use an alternative configuration file instead of the default configuration '
                        + 'given by the environment variable %s' % BIGDFT_CFG)

    parser.add_group(mutually_exclusive=True)
    parser.group_option("-v", "--verbose", action="store_true",help=
                        'Verbose output, default from a development branch')
    parser.group_option("-q", "--quiet", action="store_true",help=
                        'Verbosity disabled output, default from a development branch')

    parser.option('-d','--debug',action='store_true',
                  help='Verbose output, default from a development branch')
    parser.option('-a','--add-condition',
                  help='Add a condition for the compilation of the suite.')
    parser.option('-y','--yes',action='store_true',
                  help='Answer yes to dialog questions')
    parser.option('-c','--configure-line',remainder=True,
                  help='Specify the configure line to be passed (set BIGDFT_CONFIGURE_FLAGS variable)')

    args = parser.args()

    conds=None if args.add_condition is None else args.add_condition.split(',')

    if args.configure_line is not None:
      cfg=''
      for i in args.configure_line:
        if i is not None: cfg+=i+' '
      #scratch the BIGDFT_CFG environment variable
      import os
      os.environ[BIGDFT_CFG]=cfg

    if args.action=='help':
        print("Quick overview of the BigDFT suite Installer program")
        print(50*'-')
        print("USAGE: Installer.py <action> <package>")
        print(50*'-'+'Available actions')
        actions=list(ACTIONS.keys())
        actions.sort()
        for a in actions:
            print(a,':')
            print('     ',ACTIONS[a])
        print(50*'-')
        print('Available packages:',CHECKMODULES)
        print(50*'-')
        print(10*"QIFI-"+' (Quick Instructions For the Impatient)')
        print('Ideally, there are two different policies:')
        print('Developer: From a development branch, start by "autogen", then "build"')
        print('     User: From a tarball, start by "build"')
        print('Perform the "dry_run" command to have a graphical overview of the building procedure')
    elif args.action=='update':
        yes=args.yes
        for action in ['clean','autogen','build']:
            BigDFTInstaller(action,args.package,args.file,\
                            conds,args.verbose,args.quiet,yes)
            yes=True
    else:
        BigDFTInstaller(args.action,args.package,args.file,\
                        conds,args.verbose,args.quiet,args.yes)
