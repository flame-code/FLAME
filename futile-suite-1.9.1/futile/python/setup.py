from distutils.core import setup
import setuptools
# This is a list of files to install, and where
# (relative to the 'root' dir, where setup.py is)
# You could be more specific.
files = ["things/*"]

setup(name="futile",
      version="1.9",
      description="Python module for futile specs and data analysis",
      author="Luigi Genovese",
      author_email="luigi.genovese@cea.fr",
      url="www.bigdft.org",
      license='GNU-LGPL',
      packages=setuptools.find_packages(),
      install_requires=['matplotlib', 'numpy'],
      # setup_requires = ['pytest-runner'],
      tests_require=['nbval'],
      # 'package' package must contain files (see list above)
      # I called the package 'package' thus cleverly confusing
      # the whole issue...
      # This dict maps the package name =to=> directories
      # It says, package *needs* these files.
      # package_data = {'package' : files },
      include_package_data=True,
      # 'runner' is in the root.
      # scripts = ["runner"],
      long_description="""
Welcome to the PyFUTILE module. This module contains python programs which
are useful for the specification of the Fortran FUTILE package.
""",
      # This next part it for the Cheese Shop, look a little down the page.
      classifiers=[
      # How mature is this project? Common values are
      # 3 - Alpha
      # 4 - Beta
      # 5 - Production/Stable
      'Development Status :: 4 - Beta',

      # Indicate who your project is intended for
      'Intended Audience :: End Users/Developers',
      'Topic :: Software Development :: Analysis and Running Tools',

      # Pick your license as you wish (should match "license" above)
      'License :: OSI Approved :: GNU Lesser General Public License (LGPL)',

      # Specify the Python versions you support here. In particular, ensure
      # that you indicate whether you support Python 2, Python 3 or both.
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      ]
      )
