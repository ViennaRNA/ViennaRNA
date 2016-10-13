#!/usr/bin/env python
# Run 'make' or 'make install' to call this script.
# The result of 'make install' will generally be :
#    on Linux : /usr/lib/python2.4/site-packages/g2.so
#    on Windows : C:\Python24\lib\site-packages\g2.pyd
# Then in Python say : import g2

# For a module based on libg2.so, make sure libg2.so is in the link path
# (you may have to add a link like libg2.so -> libg2.so.0.0.70).
# For a stand alone module (not recommended), make sure libg2.a is in
# the link path (very likely, given -L../), and libg2.so is not.
# Alternatively, add '-static' to link_args below.

from distutils.core import setup, Extension

import sys

error_string = '\n No %s command line argument.' \
               '\n It should be a string of %s options.' \
               '\n Specify an empty string (\'\') if there are none.'

# set compile_args and link_args
for (pos, step) in (('1st', 'compile'), ('2nd', 'link')):
    try:
        globals()[step + '_args'] = sys.argv.pop(1).split()
    except IndexError:
        print error_string % (pos, step)
        sys.exit(1)

import platform

build_ext_opts = {}

# called without arguments, just build the module in directory g2_python
if len(sys.argv) == 1:
    build_ext_opts['force'] = True
    sys.argv.append('build_ext')

# on Linux, strip debugging info from lib
if platform.system() == 'Linux':
    link_args += ['-s']

if platform.python_compiler()[0:3] == 'GCC' and 'sparc' in platform.machine():
    link_args += ['-fPIC', '-Wl,-O2', '-Wl,-Map=g2.map']

g2base = Extension( 'g2',
                    sources = ['g2module.c'],
                    include_dirs = ['../src'],
                    library_dirs = ['../'],
                    libraries = ['g2'],
                    extra_compile_args = compile_args,
                    extra_link_args = link_args)

setup( name = 'PythonG2',
       author = 'Tijs Michels',
       author_email = 'tijs@users.sourceforge.net',
       url = 'http://g2.sourceforge.net/',
       license = 'GNU Lesser General Public License, version 2.1 or above',
       version = '0.1',
       description = 'Easy to use, portable yet powerful 2D graphics library',
       ext_modules = [g2base],
       options = { 'build_ext' : build_ext_opts })
