#!/bin/sh
# 
# Run this before configure
#
# This file blatantly ripped off from subversion.
#
# Note: this dependency on Perl is fine: only SVN developers use autogen.sh
#       and we can state that dev people need Perl on their machine
#

# [ -d ./Kinfold ] ||  mkdir Kinfold
# [ -d ./RNAforester ] ||  mkdir RNAforester

set -e

#automake --version | perl -ne 'if (/\(GNU automake\) ([0-9].[0-9])/) {print;  if ($1 < 1.4) {exit 1;}}'

#if [ $? -ne 0 ]; then
#    echo "Error: you need automake 1.4 or later.  Please upgrade."
#    exit 1
#fi

# Produce aclocal.m4, so autoconf gets the automake macros it needs
# echo "Creating aclocal.m4..."
#aclocal -I m4

# autoheader

# Produce all the `Makefile.in's, verbosely, and create neat missing things
# like `libtool', `install-sh', etc.
# automake --add-missing --gnu

# If there's a config.cache file, we may need to delete it.  
# If we have an existing configure script, save a copy for comparison.
if [ -f config.cache ] && [ -f configure ]; then
  cp configure configure.$$.tmp
fi

# Produce ./configure
echo "Creating configure..."
autoreconf -vi -I config

echo ""
echo "You can run ./configure now."
echo ""

