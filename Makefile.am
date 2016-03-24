AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS = --with-cluster --without-forester --with-python3 --with-kinwalker --disable-silent-rules PYTHON_VERSION=2.7

SUBDIRS = src man misc doc interfaces @subdirs@ tests packaging/macosx

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = RNAlib2.pc

EXTRA_DIST =  INSTALL.configure \
              examples
