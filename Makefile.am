ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS = --with-cluster --with-check PYTHON_VERSION=2.7 --without-forester

SUBDIRS = src man misc doc interfaces @subdirs@ @CHECK_DIR@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = RNAlib2.pc

EXTRA_DIST = INSTALL.configure examples
