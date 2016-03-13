ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS = --with-cluster --with-check --without-forester --with-kinwalker --disable-silent-rules PYTHON_VERSION=2.7

SUBDIRS = src man misc doc interfaces @subdirs@ @CHECK_DIR@ @MACOSX_INSTALLER@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = RNAlib2.pc

EXTRA_DIST = INSTALL.configure examples
