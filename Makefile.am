AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS = --with-cluster --without-forester --with-python3 --with-kinwalker --disable-silent-rules PYTHON_VERSION=2.7

SUBDIRS = src man misc doc interfaces @subdirs@ tests @MACOSX_INSTALLER@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = RNAlib2.pc


AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS) @FLOAT_PF_FLAG@
LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la @CHECK_LIBS@

EXTRA_DIST =  INSTALL.configure \
              examples
