AUTOMAKE_OPTIONS = foreign subdir-objects

ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS = --with-cluster --without-forester --with-python3 --with-kinwalker --disable-silent-rules

SUBDIRS = src man misc doc RNA-Tutorial interfaces @subdirs@ tests examples packaging/macosx

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = RNAlib2.pc

EXTRA_DIST =  packaging/debian \
              packaging/make_windows_installer.sh \
              license.txt \
              README.md \
              CHANGELOG.md

