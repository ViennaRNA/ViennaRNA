AUTOMAKE_OPTIONS = foreign subdir-objects

ACLOCAL_AMFLAGS = -I m4

## Make distcheck is broken in RNAforester due to g2 Makefile.in
AM_DISTCHECK_CONFIGURE_FLAGS = \
    --with-cluster \
    --without-forester \
    --with-kinwalker \
    --disable-silent-rules

SUBDIRS = \
    src \
    man \
    misc \
    doc/doxygen \
    interfaces \
    doc \
    RNA-Tutorial \
    @subdirs@ \
    tests \
    examples \
    packaging/macosx

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = RNAlib2.pc

EXTRA_DIST =  packaging/debian \
              packaging/make_windows_installer.sh \
              license.txt \
              README.md \
              CHANGELOG.md \
              MANIFEST.in

## remove all generated content
maintainer-clean-local:
	-rm -rf dist/ ViennaRNA.egg-info

distclean-local:
	-rm -f pyproject.toml setup.py setup.cfg
