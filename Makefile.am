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


############################################
# Override 'distdir' automake rule to make
# sure that we build the reference manual
# including the Python documentation whenever
# we run 'make dist' and alike.
# For that, sphinx requires a working Python
# wrapper to RNAlib, which we need to build
# first. After that, we clean up the source
# tree again and continue with the regular
# 'distdir-am' rule.
############################################
distdir: $(BUILT_SOURCES)
	((test -f doc/refman-html && test -f doc/RNAlib-@VERSION@.pdf) || \
  (cd src/ViennaRNA && $(MAKE) $(AM_MAKEFLAGS) && \
  cd ../../man && $(MAKE) $(AM_MAKEFLAGS) && \
  cd ../doc/doxygen && $(MAKE) $(AM_MAKEFLAGS) && \
  cd ../source/man && $(MAKE) $(AM_MAKEFLAGS) && \
  cd ../../../interfaces/Python && $(MAKE) $(AM_MAKEFLAGS) && \
  cd ../../doc && $(MAKE) $(AM_MAKEFLAGS) && \
  cd ../interfaces/Python && $(MAKE) $(AM_MAKEFLAGS) clean && \
  cd ../../src/ViennaRNA && $(MAKE) $(AM_MAKEFLAGS) clean))
	$(MAKE) $(AM_MAKEFLAGS) distdir-am

## remove all generated content
maintainer-clean-local:
	-rm -rf dist/ ViennaRNA.egg-info

distclean-local:
	-rm -f pyproject.toml setup.py setup.cfg
