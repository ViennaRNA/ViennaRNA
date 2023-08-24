include $(top_srcdir)/silent_rules.mk

pkgdata_DATA = \
    DoxygenLayout.xml

if WITH_REFERENCE_MANUAL

REFERENCE_MANUAL_FILES_XML = xml/*

##--------------------------------------------------##
## prepare all rules in case reference manuals are  ##
## going to be (re)generated                        ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_BUILD

doxygen-xml: $(pkginclude_HEADERS)
	$(doxygen_verbose)@doxygen@ @DOXYGEN_CONF@ >>doxyref.log 2>&1; \
  touch doxygen-xml

##--------------------------------------------------##
## In case the XML manual should be created, here   ##
## is the rule how to do so                         ##
##--------------------------------------------------##
$(REFERENCE_MANUAL_FILES_XML): doxygen-xml

endif WITH_REFERENCE_MANUAL_BUILD

endif WITH_REFERENCE_MANUAL

noinst_DATA = $(REFERENCE_MANUAL_FILES_XML)

EXTRA_DIST =  DoxygenLayout.xml \
              doxygen-xml \
              xml

# remove all generated content
maintainer-clean-local:
	-rm -rf *.log xml doxygen.warnings doxygen-xml

distclean-local:
	-rm -rf *.log xml doxygen.warnings doxygen-xml
