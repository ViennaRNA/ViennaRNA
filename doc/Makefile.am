pkgdata_DATA = viennarna.bib

if WITH_REFERENCE_MANUAL

##--------------------------------------------------##
## prepare variables in case HTML reference manual  ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_HTML
    REFERENCE_MANUAL_FILES_HTML = html/* 
endif WITH_REFERENCE_MANUAL_HTML

##--------------------------------------------------##
## prepare variables in case PDF reference manual   ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_PDF
    REFERENCE_MANUAL_FILES_PDF = @REFERENCE_MANUAL_PDF_NAME@
endif WITH_REFERENCE_MANUAL_PDF

##--------------------------------------------------##
## prepare all rules in case reference manuals are  ##
## going to be (re)generated                        ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_BUILD

dox:  @DOXYGEN_CONF@ $(pkginclude_HEADERS)
	$(AM_V_GEN)@doxygen@ @DOXYGEN_CONF@

##--------------------------------------------------##
## In case the HTML manual should be created, here  ##
## is the rule how to do so                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_HTML

$(REFERENCE_MANUAL_FILES_HTML): dox

endif WITH_REFERENCE_MANUAL_HTML

##--------------------------------------------------##
## In case the PDF manual should be created, here   ##
## is the rule how to do so                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_PDF

$(REFERENCE_MANUAL_FILES_PDF): dox
	$(AM_V_GEN)cd latex; make; mv refman.pdf ../$(REFERENCE_MANUAL_FILES_PDF)

endif WITH_REFERENCE_MANUAL_PDF

endif WITH_REFERENCE_MANUAL_BUILD

##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
  html_DATA = $(REFERENCE_MANUAL_FILES_HTML)
  pdf_DATA = $(REFERENCE_MANUAL_FILES_PDF)

endif WITH_REFERENCE_MANUAL


if WITH_REFERENCE_MANUAL_HTML

vrna_htmlrefman_inst: installdirs html
	cp -r html/* $(DESTDIR)$(htmldir)

vrna_htmlrefman_uninst:
	rm -rf $(DESTDIR)$(htmldir)

else

vrna_htmlrefman_inst:
vrna_htmlrefman_uninst:

endif

install-htmlDATA: vrna_htmlrefman_inst
uninstall-htmlDATA: vrna_htmlrefman_uninst

EXTRA_DIST = $(pkgdata_DATA) refman.include html latex @REFERENCE_MANUAL_PDF_NAME@

# remove all generated content
maintainer-clean-local:
	-rm -f *.tag *.pdf html/ latex/

