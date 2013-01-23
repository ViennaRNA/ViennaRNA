pkgdata_DATA = viennarna.bib

if WITH_REFERENCE_MANUAL

##--------------------------------------------------##
## prepare variables in case HTML reference manual  ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_HTML
    REFERENCE_MANUAL_FILES_HTML = html/*
    docdir_htmldir = $(docdir)/html
    docdir_html_DATA = $(REFERENCE_MANUAL_FILES_HTML)
    dist_docdir_htmldir = $(docdir)/html
    dist_docdir_html_DATA = $(REFERENCE_MANUAL_FILES_HTML)

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

##--------------------------------------------------##
## the tagfile is used as indicator to decide       ##
## whether the reference manual should be generated ##
##--------------------------------------------------##
  REFERENCE_MANUAL_FILES_TAG = @REFERENCE_MANUAL_TAGFILE@

##--------------------------------------------------##
## This is the rule how to generate the tagfile     ##
## and the actual reference manual sources          ##
##                                                  ##
## here we need to chdir into the output directory  ##
## and get the bibtex file since otherwise doxygen  ##
## would not link the references correctly          ##
##--------------------------------------------------##
$(REFERENCE_MANUAL_FILES_TAG):  @DOXYGEN_CONF@ $(pkginclude_HEADERS)
	@doxygen@ @DOXYGEN_CONF@

##--------------------------------------------------##
## In case the HTML manual should be created, here  ##
## is the rule how to do so                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_HTML

$(REFERENCE_MANUAL_FILES_HTML):  $(REFERENCE_MANUAL_FILES_TAG)

endif WITH_REFERENCE_MANUAL_HTML

##--------------------------------------------------##
## In case the PDF manual should be created, here   ##
## is the rule how to do so                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_PDF

$(REFERENCE_MANUAL_FILES_PDF): $(REFERENCE_MANUAL_FILES_TAG)
	cd latex; \
	make; \
	mv refman.pdf ../$(REFERENCE_MANUAL_FILES_PDF)

endif WITH_REFERENCE_MANUAL_PDF

endif WITH_REFERENCE_MANUAL_BUILD

##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
  doc_DATA = $(REFERENCE_MANUAL_FILES_TAG) $(REFERENCE_MANUAL_FILES_PDF)
  dist_doc_DATA = $(REFERENCE_MANUAL_FILES_TAG) $(REFERENCE_MANUAL_FILES_PDF)

endif WITH_REFERENCE_MANUAL




EXTRA_DIST = $(pkgdata_DATA) mainpage.h

