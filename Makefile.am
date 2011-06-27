ACLOCAL_AMFLAGS = -I m4

SUBDIRS = lib Progs man Utils H @subdirs@

pkgdata_DATA = rna_turner1999.par rna_turner2004.par rna_andronescu2007.par dna_mathews1999.par dna_mathews2004.par

##include aminclude.am
##
##
##if DX_COND_doc
##
##if DX_COND_html
##  DOX_HTML_FILES = $(DX_DOCDIR)/html/*
##  $(DOX_HTML_FILES): doxygen-doc
##  docdir_htmldir = $(docdir)/html
##  dist_docdir_html_DATA = $(DOX_HTML_FILES)
##endif
##
##if DX_COND_pdf
##  DOX_PDF = $(DX_DOCDIR)/$(DX_PROJECT).pdf
##endif
##
##if DX_COND_ps
##  DOX_PS = $(DX_DOCDIR)/$(DX_PROJECT).ps
##endif
##
##DOX_TAG = $(DX_DOCDIR)/$(DX_PROJECT).tag
##dist_doc_DATA = $(DOX_PDF) $(DOX_PS) $(DOX_TAG)
##
##endif

EXTRA_DIST = $(pkgdata_DATA) INSTALL.configure Readseq mainpage.h @SVM_SOURCE_ARCHIVE@

## remove all generated reference documents
#maintainer-clean-local:
#	rm -rf $(DX_DOCDIR)
#
