ACLOCAL_AMFLAGS = -I m4

SUBDIRS = lib Progs man Utils H @subdirs@

pkgdata_DATA = rna_turner1999.par rna_turner2004.par rna_andronescu2007.par dna_mathews1999.par dna_mathews2004.par

include aminclude.am


if DX_COND_doc

if DX_COND_html
  DOX_HTML_FILES = $(DX_DOCDIR)/html/*
  $(DOX_HTML_FILES): doxygen-doc
  docdir_htmldir = $(docdir)/html
  dist_docdir_html_DATA = $(DOX_HTML_FILES)
endif

if DX_COND_pdf
  DOX_PDF = $(DX_DOCDIR)/$(DX_PROJECT).pdf
##  $(DOX_PDF): doxygen-doc
endif

if DX_COND_ps
  DOX_PS = $(DX_DOCDIR)/$(DX_PROJECT).ps
##  $(DOX_PS): doxygen-doc
endif

DOX_TAG = $(DX_DOCDIR)/$(DX_PROJECT).tag
dist_doc_DATA = $(DOX_PDF) $(DOX_PS) $(DOX_TAG)
endif

EXTRA_DIST = $(pkgdata_DATA) INSTALL.configure Readseq $(DX_CONFIG) mainpage.h @SVM_SOURCE_ARCHIVE@

# remove all generated reference documents
maintainer-clean-local:
	rm -rf $(DX_DOCDIR)

#install-data-hook:
#	if test x"$(DX_FLAG_html)" = x1; then\
#	  printf "Installing html documentation";\
#	  $(MKDIR_P) $(DEST_DIR)/$(docdir)/html; \
#	  for i in `ls $(DX_DOCDIR)/html`; do\
#	    mv $(DEST_DIR)/$(docdir)/$$i $(DEST_DIR)/$(docdir)/html/; printf ".";\
#	  done;\
#	  printf "done\n";\
#	fi

#uninstall-hook:
#	if test x"$(DX_FLAG_html)" = x1; then\
#	  printf "Uninstalling html documentation...";\
#	  rm -rf $(DEST_DIR)/$(docdir)/html;\
#	  printf "done\n";\
#	fi
