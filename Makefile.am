ACLOCAL_AMFLAGS = -I m4

SUBDIRS = lib Progs man Utils H @subdirs@

pkgdata_DATA = rna_turner1999.par rna_turner2004.par rna_andronescu2007.par dna_mathews1999.par dna_mathews2004.par

include aminclude.am

MOSTLYCLEANFILES = $(DX_CLEANFILES)

if DX_COND_doc

if DX_COND_html
  DOX_HTML = $(DX_DOCDIR)/html
  DOX_HTML_FILES = $(DOX_HTML)/*
  $(DOX_HTML)/* : doxygen-doc
endif

if DX_COND_pdf
  DOX_PDF = $(DX_DOCDIR)/$(DX_PROJECT).pdf
  $(DOX_PDF): doxygen-doc
endif

if DX_COND_ps
  DOX_PS = $(DX_DOCDIR)/$(DX_PROJECT).ps
  $(DOX_PS): doxygen-doc
endif

endif


dist_doc_DATA = $(DOX_HTML_FILES) $(DOX_PDF) $(DOX_PS)


EXTRA_DIST = $(pkgdata_DATA) INSTALL.configure Readseq $(DX_CONFIG) mainpage.h @SVM_SOURCE_ARCHIVE@

install-data-hook:
	if test x"$(DX_FLAG_html)" = x1; then\
	  printf "Installing html documentation";\
	  mkdir -p $(DEST_DIR)/$(docdir)/html; \
	  for i in `ls $(DX_DOCDIR)/html`; do\
	    mv $(DEST_DIR)/$(docdir)/$$i $(DEST_DIR)/$(docdir)/html/; printf ".";\
	  done;\
	  printf "done\n";\
	fi
