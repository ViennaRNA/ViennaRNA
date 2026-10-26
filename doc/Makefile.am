include $(top_srcdir)/silent_rules.mk

SUBDIRS = CLA source/usage/man

pkgdata_DATA = \
    viennarna.bib

clean_latex = \
{ \
  test ! -d build/latex \
  ||  { cd build/latex \
        && rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl viennarna.pdf \
        && cd ../.. ; \
      } \
}

REFMAN_MANPAGES = \
    source/usage/man/RNA2Dfold.rst \
    source/usage/man/RNAaliduplex.rst \
    source/usage/man/RNAalifold.rst \
    source/usage/man/RNAcofold.rst \
    source/usage/man/RNAdistance.rst \
    source/usage/man/RNAdos.rst \
    source/usage/man/RNAduplex.rst \
    source/usage/man/RNAeval.rst \
    source/usage/man/RNAfold.rst \
    source/usage/man/RNAheat.rst \
    source/usage/man/RNAinverse.rst \
    source/usage/man/RNALalifold.rst \
    source/usage/man/RNALfold.rst \
    source/usage/man/RNAmultifold.rst \
    source/usage/man/RNApaln.rst \
    source/usage/man/RNAparconv.rst \
    source/usage/man/RNApdist.rst \
    source/usage/man/RNAPKplex.rst \
    source/usage/man/RNAplex.rst \
    source/usage/man/RNAplfold.rst \
    source/usage/man/RNAplot.rst \
    source/usage/man/RNApvmin.rst \
    source/usage/man/RNAsnoop.rst \
    source/usage/man/RNAsubopt.rst \
    source/usage/man/RNAup.rst

if WITH_REFERENCE_MANUAL

##--------------------------------------------------##
## prepare variables in case HTML reference manual  ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_HTML
    REFERENCE_MANUAL_FILES_HTML = build/html/*
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

refman-html: $(pkginclude_HEADERS) doxygen/doxygen-xml $(REFMAN_MANPAGES) CLA/@CLA_I_PDF_NAME@ CLA/@CLA_E_PDF_NAME@
	$(sphinx_verbose)(rm -f refman-html && \
    $(SPHINXBUILD) -M html source build >sphinx.log 2>&1 && \
    touch refman-html) || cat sphinx.log

##--------------------------------------------------##
## In case the HTML manual should be created, here  ##
## is the rule how to do so                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_HTML

$(REFERENCE_MANUAL_FILES_HTML): refman-html

endif WITH_REFERENCE_MANUAL_HTML

##--------------------------------------------------##
## In case the PDF manual should be created, here   ##
## is the rule how to do so                         ##
##--------------------------------------------------##
if WITH_REFERENCE_MANUAL_PDF

$(REFERENCE_MANUAL_FILES_PDF): doxygen/doxygen-xml
	$(sphinx_verbose)$(clean_latex); \
  $(SPHINXBUILD) -M latex source build >sphinx_latex.log 2>&1 && \
  cd build/latex; $(LATEX_CMD) -interaction=nonstopmode viennarna >>../../sphinx_latex.log 2>&1; \
  $(MAKEINDEX_CMD) -s python.ist viennarna.idx >>../../sphinx_latex.log 2>&1; \
  $(BIBTEX_CMD) viennarna >>../../sphinx_latex.log 2>&1; \
  $(LATEX_CMD) -interaction=nonstopmode viennarna >>../../sphinx_latex.log 2>&1; \
  $(LATEX_CMD) -interaction=nonstopmode viennarna >>../../sphinx_latex.log 2>&1; \
  latex_count=8 ; \
  while grep -E -s 'Rerun (LaTeX|to get cross-references right)' viennarna.log && [ $$latex_count -gt 0 ] ;\
      do \
        echo "Rerunning latex...." ;\
        $(LATEX_CMD) -interaction=nonstopmode viennarna >>../../sphinx_latex.log 2>&1;\
        latex_count=`expr $$latex_count - 1` ;\
      done; \
  $(MAKEINDEX_CMD) -s python.ist viennarna.idx >>../../sphinx_latex.log 2>&1;\
  $(LATEX_CMD) -interaction=nonstopmode viennarna >>../../sphinx_latex.log 2>&1; \
  mv viennarna.pdf "../../$(REFERENCE_MANUAL_FILES_PDF)"; \
  cd ../..; \
  $(clean_latex);

endif WITH_REFERENCE_MANUAL_PDF

endif WITH_REFERENCE_MANUAL_BUILD

endif WITH_REFERENCE_MANUAL


##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
html_DATA = $(REFERENCE_MANUAL_FILES_HTML)
pdf_DATA = $(REFERENCE_MANUAL_FILES_PDF)


if WITH_REFERENCE_MANUAL_HTML

vrna_htmlrefman_inst: installdirs refman-html
	@$(NORMAL_INSTALL)
	test -d $(builddir)/build/html || ($(MKDIR_P) $(builddir)/build && $(LN_S) $(srcdir)/build/html $(builddir)/build/html)
	@for item in $$(find $(builddir)/build/html/ -type d) ; do install -d -m 755 $(DESTDIR)$(htmldir)/$${item#$(builddir)/build/html/} ; done
	@for item in $$(find $(builddir)/build/html/ -type f) ; do install -m 644 $${item} $(DESTDIR)$(htmldir)/$${item#$(builddir)/build/html/} ; done
	{ test ! -h $(builddir)/build/html || (rm -f $(builddir)/build/html && rm -rf $(builddir)/build) }

vrna_htmlrefman_uninst:
	@$(NORMAL_UNINSTALL)
	rm -rf "$(DESTDIR)$(htmldir)"

else

vrna_htmlrefman_inst:
vrna_htmlrefman_uninst:

endif

if WITH_REFERENCE_MANUAL_PDF

vrna_pdfrefman_inst: installdirs $(REFERENCE_MANUAL_FILES_PDF)
	@$(NORMAL_INSTALL)
	test -f $(builddir)/$(REFERENCE_MANUAL_FILES_PDF) || $(LN_S) $(srcdir)/$(REFERENCE_MANUAL_FILES_PDF) $(builddir)/$(REFERENCE_MANUAL_FILES_PDF)
	$(INSTALL_DATA) $(builddir)/$(REFERENCE_MANUAL_FILES_PDF) $(DESTDIR)$(pdfdir)
	{ test ! -h $(builddir)/$(REFERENCE_MANUAL_FILES_PDF) || rm -f $(builddir)/$(REFERENCE_MANUAL_FILES_PDF); }

vrna_pdfrefman_uninst:
	@$(NORMAL_UNINSTALL)
	rm -f "$(DESTDIR)$(pdfdir)/$(REFERENCE_MANUAL_FILES_PDF)"

else

vrna_pdfrefman_inst:
vrna_pdfrefman_uninst:

endif


install-htmlDATA: vrna_htmlrefman_inst
uninstall-htmlDATA: vrna_htmlrefman_uninst

install-pdfDATA: vrna_pdfrefman_inst
uninstall-pdfDATA: vrna_pdfrefman_uninst

EXTRA_DIST =  $(pkgdata_DATA) \
              source \
              sphinx-extensions \
              man2rst.py \
              requirements.txt \
              refman-html \
              build \
              @REFERENCE_MANUAL_PDF_NAME@

clean-local:
	$(clean_latex)

# remove all generated content
maintainer-clean-local:
	-rm -rf *.log *.pdf build refman-html

distclean-local:
	-rm -rf *.log *.pdf build refman-html
