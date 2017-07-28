doxygen_verbose = $(doxygen_verbose_@AM_V@)
doxygen_verbose_ = $(doxygen_verbose_@AM_DEFAULT_V@)
doxygen_verbose_0 = @echo "  DOC      $@";

tutorial_verbose = $(tutorial_verbose_@AM_V@)
tutorial_verbose_ = $(tutorial_verbose_@AM_DEFAULT_V@)
tutorial_verbose_0 = @echo "  TUT      $@";

pkgdata_DATA = viennarna.bib

clean_latex = \
{ \
  test ! -d latex \
  ||  { cd latex \
        && rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl refman.pdf doxyref.log \
        && cd .. ; \
      } \
}

clean_tutorial = \
{ \
  test ! -d tutorial \
  || { cd tutorial \
       && rm -f *.ps *.dvi *.idv *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl tutorial.pdf \
       && cd .. ; \
     } \
}

clean_tutorial_html = \
{ \
  test ! -d tutorial \
  || { cd tutorial \
       && rm -f *.html *.css *.4ct *.4tc *.tmp *.lg *.xref \
       && cd .. ; \
     } \
}

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

refman-html: $(pkginclude_HEADERS)
	$(doxygen_verbose)@doxygen@ @DOXYGEN_CONF@ >>doxyref.log 2>&1; \
  touch refman-html

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
if ! WITH_REFERENCE_MANUAL_PDFLATEX
  CONVERT_TO_PDF = $(DOXYGEN_CMD_DVIPDF) refman.dvi;
endif ! WITH_REFERENCE_MANUAL_PDFLATEX

if WITH_REFERENCE_MANUAL_PDF

$(REFERENCE_MANUAL_FILES_PDF): refman-html
	$(doxygen_verbose)$(clean_latex); \
  cd latex; $(DOXYGEN_CMD_LATEX) refman >>doxyref.log 2>&1; \
  $(DOXYGEN_CMD_MAKEINDEX) refman.idx >>doxyref.log 2>&1; \
  $(DOXYGEN_CMD_BIBTEX) refman >>doxyref.log 2>&1; \
  $(DOXYGEN_CMD_LATEX) refman >>doxyref.log 2>&1; \
  $(DOXYGEN_CMD_LATEX) refman >>doxyref.log 2>&1; \
  latex_count=8 ; \
  while egrep -s 'Rerun (LaTeX|to get cross-references right)' refman.log && [ $$latex_count -gt 0 ] ;\
      do \
        echo "Rerunning latex...." ;\
        $(DOXYGEN_CMD_LATEX) refman >>doxyref.log 2>&1;\
        latex_count=`expr $$latex_count - 1` ;\
      done; \
  $(DOXYGEN_CMD_MAKEINDEX) refman.idx >>doxyref.log 2>&1;\
  $(DOXYGEN_CMD_LATEX) refman >>doxyref.log 2>&1; \
  $(CONVERT_TO_PDF) \
  mv refman.pdf "../$(REFERENCE_MANUAL_FILES_PDF)"; \
  cd ..; \
  $(clean_latex);

endif WITH_REFERENCE_MANUAL_PDF

endif WITH_REFERENCE_MANUAL_BUILD

endif WITH_REFERENCE_MANUAL


######################################################
## Tutorial related instructions below              ##
######################################################
if WITH_TUTORIAL

##--------------------------------------------------##
## prepare variables in case HTML tutorial          ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_TUTORIAL_HTML
    TUTORIAL_FILES_HTML = tutorial/html/* 
endif WITH_TUTORIAL_HTML

##--------------------------------------------------##
## prepare variables in case PDF tutorial           ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_TUTORIAL_PDF
    TUTORIAL_FILES_PDF = @TUTORIAL_PDF_NAME@
endif WITH_TUTORIAL_PDF

##--------------------------------------------------##
## instructions on how to build the tutorial        ##
##--------------------------------------------------##
if WITH_TUTORIAL_BUILD

##--------------------------------------------------##
## In case the HTML tutorial must be generated,     ##
## here is the rule how to do so                    ##
##--------------------------------------------------##
if WITH_TUTORIAL_HTML

tutorial-html: tutorial/tutorial.tex tutorial/ViennaRNA-homepage.cfg
	$(tutorial_verbose)$(clean_tutorial_html); \
  cd tutorial; \
  $(TUTORIAL_CMD_HTLATEX) tutorial 'ViennaRNA-homepage, mathml,charset=utf8' ' -css, -cmozhtf, -cunihtf' >>tutorial_html.log 2>&1; \
  $(MKDIR_P) html; \
  ./cleanupHTML.pl -i tutorial.html > html/index.html; \
  cp tutorial.css html/; \
  $(MKDIR_P) html/Figures; \
  mv Figures/*.png html/Figures; \
  cd ..; \
  touch tutorial-html; \
  $(clean_tutorial_html);

$(TUTORIAL_FILES_HTML): tutorial-html


endif WITH_TUTORIAL_HTML

##--------------------------------------------------##
## In case the PDF tutorial must be generated,      ##
## here is the rule how to do so                    ##
##--------------------------------------------------##
if WITH_TUTORIAL_PDF

if ! WITH_TUTORIAL_PDFLATEX
  CONVERT_TUT_TO_PDF = $(TUTORIAL_CMD_DVIPDF) tutorial.dvi;
endif ! WITH_TUTORIAL_PDFLATEX

$(TUTORIAL_FILES_PDF):
	$(tutorial_verbose)$(clean_tutorial); \
  cd tutorial; \
  $(TUTORIAL_CMD_LATEX) tutorial >>tutorial.log 2>&1; \
  $(TUTORIAL_CMD_LATEX) tutorial >>tutorial.log 2>&1; \
  $(TUTORIAL_CMD_LATEX) tutorial >>tutorial.log 2>&1; \
  $(CONVERT_TUT_TO_PDF) \
  mv tutorial.pdf "../$(TUTORIAL_FILES_PDF)"; \
  cd ..; \
  $(clean_tutorial);

endif WITH_TUTORIAL_PDF

endif WITH_TUTORIAL_BUILD

endif WITH_TUTORIAL


##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
html_DATA = $(REFERENCE_MANUAL_FILES_HTML) $(TUTORIAL_FILES_HTML)
pdf_DATA = $(REFERENCE_MANUAL_FILES_PDF) $(TUTORIAL_FILES_PDF)


if WITH_REFERENCE_MANUAL_HTML

vrna_htmlrefman_inst: installdirs html
	@$(NORMAL_INSTALL)
	test -d $(builddir)/html || $(LN_S) $(srcdir)/html $(builddir)/html
	{ for i in $(builddir)/html/*; do cp -r "$$i" "$(DESTDIR)$(htmldir)/"; done }
	{ for i in "$(DESTDIR)$(htmldir)"/*; do test -f "$$i" && chmod 644 "$$i"; done }
	chmod 755 "$(DESTDIR)$(htmldir)"
	{ for i in "$(DESTDIR)$(htmldir)"/search/*; do test -f "$$i" && chmod 644 "$$i"; done }
	chmod 755 "$(DESTDIR)$(htmldir)"/search
	{ test ! -h $(builddir)/html || rm -f $(builddir)/html; }

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


vrna_htmltutorial_inst:
vrna_htmltutorial_uninst:


if WITH_TUTORIAL_PDF

tut_datadir         = $(pkgdatadir)/tutorial
tut_programsdir     = $(tut_datadir)/programs
tut_alifoldzdir     = $(tut_datadir)/alifoldz
tut_qrnadir         = $(tut_datadir)/qrna
tut_rnazdir         = $(tut_datadir)/rnaz
tut_rnaz_humandir   = $(tut_rnazdir)/human
tut_rnaz_humanxdir  = $(tut_rnaz_humandir)/x
tut_rnaz_herpesdir  = $(tut_rnazdir)/herpes
tut_rnaz_herpesxdir = $(tut_rnaz_herpesdir)/x

tut_data_DATA = $(srcdir)/tutorial/Data/5S.seq \
                $(srcdir)/tutorial/Data/duplex.seq \
                $(srcdir)/tutorial/Data/FOO.seq \
                $(srcdir)/tutorial/Data/four.seq \
                $(srcdir)/tutorial/Data/RNaseP.seq \
                $(srcdir)/tutorial/Data/RybB.seq \
                $(srcdir)/tutorial/Data/t.seq \
                $(srcdir)/tutorial/Data/dot.emacs \
                $(srcdir)/tutorial/Data/RNaseP.aln \
                $(srcdir)/tutorial/Data/RNaseP.dnd

tut_programs_DATA = $(srcdir)/tutorial/Data/programs/alifoldz_adopted.tar.gz

tut_alifoldz_DATA = $(srcdir)/tutorial/Data/alifoldz/miRNA.aln \
                    $(srcdir)/tutorial/Data/alifoldz/unknown.aln

tut_qrna_DATA = $(srcdir)/tutorial/Data/qrna/coding.fa \
                $(srcdir)/tutorial/Data/qrna/tRNA.fa

tut_rnaz_DATA = $(srcdir)/tutorial/Data/rnaz/*.aln \
                $(srcdir)/tutorial/Data/rnaz/*.maf

tut_rnaz_human_DATA = $(srcdir)/tutorial/Data/rnaz/human/annotation.bed

tut_rnaz_humanx_DATA = $(srcdir)/tutorial/Data/rnaz/human/x/*.bed \
                       $(srcdir)/tutorial/Data/rnaz/human/x/*.maf \
                       $(srcdir)/tutorial/Data/rnaz/human/x/*.dat \
                       $(srcdir)/tutorial/Data/rnaz/human/x/*.out

tut_rnaz_herpes_DATA = $(srcdir)/tutorial/Data/rnaz/herpes/annotation.bed \
                       $(srcdir)/tutorial/Data/rnaz/herpes/NC_004367 \
                       $(srcdir)/tutorial/Data/rnaz/herpes/NC_006146 \
                       $(srcdir)/tutorial/Data/rnaz/herpes/NC_007605

tut_rnaz_herpesx_DATA = $(srcdir)/tutorial/Data/rnaz/herpes/x/*.bed \
                        $(srcdir)/tutorial/Data/rnaz/herpes/x/*.maf \
                        $(srcdir)/tutorial/Data/rnaz/herpes/x/*.dat \
                        $(srcdir)/tutorial/Data/rnaz/herpes/x/*.out


vrna_pdftutorial_inst: installdirs $(TUTORIAL_FILES_PDF)
	@$(NORMAL_INSTALL)
	test -f $(builddir)/$(TUTORIAL_FILES_PDF) || $(LN_S) $(srcdir)/$(TUTORIAL_FILES_PDF) $(builddir)/$(TUTORIAL_FILES_PDF)
	$(INSTALL_DATA) $(builddir)/$(TUTORIAL_FILES_PDF) $(DESTDIR)$(pdfdir)
	{ test ! -h $(builddir)/$(TUTORIAL_FILES_PDF) || rm -f $(builddir)/$(TUTORIAL_FILES_PDF); }

vrna_pdftutorial_uninst:
	@$(NORMAL_UNINSTALL)
	rm -f "$(DESTDIR)$(pdfdir)/$(TUTORIAL_FILES_PDF)"

else

vrna_pdftutorial_inst:
vrna_pdftutorial_uninst:

endif

install-htmlDATA: vrna_htmlrefman_inst vrna_htmltutorial_inst
uninstall-htmlDATA: vrna_htmlrefman_uninst vrna_htmltutorial_uninst

install-pdfDATA: vrna_pdfrefman_inst vrna_pdftutorial_inst
uninstall-pdfDATA: vrna_pdfrefman_uninst vrna_pdftutorial_uninst

EXTRA_DIST =  $(pkgdata_DATA) \
              refman.include \
              refman.images \
              refman-html \
              html \
              latex \
              @REFERENCE_MANUAL_PDF_NAME@ \
              tutorial \
              @TUTORIAL_PDF_NAME@

clean-local:
	$(clean_latex)
	$(clean_tutorial)
	$(clean_tutorial_html)

# remove all generated content
maintainer-clean-local:
	-rm -rf *.log *.pdf html latex doxygen.warnings refman-html tutorial/html tutorial-html

distclean-local:
	-rm -rf *.log *.pdf html latex doxygen.warnings refman-html tutorial/html tutorial-html
