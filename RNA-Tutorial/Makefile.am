tutorial_verbose = $(tutorial_verbose_@AM_V@)
tutorial_verbose_ = $(tutorial_verbose_@AM_DEFAULT_V@)
tutorial_verbose_0 = @echo "  TUT      $@";

clean_tutorial = rm -f *.ps *.dvi *.idv *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl tutorial.pdf

clean_tutorial_html = rm -f *.html *.css *.4ct *.4tc *.tmp *.lg *.xref

######################################################
## Tutorial related instructions below              ##
######################################################
if WITH_TUTORIAL

##--------------------------------------------------##
## prepare variables in case HTML tutorial          ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_TUTORIAL_HTML
    TUTORIAL_FILES_HTML = html/* 
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

tutorial-html: tutorial.tex ViennaRNA-homepage.cfg
	$(tutorial_verbose)$(clean_tutorial_html); \
  $(TUTORIAL_CMD_HTLATEX) tutorial 'ViennaRNA-homepage, mathml,charset=utf8' ' -css, -cmozhtf, -cunihtf' >>tutorial_html.log 2>&1; \
  $(MKDIR_P) html; \
  ./cleanupHTML.pl -i tutorial.html > html/index.html; \
  cp tutorial.css html/; \
  $(MKDIR_P) html/Figures; \
  mv Figures/*.png html/Figures; \
  touch tutorial-html; \
  $(clean_tutorial_html);

$(TUTORIAL_FILES_HTML): tutorial-html


endif WITH_TUTORIAL_HTML

##--------------------------------------------------##
## In case the PDF tutorial must be generated,      ##
## here is the rule how to do so                    ##
##--------------------------------------------------##
if WITH_TUTORIAL_PDF

$(TUTORIAL_FILES_PDF):
	$(tutorial_verbose)$(clean_tutorial) && \
  $(TUTORIAL_CMD_LATEX) tutorial >/dev/null 2>&1 && \
  $(TUTORIAL_CMD_LATEX) tutorial >/dev/null 2>&1 && \
  $(TUTORIAL_CMD_LATEX) tutorial >/dev/null 2>&1 && \
  mv tutorial.pdf $(TUTORIAL_FILES_PDF) && \
  $(clean_tutorial);

endif WITH_TUTORIAL_PDF

endif WITH_TUTORIAL_BUILD

endif WITH_TUTORIAL


##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
html_DATA = $(TUTORIAL_FILES_HTML)
pdf_DATA = $(TUTORIAL_FILES_PDF)


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

tut_data_DATA = $(srcdir)/Data/5S.seq \
                $(srcdir)/Data/duplex.seq \
                $(srcdir)/Data/FOO.seq \
                $(srcdir)/Data/four.seq \
                $(srcdir)/Data/RNaseP.seq \
                $(srcdir)/Data/RybB.seq \
                $(srcdir)/Data/t.seq \
                $(srcdir)/Data/dot.emacs \
                $(srcdir)/Data/RNaseP.aln \
                $(srcdir)/Data/RNaseP.dnd

tut_programs_DATA = $(srcdir)/Data/programs/alifoldz_adopted.tar.gz

tut_alifoldz_DATA = $(srcdir)/Data/alifoldz/miRNA.aln \
                    $(srcdir)/Data/alifoldz/unknown.aln

tut_qrna_DATA = $(srcdir)/Data/qrna/coding.fa \
                $(srcdir)/Data/qrna/tRNA.fa

tut_rnaz_DATA = $(srcdir)/Data/rnaz/*.aln \
                $(srcdir)/Data/rnaz/*.maf

tut_rnaz_human_DATA = $(srcdir)/Data/rnaz/human/annotation.bed

tut_rnaz_humanx_DATA = $(srcdir)/Data/rnaz/human/x/*.bed \
                       $(srcdir)/Data/rnaz/human/x/*.maf \
                       $(srcdir)/Data/rnaz/human/x/*.dat \
                       $(srcdir)/Data/rnaz/human/x/*.out

tut_rnaz_herpes_DATA = $(srcdir)/Data/rnaz/herpes/annotation.bed \
                       $(srcdir)/Data/rnaz/herpes/NC_004367 \
                       $(srcdir)/Data/rnaz/herpes/NC_006146 \
                       $(srcdir)/Data/rnaz/herpes/NC_007605

tut_rnaz_herpesx_DATA = $(srcdir)/Data/rnaz/herpes/x/*.bed \
                        $(srcdir)/Data/rnaz/herpes/x/*.maf \
                        $(srcdir)/Data/rnaz/herpes/x/*.dat \
                        $(srcdir)/Data/rnaz/herpes/x/*.out


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

install-htmlDATA: vrna_htmltutorial_inst
uninstall-htmlDATA: vrna_htmltutorial_uninst

install-pdfDATA: vrna_pdftutorial_inst
uninstall-pdfDATA: vrna_pdftutorial_uninst

EXTRA_DIST =  Data \
              Figures \
              tutorial.tex \
              ViennaRNA-homepage.cfg \
              keystroke.4ht \
              Makefile.standalone \
              pgfsys-tex4ht.def \
              @TUTORIAL_PDF_NAME@

clean-local:
	$(clean_tutorial)
	$(clean_tutorial_html)

# remove all generated content
maintainer-clean-local:
	-rm -rf *.log *.pdf html tutorial-html

distclean-local:
	-rm -rf *.log *.pdf html tutorial-html
