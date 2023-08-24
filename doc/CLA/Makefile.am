include $(top_srcdir)/silent_rules.mk

clean_cla = rm -f *.dvi *.idv *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl CLA-*.pdf

######################################################
## Tutorial related instructions below              ##
######################################################
if WITH_CLA

##--------------------------------------------------##
## prepare variables in case PDF cla           ##
## is going to be installed                         ##
##--------------------------------------------------##
if WITH_CLA_PDF
    CLA_FILES_PDF = @CLA_I_PDF_NAME@ @CLA_E_PDF_NAME@
endif WITH_CLA_PDF

##--------------------------------------------------##
## instructions on how to build the cla        ##
##--------------------------------------------------##
if WITH_CLA_BUILD

##--------------------------------------------------##
## In case the PDF cla must be generated,      ##
## here is the rule how to do so                    ##
##--------------------------------------------------##
if WITH_CLA_PDF

@CLA_I_PDF_NAME@:
	$(cla_verbose)$(clean_cla) && \
  $(CLA_CMD_LATEX) CLA-Individual >/dev/null 2>&1 && \
  $(CLA_CMD_LATEX) CLA-Individual >/dev/null 2>&1 && \
  $(CLA_CMD_LATEX) CLA-Individual >/dev/null 2>&1 && \
  mv CLA-Individual.pdf @CLA_I_PDF_NAME@ ;

@CLA_E_PDF_NAME@:
	$(cla_verbose)$(clean_cla) && \
  $(CLA_CMD_LATEX) CLA-Entity >/dev/null 2>&1 && \
  $(CLA_CMD_LATEX) CLA-Entity >/dev/null 2>&1 && \
  $(CLA_CMD_LATEX) CLA-Entity >/dev/null 2>&1 && \
  mv CLA-Entity.pdf @CLA_E_PDF_NAME@ ;

endif WITH_CLA_PDF

endif WITH_CLA_BUILD

endif WITH_CLA


##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
pdf_DATA = $(CLA_FILES_PDF)


if WITH_CLA_PDF

cla_datadir         = $(pkgdatadir)/cla


vrna_pdfcla_inst: installdirs @CLA_I_PDF_NAME@ @CLA_E_PDF_NAME@
	@$(NORMAL_INSTALL)
	test -f $(builddir)/@CLA_I_PDF_NAME@ || $(LN_S) $(srcdir)/@CLA_I_PDF_NAME@ $(builddir)/@CLA_I_PDF_NAME@
	$(INSTALL_DATA) $(builddir)/@CLA_I_PDF_NAME@ $(DESTDIR)$(pdfdir)
	{ test ! -h $(builddir)/@CLA_I_PDF_NAME@ || rm -f $(builddir)/@CLA_I_PDF_NAME@; }
	test -f $(builddir)/@CLA_E_PDF_NAME@ || $(LN_S) $(srcdir)/@CLA_E_PDF_NAME@ $(builddir)/@CLA_E_PDF_NAME@
	$(INSTALL_DATA) $(builddir)/@CLA_E_PDF_NAME@ $(DESTDIR)$(pdfdir)
	{ test ! -h $(builddir)/@CLA_E_PDF_NAME@ || rm -f $(builddir)/@CLA_E_PDF_NAME@; }

vrna_pdfcla_uninst:
	@$(NORMAL_UNINSTALL)
	rm -f "$(DESTDIR)$(pdfdir)/@CLA_I_PDF_NAME@"
	rm -f "$(DESTDIR)$(pdfdir)/@CLA_E_PDF_NAME@"

else

vrna_pdfcla_inst:
vrna_pdfcla_uninst:

endif

install-pdfDATA: vrna_pdfcla_inst
uninstall-pdfDATA: vrna_pdfcla_uninst

EXTRA_DIST =  CLA-Individual.tex \
              CLA-Entity.tex \
              @CLA_I_PDF_NAME@ \
              @CLA_E_PDF_NAME@

clean-local:
	$(clean_cla)

# remove all generated content
maintainer-clean-local:
	-rm -rf *.log *.pdf

distclean-local:
	-rm -rf *.log *.pdf
