manpages_verbose = $(manpages_verbose_@AM_V@)
manpages_verbose_ = $(manpages_verbose_@AM_DEFAULT_V@)
manpages_verbose_0 = @echo "  MAN      $@";

info_TEXINFOS = RNAlib.texi

main_manpages = RNA2Dfold.1 \
                RNAaliduplex.1 \
                RNAalifold.1 \
                RNAcofold.1 \
                RNAduplex.1 \
                RNAeval.1 \
                RNAfold.1 \
                RNAheat.1 \
                RNAinverse.1 \
                RNALalifold.1 \
                RNALfold.1 \
                RNApaln.1 \
                RNAparconv.1 \
                RNAplfold.1 \
                RNAplot.1 \
                RNApvmin.1 \
                RNAsubopt.1 \
                RNAup.1 \
                RNAPKplex.1 \
                RNAsnoop.1 \
                RNAplex.1 \
                RNAdistance.1 \
                RNApdist.1 \
                RNAdos.1

dist_man_MANS = $(main_manpages) \
                ct2db.1

SUFFIXES = .1 .ggo

vpath %.ggo = ../src/bin

if VRNA_AM_SWITCH_BUILD_MANPAGES

.ggo.1 :
	$(manpages_verbose)./help2manWrapper.sh \
                "$(HELP2MAN)" $< | perl -pe 's/\\e\\e/\\e/g' > $@

ct2db.1:  ../src/Utils/ct2db.ggo
	$(manpages_verbose)$(HELP2MAN) \
                -N \
                --help-option=--detailed-help \
                --opt-include=./include/ref_package.inc \
                "./cmdlopt.sh ../src/Utils/ct2db" > ct2db.1

endif

EXTRA_DIST =  include

# remove all generated man pages
maintainer-clean-local:
	-rm -f $(dist_man_MANS)

