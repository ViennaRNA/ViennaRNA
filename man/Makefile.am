manpages_verbose = $(manpages_verbose_@AM_V@)
manpages_verbose_ = $(manpages_verbose_@AM_DEFAULT_V@)
manpages_verbose_0 = @echo "  MAN      $@";

info_TEXINFOS = RNAlib.texi

main_manpages = RNA2Dfold.1 RNAaliduplex.1 RNAalifold.1 RNAcofold.1 RNAduplex.1 RNAeval.1\
                RNAfold.1 RNAheat.1 RNAinverse.1 RNALalifold.1 RNALfold.1\
                RNApaln.1 RNAparconv.1 RNAplfold.1 RNAplot.1 RNApvmin.1 RNAsubopt.1 RNAup.1\
                RNAPKplex.1 RNAsnoop.1 RNAplex.1 RNAdistance.1 RNApdist.1

dist_man_MANS = $(main_manpages) \
                ct2db.1

SUFFIXES = .1 .ggo

vpath %.ggo = ../src/bin

.ggo.1 :
	$(manpages_verbose)./help2manWrapper.sh $< > $@
##  $(eval $@_PROG_DIR=$(dir $< ))
##	$(eval $@_PROG=$(strip $(basename $(notdir $<) .ggo)))
##	$(manpages_verbose)help2man -N --help-option=--detailed-help \
##  --opt-include=./include/ref_package.inc \
##                --opt-include=./include/$($@_PROG).inc \
##                --opt-include=./include/ref_energy_par.inc \
##                "./cmdlopt.sh $($@_PROG_DIR)$($@_PROG)" > $@

ct2db.1:  ../src/Utils/ct2db.ggo
	$(manpages_verbose)help2man -N --help-option=--detailed-help \
                --opt-include=./include/ref_package.inc \
                "./cmdlopt.sh ../src/Utils/ct2db" > ct2db.1

EXTRA_DIST =  include

# remove all generated man pages
maintainer-clean-local:
	-rm -f $(dist_man_MANS)

