info_TEXINFOS = RNAlib.texinfo

man_MANS = RNAfold.man RNAcofold.man

EXTRA_DIST = $(man_MANS)

RNAfold.man:  ../Progs/RNAfold.ggo include/RNAfold.inc
	./help2man -N --help-option=--detailed-help --opt-include=./include/RNAfold.inc "./cmdlopt.sh ../Progs/RNAfold.ggo" > RNAfold.man

RNAcofold.man:  ../Progs/RNAcofold.ggo
	./help2man -N --help-option=--detailed-help "./cmdlopt.sh ../Progs/RNAcofold.ggo" > RNAcofold.man
