info_TEXINFOS = RNAlib.texinfo

man_MANS = RNAfold.man

EXTRA_DIST = $(man_MANS)

RNAfold.man:  ../Progs/RNAfold.ggo maninclude/RNAfold.examples
	./help2man -N --help-option=--detailed-help --opt-include=./include/RNAfold.inc "./cmdlopt.sh ../Progs/RNAfold.ggo" > RNAfold.man
