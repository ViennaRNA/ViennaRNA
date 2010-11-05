info_TEXINFOS = RNAlib.texinfo

dist_man_MANS = RNAfold.1 RNAcofold.1

RNAfold.1:  ../Progs/RNAfold.ggo
	./help2man -N --help-option=--detailed-help --opt-include=./include/RNAfold.inc "./cmdlopt.sh ../Progs/RNAfold.ggo" > RNAfold.1

RNAcofold.1:  ../Progs/RNAcofold.ggo
	./help2man -N --help-option=--detailed-help "./cmdlopt.sh ../Progs/RNAcofold.ggo" > RNAcofold.1
