# File : Makefile.PL
# Last changed Time-stamp: <1999-10-27 12:32:33 ivo>
use ExtUtils::MakeMaker;
WriteMakefile(
	      NAME     =>	"RNA",
	      LIBS     =>	["-lm"],
	      MYEXTLIB =>       "../lib/libRNA.a",
	      OBJECT   =>	"RNA_wrap.o",
	      VERSION  =>	0.2
	      );
