# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(
	      MAKEFILE =>	"Makefile.perl",
	      NAME    =>	"RNA",
	      LIBS    =>	["-lm"],
	      MYEXTLIB =>       "../lib/libRNA.a",
	      OBJECT  =>	"RNA_wrap.o",
	      VERSION =>      0.2
);
