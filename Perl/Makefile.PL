# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(	      
	      MAKEFILE =>	"Makefile.perl",
	      NAME    =>	"RNA",
	      LIBS    =>	["-lm"],
	      MYEXTLIB =>       "../lib/libRNA.a",
	      PM       =>       {"RNA.pm", '$(INST_LIBDIR)/RNA.pm'},
	      OBJECT  =>	"RNA_wrap.o",
	      AUTHOR  =>        'Ivo Hofacker <ivo@tbi.univie.ac.at>',
	      VERSION =>      0.2
);
