# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(
	NAME	=>	"RNA",
	LIBS	=>	["-L/usr/local/lib -lRNA"],
	OBJECT  =>	"RNA_wrap.o"
);
