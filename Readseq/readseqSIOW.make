#  Macintosh MPW-C Makefile
#  using Simple Input/Output Window library
#
#   File:       ReadseqSIOW.make
#   Target:     ReadseqSIOW
#   Sources:    readseq.c ureadseq.c ureadasn.c macinit.c
#   Created:    Wednesday, November 13, 1991 8:23:00 PM


#OBJECTS = macinit.c.o readseq.c.o ureadseq.c.o
#COptions =  -D SIOW  # -r

#if NCBI is available, set path here to NCBI toolkit:
NCBI = "{Boot}@molbio:ncbi:"
OBJECTS = macinit.c.o readseq.c.o ureadseq.c.o ureadasn.c.o
COptions =  -D SIOW -d NCBI -i "{NCBI}"include:  
NCBILIBS = "{NCBI}"lib:libncbi.o "{NCBI}"lib:libncbiobj.o "{NCBI}"lib:libvibrant.o
#endif NCBI

ReadseqSIOW ÄÄ ReadseqSIOW.make {OBJECTS}
	Link -d -c '????' -t APPL ¶
		{OBJECTS} ¶
		"{CLibraries}"StdClib.o ¶
		"{MPW}"Libraries:Libraries:SIOW.o ¶
		"{Libraries}"Runtime.o ¶
		"{Libraries}"Interface.o ¶
#if NCBI
		{NCBILIBS} ¶
		"{CLibraries}"CSANELib.o ¶
		"{CLibraries}"Math.o ¶
#endif NCBI
		-o ReadseqSIOW
		
readseq.c.o Ä ReadseqSIOW.make readseq.c
ureadseq.c.o Ä ReadseqSIOW.make ureadseq.c
macinit.c.o Ä ReadseqSIOW.make macinit.c
#if NCBI
ureadasn.c.o Ä ReadseqSIOW.make ureadasn.c
#endif NCBI

ReadseqSIOW ÄÄ macinit.r
	Rez -a macinit.r -o ReadseqSIOW
