Howto build perl library

perl Makefile.PL
swig -perl5 -o g2_wrap.c -I../src -DDO_PS -DDO_GD ../src/g2.h
"c:\Program Files\SWIG-1.3.13\swig" -perl5 -o g2_wrap.c -I../src -DDO_PS -DDO_GD -DDO_WIN32 ../src/g2.h
nmake
nmake install

To make ActiveState PPD:

tar cvf g2.tar blib
gzip --best g2.tar
nmake ppd
copy tar.gz file into subdir x86
You have to edit the resulting PPD file and add the location of the package archive into <CODEBASE HREF="x86/g2.tar.gz" />. The location is relative to the PPD file.
To install do:
ppm install g2.ppd

