!define VERSION "2.1.4"
!define PACKAGE "ViennaRNA Package"
!define MULTIUSER_INSTALLMODE_INSTDIR "${PACKAGE}"
BrandingText "${PACKAGE} - ${VERSION} [ TBI / Institute for Theoretical Chemistry / University of Vienna ]"

!define MULTIUSER_EXECUTIONLEVEL Highest
!define MULTIUSER_MUI
!define MULTIUSER_INSTALLMODE_COMMANDLINE
!include MultiUser.nsh
!include MUI2.nsh

!insertmacro MUI_LANGUAGE "English"

Function .onInit
  !insertmacro MULTIUSER_INIT
FunctionEnd

Function un.onInit
  !insertmacro MULTIUSER_UNINIT
FunctionEnd

# Name of the software
name "${PACKAGE}"
# Name of the installer
outfile "Install-ViennaRNA-${VERSION}_32bit.exe"

# define the directory to install to, the desktop in this case as specified  
# by the predefined $DESKTOP variable
#installDir "$PROGRAMFILES\${PACKAGE}"

#Page license
!insertmacro MUI_PAGE_LICENSE "licence.txt"

!insertmacro MULTIUSER_PAGE_INSTALLMODE
#Page directory
!insertmacro MUI_PAGE_DIRECTORY
#Page instfiles

!insertmacro MUI_PAGE_INSTFILES 

#Page components
#ComponentText bla

#Section "ViennaRNA Base" s1 
#SectionEnd
#Section Kinfold s2
#SectionEnd
#Section /o RNAforester s3
#SectionEnd
#Section Readseq s4
#SectionEnd
#Section /o "Perl interface" s5
#SectionEnd
#Section Utils s6
#SectionEnd


UninstPage uninstConfirm
UninstPage instfiles
 
# default section
section "ViennaRNA Package"
 
  # define the output path for this file
  setOutPath $INSTDIR
   
  # define what to install and place it in the output path

  # first all executable programs
  File "Progs\RNAeval.exe"
  File "Progs\RNAfold.exe"
  File "Progs\RNAheat.exe"
  File "Progs\RNApdist.exe"
  File "Progs\RNAdistance.exe"
  File "Progs\RNAinverse.exe"
  File "Progs\RNAplot.exe"
  File "Progs\RNAsubopt.exe"
  File "Progs\RNALfold.exe"
  File "Progs\RNAcofold.exe"
  File "Progs\RNAduplex.exe"
  File "Progs\RNApaln.exe"
  File "Progs\RNAalifold.exe"
  File "Progs\RNAplfold.exe"
  File "Progs\RNAup.exe"
  File "Progs\RNAaliduplex.exe"
  File "Progs\RNAparconv.exe"
  File "Progs\RNAPKplex.exe"
  File "Progs\RNALalifold.exe"
  File "Progs\RNA2Dfold.exe"
  File "Progs\RNAplex.exe"
  File "Progs\RNAsnoop.exe"
  File "Kinfold\Kinfold.exe"

  # create the Utils directory and add programs there
  CreateDirectory "$INSTDIR\Utils"
  File /oname=Utils\Fold Utils/Fold
  File /oname=Utils\dpzoom.pl Utils/dpzoom.pl
  File /oname=Utils\ct2b.pl Utils/ct2b.pl
  File /oname=Utils\b2mt.pl Utils/b2mt.pl
  File /oname=Utils\relplot.pl Utils/relplot.pl
  File /oname=Utils\mountain.pl Utils/mountain.pl
  File /oname=Utils\rotate_ss.pl Utils/rotate_ss.pl
  File /oname=Utils\colorrna.pl Utils/colorrna.pl
  File /oname=Utils\cmount.pl Utils/cmount.pl
  File /oname=Utils\coloraln.pl Utils/coloraln.pl
  File /oname=Utils\switch.pl Utils/switch.pl
  File /oname=Utils\refold.pl Utils/refold.pl
  File /oname=Utils\b2ct.exe Utils/b2ct.exe
  File /oname=Utils\popt.exe Utils/popt.exe

  # create the misc directory and add files
  CreateDirectory "$INSTDIR\Misc"
  File /oname=Misc\rna_turner2004.par misc/rna_turner2004.par
  File /oname=Misc\rna_turner1999.par misc/rna_turner1999.par
  File /oname=Misc\rna_andronescu2007.par misc/rna_andronescu2007.par
  File /oname=Misc\dna_mathews2004.par misc/dna_mathews2004.par
  File /oname=Misc\dna_mathews1999.par misc/dna_mathews1999.par
  File /oname=Misc\2Dlandscape_pf.gri misc/2Dlandscape_pf.gri
  File /oname=Misc\2Dlandscape_mfe.gri misc/2Dlandscape_mfe.gri


  # then all the necessary libraries
  File "/usr/i686-w64-mingw32/sys-root/mingw/bin/libgomp-1.dll"
  File "/usr/i686-w64-mingw32/sys-root/mingw/bin/pthreadGC2.dll"
  File "/usr/i686-w64-mingw32/sys-root/mingw/bin/libgcc_s_sjlj-1.dll"
  #File "/usr/x86_64-w64-mingw32/sys-root/mingw/bin/libgomp-1.dll"
  #File "/usr/x86_64-w64-mingw32/sys-root/mingw/bin/pthreadGC2.dll"
  File ".local"

  # we also want an uninstall to be installed
  writeUninstaller "$INSTDIR\Uninstall-${PACKAGE}.exe"

  # Start Menu
  createDirectory "$SMPROGRAMS\${PACKAGE}"
  createShortCut "$SMPROGRAMS\${PACKAGE}\Uninstall-${PACKAGE}.lnk" "$INSTDIR\Uninstall-${PACKAGE}.exe"

sectionEnd

# create a section to define what the uninstaller does.
# the section will always be named "Uninstall"
section "Uninstall"
 
  # Remove Start Menu launcher
  delete "$SMPROGRAMS\${PACKAGE}\Uninstall-${PACKAGE}.lnk"
  # Try to remove the Start Menu folder - this will only happen if it is empty
  rmDir "$SMPROGRAMS\${PACKAGE}"

  delete $INSTDIR\uninstaller.exe

  # now delete installed file
  delete $INSTDIR\RNAeval.exe
  delete $INSTDIR\RNAfold.exe
  delete $INSTDIR\RNAheat.exe
  delete $INSTDIR\RNApdist.exe
  delete $INSTDIR\RNAdistance.exe
  delete $INSTDIR\RNAinverse.exe
  delete $INSTDIR\RNAplot.exe
  delete $INSTDIR\RNAsubopt.exe
  delete $INSTDIR\RNALfold.exe
  delete $INSTDIR\RNAcofold.exe
  delete $INSTDIR\RNAduplex.exe
  delete $INSTDIR\RNApaln.exe
  delete $INSTDIR\RNAalifold.exe
  delete $INSTDIR\RNAplfold.exe
  delete $INSTDIR\RNAup.exe
  delete $INSTDIR\RNAaliduplex.exe
  delete $INSTDIR\RNAparconv.exe
  delete $INSTDIR\RNAPKplex.exe
  delete $INSTDIR\RNALalifold.exe
  delete $INSTDIR\RNA2Dfold.exe
  delete $INSTDIR\RNAplex.exe
  delete $INSTDIR\RNAsnoop.exe
  delete $INSTDIR\Kinfold.exe

  delete $INSTDIR\Utils\dpzoom.pl
  delete $INSTDIR\Utils\Fold
  delete $INSTDIR\Utils\dpzoom.pl
  delete $INSTDIR\Utils\ct2b.pl
  delete $INSTDIR\Utils\b2mt.pl
  delete $INSTDIR\Utils\relplot.pl
  delete $INSTDIR\Utils\mountain.pl
  delete $INSTDIR\Utils\rotate_ss.pl
  delete $INSTDIR\Utils\colorrna.pl
  delete $INSTDIR\Utils\cmount.pl
  delete $INSTDIR\Utils\coloraln.pl
  delete $INSTDIR\Utils\switch.pl
  delete $INSTDIR\Utils\refold.pl
  delete $INSTDIR\Utils\b2ct.exe
  delete $INSTDIR\Utils\popt.exe
  rmDir  $INSTDIR\Utils

  delete $INSTDIR\Misc\rna_turner2004.par
  delete $INSTDIR\Misc\rna_turner1999.par
  delete $INSTDIR\Misc\rna_andronescu2007.par
  delete $INSTDIR\Misc\dna_mathews2004.par
  delete $INSTDIR\Misc\dna_mathews1999.par
  delete $INSTDIR\Misc\2Dlandscape_pf.gri
  delete $INSTDIR\Misc\2Dlandscape_mfe.gri

  rmDir  $INSTDIR\Misc

  delete $INSTDIR\libgomp-1.dll
  delete $INSTDIR\pthreadGC2.dll
  delete $INSTDIR\libgcc_s_sjlj-1.dll
  delete $INSTDIR\.local

  rmDir  $INSTDIR

sectionEnd
