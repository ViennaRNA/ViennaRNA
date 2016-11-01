# Microsoft Developer Studio Project File - Name="g2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=g2 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "g2.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "g2.mak" CFG="g2 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "g2 - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "g2 - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "g2 - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "../gdwin32" /I "./src" /I "../src" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "DO_PS" /D "DO_GIF" /D "DO_WIN32" /D "DO_FIG" /YX /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "g2 - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "../gdwin32" /I "./src" /I "../src" /D "_DEBUG" /D "_WINDOWS" /D "WIN32" /D "DO_PS" /D "DO_GIF" /D "DO_WIN32" /D "DO_FIG" /FR /YX /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "g2 - Win32 Release"
# Name "g2 - Win32 Debug"
# Begin Source File

SOURCE=.\src\g2_control_pd.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_device.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_fif.c
# End Source File
# Begin Source File

SOURCE=.\src\FIG\g2_FIG.c
# End Source File
# Begin Source File

SOURCE=.\src\gd\g2_gd.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_graphic_pd.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_physical_device.c
# End Source File
# Begin Source File

SOURCE=.\src\PS\g2_PS.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_splines.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_ui_control.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_ui_device.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_ui_graphic.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_ui_virtual_device.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_util.c
# End Source File
# Begin Source File

SOURCE=.\src\g2_virtual_device.c
# End Source File
# Begin Source File

SOURCE=.\src\Win32\g2_win32.c
# End Source File
# Begin Source File

SOURCE=.\src\Win32\g2_win32_thread.c
# End Source File
# End Target
# End Project
