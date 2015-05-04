/*
  macinit.c
  -- Macintosh initializations, then call real main

Note: compile this segment as Main for generic 68000 processor, so it won't
 fail on generic mac

*/

#pragma segment Main

#include <Types.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <StdLib.h>

#include <Quickdraw.h>
#include <Memory.h>
#include <OSUtils.h>
#include <ToolUtils.h>
#include <Windows.h>
#include <Palettes.h>
#include  <dialogs.h>
#include  <StandardFile.h>
#include  <Events.h>
//#include <Menus.h>
//#include <Fonts.h>


Boolean StopKey()
{
EventRecord ev;

  if (EventAvail( keyDownMask+autoKeyMask, &ev)) {
    if (  (ev.modifiers & cmdKey)
     && ((char)(ev.message & charCodeMask) == '.') ) {
        SysBeep(1);
      (void) GetNextEvent( keyDownMask+autoKeyMask, &ev);
      return true;
      }
    }
  return false;
}

Boolean cmdKeyIsDown()
{ KeyMap  kmap;
  GetKeys(&kmap);
  return BitTst(kmap, (sizeof(KeyMap)*8) - 55);
}

Boolean shiftKeyIsDown()
{   KeyMap  kmap;
  GetKeys(&kmap);
  return BitTst(kmap, (sizeof(KeyMap)*8) - 56);
}

Boolean capsLockIsDown()
{ KeyMap  kmap;
  GetKeys(&kmap);
  return BitTst(kmap, (sizeof(KeyMap)*8) - 57);
}

Boolean optionKeyIsDown()
{ KeyMap  kmap;
  GetKeys(&kmap);
  return BitTst(kmap, (sizeof(KeyMap)*8) - 58);
}

Boolean MouseButton()
{
  return Button();
}

Boolean Keypress()
{ EventRecord ev;
  return EventAvail( keyDownMask+keyUpMask+autoKeyMask, &ev);
}



char *StdGetFile(
  char*  prompt,
  OSType fileTypes[],
  int    nFileTypes)
{
  Point       wher;       /*where to display dialog*/
  SFReply     reply;      /*reply record*/
  short       len;
  static char filename[80] = "\0";

  wher.h = 80;
  wher.v = 90;
  if (optionKeyIsDown()) nFileTypes=0;

  SFGetFile(wher, prompt, nil, nFileTypes, fileTypes, nil, &reply);

  if (reply.good) {
      len = SetVol(nil, reply.vRefNum);
    len = reply.fName[0];
    strncpy(filename, (char *)(&reply.fName[1]), len);
    filename[len]= '\0';
    return filename;
    }
  else
    return NULL;
}


int readCmdOptions(FILE *cl, char *progname, char ***argv)
/* command line reader for Mac/MPW  -- dgg */
{
#define MAXS  255
#define addarg(sptr)  if (strlen(sptr)>0) { \
  targv = (char **) realloc( targv, (argc+1) * sizeof(char *)); \
  targv[argc] = (char *) malloc(1+strlen(sptr) * sizeof(char)); \
  strcpy( targv[argc], sptr);  \
  argc++; }

  char  *pword, st[MAXS];
  int   argc = 0;
  char  **targv;

  targv = (char **) malloc(1);
  if (progname==NULL) progname = "program";
  addarg( progname);
  fgets( st, MAXS, cl);
  if (!feof(cl) && st!=NULL && *st!=0) {
    pword = strtok( st, "\ \n");
    while (pword!=NULL) {
      addarg( pword);
      pword = strtok( NULL, "\ \n");
      }
    }

  *argv = targv;
  return argc;
}

int ccommand(char ***argv)
{
  int argc;
  char  **targv;

  argc = readCmdOptions(stdin, *argv[0], &targv);
  *argv = targv;
  return argc;
}




extern _DataInit();

//#define VERSION     curSysEnvVers
#define nocolorID   130
#define no68020     133
#define no68881     132
#define no256       134
#define nosys6      135

void MacInit()
{
SysEnvRec   theWorld;
OSErr       OSys;
DialogPtr   crashDia;
long        tick;

    UnloadSeg(_DataInit);

    InitGraf((Ptr)&qd.thePort);
    //InitFonts();
    InitWindows();
    //InitMenus();
    //TEInit();
    InitDialogs(nil);
    InitCursor();

/*______________________________________________________*/
/*            If not right Machine then stop            */
/*______________________________________________________*/
  OSys = SysEnvirons( curSysEnvVers,&theWorld);

  /*if(!theWorld.hasColorQD) {
    crashDia = GetNewDialog (nocolorID, nil, (WindowPtr) -1);
    DrawDialog (crashDia);
    Delay (300, &tick);
    ExitToShell();
    }*/
  /*if(theWorld.processor < env68020) {
    crashDia = GetNewDialog (no68020, nil, (WindowPtr) -1);
    DrawDialog (crashDia);
    Delay (300, &tick);
    ExitToShell();
    }*/
  /*if(!theWorld.hasFPU) {
    crashDia = GetNewDialog (no68881, nil, (WindowPtr) -1);
    DrawDialog (crashDia);
    Delay (300, &tick);
    ExitToShell();
    }
  if(theWorld.systemVersion < 0x0600) {
    crashDia = GetNewDialog (nosys6, nil, (WindowPtr) -1);
    DrawDialog (crashDia);
    Delay (300, &tick);
    ExitToShell();
    }*/

#ifdef UnDeFineD
/*______________________________________________________*/
/*                     Set Rects                        */
/*______________________________________________________*/
  screenRect = qd.screenBits.bounds;
  offLeft = 0;
  offTop = 0;
  offRight = screenRect.right;
  offBottom = screenRect.bottom;
  SetRect(&BaseRect, 40, 60, 472, 282);
  tempRgn = GetGrayRgn();
  HLock ((Handle) tempRgn);
  TotalRect = (**tempRgn).rgnBBox;
  SetRect(&minRect, 80, 80, (**tempRgn).rgnBBox.right - 40,
        (**tempRgn).rgnBBox.bottom - 40);
  HUnlock ((Handle) tempRgn);

/*______________________________________________________*/
/*        Open Window & set Palette & Picture           */
/*______________________________________________________*/
  theGDevice = GetMainDevice();
  HLock ((Handle) theGDevice);
  mycolors = (**(**theGDevice).gdPMap).pmTable;
  numcolor = (**(**theGDevice).gdPMap).pixelSize;
  HUnlock((Handle) theGDevice);
  switch (numcolor) {
    case 1:
      numcolor = 2;
      break;
    case 2:
      numcolor = 4;
      break;
    case 4:
      numcolor = 16;
      break;
    case 8:
      numcolor = 256;
      break;
    }

  myWindow = NewCWindow(nil, &BaseRect, "", true, zoomDocProc,
              (WindowPtr) -1, true, 150);
  SetPort((WindowPtr) myWindow);
  DrawGrowIcon (myWindow);

  srcPalette = NewPalette (numcolor, mycolors, pmCourteous, 0);
  SetPalette ((WindowPtr) myWindow, srcPalette, true);

/*______________________________________________________*/
/*                    Set menus                         */
/*______________________________________________________*/
  mymenu0 = GetMenu(appleID);
  AddResMenu(mymenu0, 'DRVR');
  InsertMenu(mymenu0,0);
  mymenu1 = newmenu(129,"File");
  appendmenu(mymenu1,"Start;Quit");
  InsertMenu(mymenu1,0);
  mymenu2 = newmenu(130,"Edit");
  InsertMenu(mymenu2,0);
  DrawMenuBar();

/*______________________________________________________*/
/*                  Init variables                      */
/*______________________________________________________*/
  DoneFlag = false;
  yieldTime = 0;
  return;
#endif

}




main(int argc, char *argv[])
{
  Boolean loop = true;
  char **myargv;
  int   myargc;

  /* MacInit();  -- SIOW library handles this */
  do {
    fprintf(stderr,"\nEnter command line for %s [cmd-Q to quit]\n", argv[0]);
    fprintf(stderr,"-> %s ",argv[0]);
    myargv = argv;
    myargc = ccommand(&myargv);

    siow_main(myargc, myargv);
    fflush(stdout);

    } while (true);
  exit(0);
}

