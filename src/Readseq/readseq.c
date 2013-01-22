/* File: readseq.c
 * main() program for ureadseq.c, ureadseq.h
 *
 * Reads and writes nucleic/protein sequence in various
 * formats. Data files may have multiple sequences.
 *
 * Copyright 1990 by d.g.gilbert
 * biology dept., indiana university, bloomington, in 47405
 * e-mail: gilbertd@bio.indiana.edu
 *
 * This program may be freely copied and used by anyone.
 * Developers are encourged to incorporate parts in their
 * programs, rather than devise their own private sequence
 * format.
 *
 * This should compile and run with any ANSI C compiler.
 * Please advise me of any bugs, additions or corrections.
 *
 */

const char *title
    = "readSeq (1Feb93), multi-format molbio sequence reader.\n";

 /*  History
  27 Feb 90.  1st release to public.
   4 Mar 90.  + Gary Olsen format
              + case change
              * minor corrections to NBRF,EMBL,others
              * output 1 file per sequence for gcg, unknown
              * define -DNOSTR for c-libraries w/o strstr
              - readseq.p, pascal version, becomes out-of-date
  24 May 90.  + Phylip 3.2 output format (no input)
  20 Jul 90.  + Phylip 3.3 output (no input yet)
              + interactive output re-direction
              + verbose progress info
              * interactive help output
              * dropped line no.s on NBRF output
              * patched in HyperGCG XCMD corrections,
                - except for seq. documentation handling
              * dropped the IG special nuc codes, as IG has
                adopted the standard IUB codes (now if only
                everyone would adopt a standard format !)
  11 Oct 90.  * corrected bug in reading/writing of EMBL format

  17 Oct 91.  * corrected bug in reading Olsen format
                (serious-deletion)
  10 Nov 91.  * corrected bug in reading some GCG format files
                (serious-last line duplicated)
              + add format name parsing (-fgb, -ffasta, ...)
              + Phylip v3.4 output format (== v3.2, sequential)
              + add checksum output to all forms that have document
              + skip mail headers in seq file
              + add pipe for standard input == seq file (with -p)
              * fold in parts of MacApp Seq object
              * strengthen format detection
              * clarify program structure
              * remove fixed sequence size limit (now dynamic, sizeof memory)
              * check and fold in accumulated bug reports:
              *   Now ANSI-C fopen(..,"w") & check open failure
              *   Define -DFIXTOUPPER for nonANSI C libraries that mess
                  up toupper/tolower
              = No command-line changes; callers of readseq main() should be okay
              - ureadseq.h functions have changed; client programs need to note.
              + added Unix and VMS Make scripts, including validation tests

   4 May 92.  + added 32 bit CRC checksum as alternative to GCG 6.5bit checksum
                (-DBIGCHECKSUM)
    Aug 92    = fixed Olsen format input to handle files w/ more sequences,
                not to mess up when more than one seq has same identifier,
                and to convert number masks to symbols.
              = IG format fix to understand ^L

  25-30 Dec 92
              * revised command-line & interactive interface.  Suggested form is now
                  readseq infile -format=genbank -output=outfile -item=1,3,4 ...
                but remains compatible with prior commandlines:
                  readseq infile -f2 -ooutfile -i3 ...
              + added GCG MSF multi sequence file format
              + added PIR/CODATA format
              + added NCBI ASN.1 sequence file format
              + added Pretty, multi sequence pretty output (only)
              + added PAUP multi seq format
              + added degap option
              + added Gary Williams (GWW, G.Williams@CRC.AC.UK) reverse-complement option.
              + added support for reading Phylip formats (interleave & sequential)
              * string fixes, dropped need for compiler flags NOSTR, FIXTOUPPER, NEEDSTRCASECMP
              * changed 32bit checksum to default, -DSMALLCHECKSUM for GCG version

   1Feb93
              = revert GenBank output to a fixed left number width which 
               other software depends on.
	      = fix for MSF input to handle symbols in names
	      = fix bug for possible memory overrun when truncating seqs for
		Phylip or Paup formats (thanks Anthony Persechini)

 */



/*
   Readseq has been tested with:
      Macintosh MPW C
      GNU gcc
      SGI cc
      VAX-VMS cc
   Any ANSI C compiler should be able to handle this.
   Old-style C compilers barf all over the source.


How do I build the readseq program if I have an Ansi C compiler?
#--------------------
# Unix ANSI C
# Use the supplied Makefile this way:
%  make CC=name-of-c-compiler
# OR do this...
% gcc readseq.c ureadseq.c -o readseq

#--------------------
$!VAX-VMS cc
$! Use the supplied Make.Com this way:
$  @make
$! OR, do this:
$ cc readseq, ureadseq
$ link readseq, ureadseq, sys$library:vaxcrtl/lib
$ readseq :== $ MyDisk:[myacct]readseq

#--------------------
# Macintosh Simple Input/Output Window application
# requires MPW-C and SIOW library (from APDA)
# also uses files macinit.c, macinit.r, readseqSIOW.make
#
Buildprogram readseqSIOW

#--------------------
#MPW-C v3 tool
C  ureadseq.c
C  readseq.c
link -w -o readseq -t MPST -c 'MPS ' ¶
   readseq.c.o Ureadseq.c.o ¶
    "{Libraries}"Interface.o ¶
    "{Libraries}"ToolLibs.o ¶
    "{Libraries}"Runtime.o ¶
    "{CLibraries}"StdClib.o
readseq -i1 ig.seq

# MPW-C with NCBI tools

set NCBI "{Boot}@molbio:ncbi:"; EXPORT NCBI
set NCBILIB1  "{NCBI}"lib:libncbi.o; export NCBILIB1
set NCBILIB2  "{NCBI}"lib:libncbiobj.o; export NCBILIB2
set NCBILIB3  "{NCBI}"lib:libncbicdr.o; export NCBILIB3
set NCBILIB4  "{NCBI}"lib:libvibrant.o; export NCBILIB4

C  ureadseq.c
C  -d NCBI -i "{NCBI}"include: ureadasn.c
C  -d NCBI -i "{NCBI}"include: readseq.c
link -w -o readseq -t MPST -c 'MPS ' ¶
   ureadseq.c.o ureadasn.c.o readseq.c.o  ¶
    {NCBILIB4} {NCBILIB2} {NCBILIB1} ¶
    "{Libraries}"Interface.o ¶
    "{Libraries}"ToolLibs.o ¶
    "{Libraries}"Runtime.o ¶
    "{CLibraries}"CSANELib.o ¶
    "{CLibraries}"Math.o ¶
    "{CLibraries}"StdClib.o

===========================================================*/



#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "ureadseq.h"

#pragma segment readseq



static char inputfilestore[256], *inputfile = inputfilestore;

const char *formats[kMaxFormat+1] = {
    " 1. IG/Stanford",
    " 2. GenBank/GB",
    " 3. NBRF",
    " 4. EMBL",
    " 5. GCG",
    " 6. DNAStrider",
    " 7. Fitch",
    " 8. Pearson/Fasta",
    " 9. Zuker (in-only)",
    "10. Olsen (in-only)",
    "11. Phylip3.2",
    "12. Phylip",
    "13. Plain/Raw",
    "14. PIR/CODATA",
    "15. MSF",
    "16. ASN.1",
    "17. PAUP/NEXUS",
    "18. Pretty (out-only)",
    "" };

#define kFormCount  30
#define kMaxFormName 15

const  struct formatTable {
  char  *name;
  short num;
  } formname[] = {
    {"ig",  kIG},
    {"stanford", kIG},
    {"genbank", kGenBank},
    {"gb", kGenBank},
    {"nbrf", kNBRF},
    {"embl", kEMBL},
    {"gcg", kGCG},
    {"uwgcg", kGCG},
    {"dnastrider", kStrider},
    {"strider", kStrider},
    {"fitch", kFitch},
    {"pearson", kPearson},
    {"fasta", kPearson},
    {"zuker", kZuker},
    {"olsen", kOlsen},
    {"phylip", kPhylip},
    {"phylip3.2", kPhylip2},
    {"phylip3.3", kPhylip3},
    {"phylip3.4", kPhylip4},
    {"phylip-interleaved", kPhylip4},
    {"phylip-sequential", kPhylip2},
    {"plain", kPlain},
    {"raw", kPlain},
    {"pir", kPIR},
    {"codata", kPIR},
    {"asn.1", kASN1},
    {"msf", kMSF},
    {"paup", kPAUP},
    {"nexus", kPAUP},
    {"pretty", kPretty},
  };

const char *kASN1headline = "Bioseq-set ::= {\nseq-set {\n";

/* GWW table for getting the complement of a nucleotide (IUB codes) */
/*                     ! "#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \]^_`abcdefghijklmnopqrstuvwxyz{|}~ */
const char compl[] = " !\"#$%&'()*+,-./0123456789:;<=>?@TVGHNNCDNNMNKNNYRYSAABWNRN[\\]^_`tvghnncdnnmnknnyrysaabwnrn{|}~";



char *formatstr( short format)
{
  if (format < 1 || format > kMaxFormat) {
    switch (format) {
      case kASNseqentry :
      case kASNseqset   : return formats[kASN1-1];
      case kPhylipInterleave:
      case kPhylipSequential: return formats[kPhylip-1];
      default: return "(unknown)";
      }
    }
  else return formats[format-1];
}

int parseformat( char *name)
{
#define kDupmatch  -2
  int   namelen, maxlen, i, match, matchat;
  char  lname[kMaxFormName+1];

  skipwhitespace(name);
  namelen = strlen(name);
  if (namelen == 0)
    return kNoformat;
  else if (isdigit(*name)) {
    i = atol( name);
    if (i < kMinFormat | i > kMaxFormat) return kNoformat;
    else return i;
    }

  /* else match character name */
  maxlen = min( kMaxFormName, namelen);
  for (i=0; i<maxlen; i++) lname[i] = to_lower(name[i]);
  lname[maxlen]=0;
  matchat = kNoformat;

  for (i=0; i<kFormCount; i++) {
    match = strncmp( lname, formname[i].name, maxlen);
    if (match == 0) {
      if (strlen(formname[i].name) == namelen) return (formname[i].num);
      else if (matchat == kNoformat) matchat = i;
      else matchat = kDupmatch; /* 2 or more partial matches */
      }
    }
  if (matchat == kNoformat || matchat == kDupmatch)
    return kNoformat;
  else
    return formname[matchat].num;
}



static void dumpSeqList(char *list, short format)
{
  long i, l, listlen;
  char s[256];

  listlen = strlen(list);
  printf("Sequences in %s  (format is %s)\n", inputfile, formatstr(format));
  for (i=0, l=0; i < listlen; i++) {
    if (list[i] == (char)NEWLINE) {
      s[l] = '\0'; l = 0;
      puts(s);
      }
    else if (l < 255)
      s[l++] = list[i];
    }
  putchar('\n');
}



void usage()
{
  short   i, midi;

  fprintf(stderr,title);
  fprintf(stderr,
  "usage: readseq [-options] in.seq > out.seq\n");
  fprintf(stderr," options\n");
/* ? add -d[igits] to allow digits in sequence data, &/or option to specify seq charset !? */
  fprintf(stderr, "    -a[ll]         select All sequences\n");
  fprintf(stderr, "    -c[aselower]   change to lower case\n");
  fprintf(stderr, "    -C[ASEUPPER]   change to UPPER CASE\n");
  fprintf(stderr, "    -degap[=-]     remove gap symbols\n");
  fprintf(stderr, "    -i[tem=2,3,4]  select Item number(s) from several\n");
  fprintf(stderr, "    -l[ist]        List sequences only\n");
  fprintf(stderr, "    -o[utput=]out.seq  redirect Output\n");
  fprintf(stderr, "    -p[ipe]        Pipe (command line, <stdin, >stdout)\n");
  fprintf(stderr, "    -r[everse]     change to Reverse-complement\n");
  fprintf(stderr, "    -v[erbose]     Verbose progress\n");
  fprintf(stderr, "    -f[ormat=]#    Format number for output,  or\n");
  fprintf(stderr, "    -f[ormat=]Name Format name for output:\n");
  midi = (kMaxFormat+1) / 2;
  for (i = kMinFormat-1; i < midi; i++)
   fprintf( stderr, "        %-20s      %-20s\n",
    formats[i], formats[midi+i]);

  /* new output format options, esp. for pretty format: */
  fprintf(stderr, "     \n");
  fprintf(stderr, "   Pretty format options: \n");
  fprintf(stderr, "    -wid[th]=#            sequence line width\n");
  fprintf(stderr, "    -tab=#                left indent\n");
  fprintf(stderr, "    -col[space]=#         column space within sequence line on output\n");
  fprintf(stderr, "    -gap[count]           count gap chars in sequence numbers\n");
  fprintf(stderr, "    -nameleft, -nameright[=#]   name on left/right side [=max width]\n");
  fprintf(stderr, "    -nametop              name at top/bottom\n");
  fprintf(stderr, "    -numleft, -numright   seq index on left/right side\n");
  fprintf(stderr, "    -numtop, -numbot      index on top/bottom\n");
  fprintf(stderr, "    -match[=.]            use match base for 2..n species\n");
  fprintf(stderr, "    -inter[line=#]        blank line(s) between sequence blocks\n");

  /******  not ready yet
  fprintf(stderr, "    -code=none,rtf,postscript,ps   code syntax\n");
  fprintf(stderr, "    -namefont=, -numfont=, -seqfont=font   font choice\n");
  fprintf(stderr, "       font suggestions include times,courier,helvetica\n");
  fprintf(stderr, "    -namefontsize=, -numfontsize=, -seqfontsize=#\n");
  fprintf(stderr, "       fontsize suggestions include 9,10,12,14\n");
  fprintf(stderr, "    -namefontstyle=, -numfontstyle=, -seqfontstyle= style  fontstyle for names\n");
  fprintf(stderr, "       fontstyle options are plain,italic,bold,bold-italic\n");
  ******/
}

void erralert(short err)
{
  switch (err) {
    case 0  :
      break;
    case eFileNotFound: fprintf(stderr, "File not found: %s\n", inputfile);
      break;
    case eFileCreate: fprintf(stderr, "Can't open output file.\n");
      break;
    case eASNerr: fprintf(stderr, "Error in ASN.1 sequence routines.\n");
      break;
    case eNoData: fprintf(stderr, "No data in file.\n");
      break;
    case eItemNotFound: fprintf(stderr, "Specified item not in file.\n");
      break;
    case eUnequalSize:  fprintf(stderr,
      "This format requires equal length sequences.\nSequence truncated or padded to fit.\n");
      break;
    case eUnknownFormat: fprintf(stderr, "Error: this format is unknown to me.\n");
      break;
    case eOneFormat: fprintf(stderr,
      "Warning: This format permits only 1 sequence per file.\n");
      break;
    case eMemFull: fprintf(stderr, "Out of storage memory. Sequence truncated.\n");
      break;
    default: fprintf(stderr, "readSeq error = %d\n", err);
      break;
    }
} /* erralert */


int chooseFormat( boolean quietly)
{
  char  sform[128];
  int   midi, i, outform;

    if (quietly)
      return kPearson;  /* default */
    else {
      midi = (kMaxFormat+1) / 2;
      for (i = kMinFormat-1; i < midi; i++)
        fprintf( stderr, "        %-20s      %-20s\n",
                        formats[i], formats[midi+i]);
      fprintf(stderr,"\nChoose an output format (name or #): \n");
      gets(sform);
      outform = parseformat(sform);
      if (outform == kNoformat) outform = kPearson;
      return outform;
      }
}



/* read paramater(s) */

boolean checkopt( boolean casesense, char *sopt, const char *smatch, short minword)
{
  long  lenopt, lenmatch;
  boolean result;
  short minmaxw;

  lenopt = strlen(sopt);
  lenmatch= strlen(smatch);
  minmaxw= max(minword, min(lenopt, lenmatch));

  if (casesense)
    result= (!strncmp( sopt, smatch, minmaxw));
  else
    result= (!Strncasecmp( sopt, smatch, minmaxw ));
  /* if (result) { */
    /* fprintf(stderr,"true checkopt(opt=%s,match=%s,param=%s)\n", sopt, smatch, *sparam); */
  /*  } */
  return result;
}


#define   kMaxwhichlist  50

/* global for readopt(), main() */
boolean   chooseall = false, quietly = false, gotinputfile = false,
          listonly = false, closeout = false, verbose = false,
          manyout = false, dolower = false, doupper = false, doreverse= false,
          askout  = true, dopipe= false, interleaved = false;
short     nfile = 0, iwhichlist=0, nwhichlist = 0;
short     whichlist[kMaxwhichlist+1];
long      whichSeq = 0, outform = kNoformat;
char      onamestore[128], *oname = onamestore;
FILE      *foo = NULL;

void resetGlobals()
/* need this when used from SIOW, as these globals are not reinited automatically
between calls to local main() */
{
  chooseall = false; quietly = false; gotinputfile = false;
  listonly = false; closeout = false; verbose = false;
  manyout = false; dolower = false; doupper = false; doreverse= false;
  askout  = true; dopipe= false; interleaved = false;
  nfile = 0; iwhichlist=0; nwhichlist = 0;
  whichSeq = 0; outform = kNoformat;
  oname = onamestore;
  foo = NULL;

  gPrettyInit(gPretty);
}


#define kOptOkay  1
#define kOptNone  0

int readopt( char *sopt)
{
  char    sparamstore[256], *sparam= sparamstore;
  short   n, slen= strlen(sopt);

  /* fprintf(stderr,"readopt( %s) == ", sopt); */

  if (*sopt == '?') {
    usage();
    return kOptNone;   /*? eOptionBad or kOptNone */
    }

  else if (*sopt == '-') {

    char *cp= strchr(sopt,'=');
    *sparam= '\0';
    if (cp) {
      strcpy(sparam, cp+1);
      *cp= 0;
      }

    if (checkopt( false, sopt, "-help", 2)) {
      usage();
      return kOptNone;
      }

    if (checkopt( false, sopt, "-all", 2)) {
      whichSeq= 1; chooseall= true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-colspace", 4)) { /* test before -c[ase] */
      n= atoi( sparam);
      gPretty.spacer = n;
      return kOptOkay;
      }

    if (checkopt( true, sopt, "-caselower", 2)) {
      dolower= true;
      return kOptOkay;
      }
    if (checkopt( true, sopt, "-CASEUPPER", 2)) {
      doupper= true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-pipe", 2)) {
      dopipe= true; askout= false;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-list", 2)) {
      listonly = true; askout = false;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-reverse", 2)) {
      doreverse = true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-verbose", 2)) {
      verbose = true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-match", 5)) {
      gPretty.domatch= true;
      if (*sparam >= ' ') gPretty.matchchar= *sparam;
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-degap", 4)) {
      gPretty.degap= true;
      if (*sparam >= ' ') gPretty.gapchar= *sparam;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-interline", 4)) {
      gPretty.interline= atoi( sparam);
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-item", 2)) {
      char  *cp = sparam;
      nwhichlist= 0;
      whichlist[0]= 0;
      if (*cp == 0) cp= sopt+2; /* compatible w/ old way */
      do {
        while (*cp!=0 && !isdigit(*cp)) cp++;
        if (*cp!=0) {
          n = atoi( cp);
          whichlist[nwhichlist++]= n;
          while (*cp!=0 && isdigit(*cp)) cp++;
          }
      } while (*cp!=0 && n>0 && nwhichlist<kMaxwhichlist);
      whichlist[nwhichlist++]= 0; /* 0 == stopsign for loop */
      whichSeq= max(1,whichlist[0]); iwhichlist= 1;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-format", 5)) {/* -format=phylip, -f2, -form=phylip */
      if (*sparam==0) { for (sparam= sopt+2; isalpha(*sparam); sparam++) ; }
      outform = parseformat( sparam);
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-f", 2)) { /* compatible w/ -fphylip prior version */
      if (*sparam==0) sparam= sopt+2;
      outform = parseformat( sparam);
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-output", 3)) {/* -output=myseq */
      if (*sparam==0) { for (sparam= sopt+3; isalpha(*sparam); sparam++) ; }
      strcpy( oname, sparam);
      foo = fopen( oname, "w");
      if (!foo) { erralert(eFileCreate); return eFileCreate; }
      closeout = true;
      askout = false;
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-o", 2)) {  /* compatible w/ -omyseq prior version */
      if (*sparam==0) sparam= sopt+2;
      strcpy( oname, sparam);
      foo = fopen( oname, "w");
      if (!foo) { erralert(eFileCreate); return eFileCreate; }
      closeout = true;
      askout = false;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-width", 2)) {
      if (*sparam==0) { for (sparam= sopt+2; !isdigit(*sparam) && *sparam!=0; sparam++) ; }
      n= atoi( sparam);
      if (n>0) gPretty.seqwidth = n;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-tab", 4)) {
      if (*sparam==0) { for (sparam= sopt+2; !isdigit(*sparam) && *sparam!=0; sparam++) ; }
      n= atoi( sparam);
      gPretty.tab = n;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-gapcount", 4)) {
      gPretty.baseonlynum = false;
      /* if (*sparam >= ' ') gPretty.gapchar= *sparam; */
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-nointerleave", 8)) {
      gPretty.noleaves = true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-nameleft", 7)) {
      if (*sparam==0) { for (sparam= sopt+2; !isdigit(*sparam) && *sparam!=0; sparam++) ; }
      n= atoi( sparam);
      if (n>0 && n<50) gPretty.namewidth =  n;
      gPretty.nameleft= true;
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-nameright", 7)) {
      if (*sparam==0) { for (sparam= sopt+2; !isdigit(*sparam) && *sparam!=0; sparam++) ; }
      n= atoi( sparam);
      if (n>0 && n<50) gPretty.namewidth =  n;
      gPretty.nameright= true;
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-nametop", 6)) {
      gPretty.nametop= true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-numleft", 6)) {
      if (*sparam==0) { for (sparam= sopt+2; !isdigit(*sparam) && *sparam!=0; sparam++) ; }
      n= atoi( sparam);
      if (n>0 && n<50) gPretty.numwidth =  n;
      gPretty.numleft= true;
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-numright", 6)) {
      if (*sparam==0) { for (sparam= sopt+2; !isdigit(*sparam) && *sparam!=0; sparam++) ; }
      n= atoi( sparam);
      if (n>0 && n<50) gPretty.numwidth =  n;
      gPretty.numright= true;
      return kOptOkay;
      }

    if (checkopt( false, sopt, "-numtop", 6)) {
      gPretty.numtop= true;
      return kOptOkay;
      }
    if (checkopt( false, sopt, "-numbottom", 6)) {
      gPretty.numbot= true;
      return kOptOkay;
      }

    else {
      usage();
      return eOptionBad;
      }
    }

  else {
    strcpy( inputfile, sopt);
    gotinputfile = (*inputfile != 0);
    nfile++;
    return kOptOkay;
    }

 /* return kOptNone; -- never here */
}




/* this program suffers some as it tries to be a quiet translator pipe
   _and_ a noisy user interactor
*/

/* return is best for SIOW, okay for others */
#ifdef SIOW
#define Exit(a)   return(a)
siow_main( int argc, char *argv[])

#else
#define Exit(a)   exit(a)

main( int argc, char *argv[])
#endif
{
boolean   closein = false;
short     ifile, nseq, atseq, format, err = 0, seqtype = kDNA,
          nlines, seqout = 0, phylvers = 2;
long      i, skiplines, seqlen, seqlen0;
unsigned long  checksum= 0, checkall= 0;
char      *seq, *cp, *firstseq = NULL, *seqlist, *progname, tempname[256];
char      seqid[256], *seqidptr = seqid;
char      stempstore[256], *stemp = stempstore;
FILE      *ftmp, *fin, *fout;
long      outindexmax= 0, noutindex= 0, *outindex = NULL;

#define exit_main(err) {        \
  if (closeout) fclose(fout);   \
  if (closein) fclose(fin);   \
  if (*tempname!=0) remove(tempname);\
  Exit(err); }

#define indexout()  if (interleaved) {\
  if (noutindex>=outindexmax) {\
    outindexmax= noutindex + 20;\
    outindex= (long*) realloc(outindex, sizeof(long)*outindexmax);\
    if (outindex==NULL) { err= eMemFull; erralert(err); exit_main(err); }\
    }\
  outindex[noutindex++]= ftell(fout);\
  }


  resetGlobals();
  foo = stdout;
  progname = argv[0];
  *oname = 0;
  *tempname = 0;
  /* initialize gPretty ?? -- done in header */

  for (i=1; i < argc; i++) {
    err= readopt( argv[i]);
    if (err <= 0) exit_main(err);
    }

                            /* pipe input from stdin !? */
  if (dopipe && !gotinputfile) {
    int c;
    tmpnam(tempname);
    inputfile = tempname;
    ftmp = fopen( inputfile, "w");
    if (!ftmp) { erralert(eFileCreate); exit_main(eFileCreate); }
    while ((c = getc(stdin)) != EOF) fputc(c, ftmp);
    fclose(ftmp);
    gotinputfile= true;
    }

  quietly = (dopipe || (gotinputfile && (listonly || whichSeq != 0)));

  if (verbose || (!quietly && !gotinputfile)) fprintf( stderr, title);
  ifile = 1;

                            /* UI: Choose output */
  if (askout && !closeout && !quietly) {
    askout = false;
    fprintf(stderr,"\nName of output file (?=help, defaults to display): \n");
    gets(oname= onamestore);
    skipwhitespace(oname);
    if (*oname == '?') { usage(); exit_main(0); }
    else if (*oname != 0) {
      closeout = true;
      foo = fopen( oname, "w");
      if (!foo) { erralert(eFileCreate); exit_main(eFileCreate); }
      }
    }

  fout = foo;
  if (outform == kNoformat) outform = chooseFormat(quietly);

                          /* set up formats ... */
  switch (outform) {
    case kPhylip2:
      interleaved= false;
      phylvers = 2;
      outform = kPhylip;
      break;

    case kPhylip4:
      interleaved= true;
      phylvers = 4;
      outform = kPhylip;
      break;

    case kMSF:
    case kPAUP:
      interleaved= true;
      break;

    case kPretty:
      gPretty.isactive= true;
      interleaved= true;
      break;

    }

  if (gPretty.isactive && gPretty.noleaves) interleaved= false;
  if (interleaved) {
    fout = ftmp = tmpfile();
    outindexmax= 30; noutindex= 0;
    outindex = (long*) malloc(outindexmax*sizeof(long));
    if (outindex==NULL) { err= eMemFull; erralert(err); exit_main(err); }
    }

                        /* big loop over all input files */
  do {
                        /* select next input file */
    gotinputfile = (*tempname != 0);
    while ((ifile < argc) && (!gotinputfile)) {
      if (*argv[ifile] != '-') {
        strcpy( inputfile, argv[ifile]);
        gotinputfile = (*inputfile != 0);
        --nfile;
        }
      ifile++;
      }

    while (!gotinputfile) {
      fprintf(stderr,"\nName an input sequence or -option: \n");
      inputfile= inputfilestore;

      gets(stemp= stempstore);
      if (*stemp==0) goto fini;  /* !! need this to finish work during interactive use */
      stemp= strtok(stempstore, " \n\r\t");
      while (stemp) {
        err= readopt( stemp); /* will read inputfile if it exists */
        if (err<0) exit_main(err);
        stemp= strtok( NULL, " \n\r\t");
        }
      }
              /* thanks to AJB@UK.AC.DARESBURY.DLVH for this PHYLIP3 fix: */
              /* head for end (interleave if needed) */
    if (*inputfile == 0) break;

    format = seqFileFormat( inputfile, &skiplines, &err);

    if (err == 0)  {
#ifdef NCBI
      if (format == kASNseqentry || format == kASNseqset)
        seqlist = listASNSeqs( inputfile, skiplines, format, &nseq, &err);
      else
#endif
        seqlist = listSeqs( inputfile, skiplines, format, &nseq, &err);
      }

    if (err != 0)
      erralert(err);

    else if (listonly) {
      dumpSeqList(seqlist,format);
      free( seqlist);
      }

    else {
                                /* choose whichSeq if needed */
      if (nseq == 1 || chooseall || (quietly && whichSeq == 0)) {
        chooseall= true;
        whichSeq = 1;
        quietly = true; /* no loop */
        }
      else if (whichSeq > nseq && quietly) {
        erralert(eItemNotFound);
        err= eItemNotFound;
        }
      else if (whichSeq > nseq || !quietly) {
        dumpSeqList(seqlist, format);
        fprintf(stderr,"\nChoose a sequence (# or All): \n");
        gets(stemp= stempstore);
        skipwhitespace(stemp);
        if (to_lower(*stemp) == 'a') {
          chooseall= true;
          whichSeq = 1;
          quietly = true; /* !? this means we don't ask for another file 
                            as well as no more whichSeqs... */
          }
        else if (isdigit(*stemp)) whichSeq= atol(stemp);
        else whichSeq= 1; /* default */
        }
      free( seqlist);

      if (false /*chooseall*/) {  /* this isn't debugged yet...*/
        fin = fopen(inputfile, "r");
        closein= true;
        }

      while (whichSeq > 0 && whichSeq <= nseq) {
                                /* need to open multiple output files ? */
        manyout = ((chooseall || nwhichlist>1) && nseq > 1
                  && (outform == kPlain || outform == kGCG));
        if (manyout) {
          if ( whichSeq == 1 ) erralert(eOneFormat);
          else if (closeout) {
            sprintf( stemp,"%s_%d", oname, whichSeq);
            freopen( stemp, "w", fout);
            fprintf( stderr,"Writing sequence %d to file %s\n", whichSeq, stemp);
            }
          }

        if (closein) {
          /* !! this fails... skips most seqs... */
          /* !! in sequential read, must count seqs already read from whichSeq ... */
          /* need major revision of ureadseq before we can do this */
          atseq= whichSeq-1;
          seqidptr= seqid;
          seq = readSeqFp( whichSeq, fin, skiplines, format,
                          &seqlen, &atseq, &err, seqidptr);
          skiplines= 0;
          }
        else {
          atseq= 0;
          seqidptr= seqid;
#ifdef NCBI
          if (format == kASNseqentry || format == kASNseqset) {
            seqidptr= NULL;
            seq = readASNSeq( whichSeq, inputfile, skiplines, format,
                     &seqlen, &atseq, &err, &seqidptr);
            }
          else
#endif
          seq = readSeq( whichSeq, inputfile, skiplines, format,
                          &seqlen, &atseq, &err, seqidptr);
          }


        if (gPretty.degap) {
          char *newseq;
          long newlen;
          newseq= compressSeq( gPretty.gapchar, seq, seqlen, &newlen);
          if (newseq) {
            free(seq); seq= newseq; seqlen= newlen;
            }
          }

        if (outform == kMSF) checksum= GCGchecksum(seq, seqlen, &checkall);
        else if (verbose) checksum= seqchecksum(seq, seqlen, &checkall);
        if (verbose)
          fprintf( stderr, "Sequence %d, length= %d, checksum= %X, format= %s, id= %s\n",
                whichSeq, seqlen, checksum, formatstr(format), seqidptr);

        if (err != 0) erralert(err);
        else {
                                  /* format fixes that writeseq doesn't do */
          switch (outform) {
            case kPIR:
              if (seqout == 0) fprintf( foo,"\\\\\\\n");
              break;
            case kASN1:
              if (seqout == 0) fprintf( foo, kASN1headline);
              break;

            case kPhylip:
              if (seqout == 0) {
                if (!interleaved) {  /*  bug, nseq is for 1st infile only */
                  if (chooseall) i= nseq; else i=1;
                  if (phylvers >= 4) fprintf(foo," %d %d\n", i, seqlen);
                  else fprintf(foo," %d %d YF\n", i, seqlen);
                  }
                seqlen0 = seqlen;
                }
              else if (seqlen != seqlen0) {
                erralert(eUnequalSize);
                if (seqlen < seqlen0) seq = (char *)realloc(seq, seqlen0);
                for (i=seqlen; i<seqlen0; i++) seq[i]= gPretty.gapchar;
                seqlen = seqlen0;
                seq[seqlen] = 0; 
                }
              break;

            case kPAUP:
              if (seqout == 0) {
                seqtype= getseqtype(seq, seqlen);
                seqlen0 = seqlen;
                }
              else if (seqlen != seqlen0) {
                erralert(eUnequalSize);
                if (seqlen < seqlen0) seq = (char *)realloc(seq, seqlen0); 
                for (i=seqlen; i<seqlen0; i++) seq[i]= gPretty.gapchar;
                seqlen = seqlen0;
                seq[seqlen] = 0; 
                }
              break;

            }

          if (doupper)
            for (i = 0; i<seqlen; i++) seq[i] = to_upper(seq[i]);
          else if (dolower)
            for (i = 0; i<seqlen; i++) seq[i] = to_lower(seq[i]);

          if (doreverse) {
            long  j, k;
            char  ctemp;
            for (j=0, k=seqlen-1; j <= k; j++, k--) {
              ctemp = compl[seq[j] - ' '];
              seq[j] = compl[seq[k] - ' '];
              seq[k] = ctemp;
              }
            }

          if ((gPretty.isactive || outform==kPAUP) && gPretty.domatch && firstseq != NULL) {
            for (i=0; i<seqlen; i++)
              if (seq[i]==firstseq[i]) seq[i]= gPretty.matchchar;
            }


          if (gPretty.isactive && gPretty.numtop && seqout == 0) {
            gPretty.numline = 1;
            indexout();
            (void) writeSeq( fout, seq, seqlen, outform, seqidptr);
            gPretty.numline = 2;
            indexout();
            (void) writeSeq( fout, seq, seqlen, outform, seqidptr);
            gPretty.numline = 0;
            }

          indexout();
          nlines = writeSeq( fout, seq, seqlen, outform, seqidptr);
          seqout++;
          }

        if ((gPretty.isactive || outform==kPAUP) && gPretty.domatch && firstseq == NULL) {
          firstseq= seq;
          seq = NULL;
          }
        else if (seq!=NULL) { free(seq); seq = NULL; }

#ifdef NCBI
       if ( (format == kASNseqentry || format == kASNseqset)
          && seqidptr && seqidptr!= seqid)
            free(seqidptr);
#endif
        if (chooseall) whichSeq++;
        else if (iwhichlist<nwhichlist) whichSeq= whichlist[iwhichlist++];
        else whichSeq= 0;
        }
      if (closein) { fclose(fin); closein= false; }
      }
    whichSeq  = 0;
  } while (nfile > 0 || !quietly);


fini:
  if (firstseq) { free(firstseq); firstseq= NULL; }
  if (err || listonly) exit_main(err);

  if (gPretty.isactive && gPretty.numbot) {
    gPretty.numline = 2;
    indexout();
    (void) writeSeq( fout, seq, seqlen, outform, seqidptr);
    gPretty.numline = 1;
    indexout();
    (void) writeSeq( fout, seq, seqlen, outform, seqidptr);
    gPretty.numline = 0;
    }

  if (outform == kMSF) {
    if (*oname) cp= oname; else cp= inputfile;
    fprintf(foo,"\n %s  MSF: %d  Type: N  January 01, 1776  12:00  Check: %d ..\n\n",
                  cp, seqlen, checkall);
    }

  if (outform == kPAUP) {
    fprintf(foo,"#NEXUS\n");
    if (*oname) cp= oname; else cp= inputfile;
    fprintf(foo,"[%s -- data title]\n\n", cp);
    /* ! now have header lines for each sequence... put them before "begin data;... */
    }

  if (outform==kPhylip && interleaved) {
    if (phylvers >= 4) fprintf(foo," %d %d\n", seqout, seqlen);
    else fprintf(foo," %d %d YF\n", seqout, seqlen);
    }

  if (interleaved) {
    /* interleave species lines in true output */
    /* nlines is # lines / sequence */
    short iline, j, leaf, iseq;
    char  *s = stempstore;

    indexout();  noutindex--; /* mark eof */

    for (leaf=0; leaf<nlines; leaf++) {
      if (outform == kMSF && leaf == 1) {
        fputs("//\n\n", foo);
        }
      if (outform == kPAUP && leaf==1) {
        switch (seqtype) {
          case kDNA     : cp= "dna"; break;
          case kRNA     : cp= "rna"; break;
          case kNucleic : cp= "dna"; break;
          case kAmino   : cp= "protein"; break;
          case kOtherSeq: cp= "dna"; break;
          }
        fprintf(foo,"\nbegin data;\n");
        fprintf(foo," dimensions ntax=%d nchar=%d;\n", seqout, seqlen);
        fprintf(foo," format datatype=%s interleave missing=%c", cp, gPretty.gapchar);
        if (gPretty.domatch) fprintf(foo," matchchar=%c", gPretty.matchchar);
        fprintf(foo,";\n  matrix\n");
        }

      for (iseq=0; iseq<noutindex; iseq++) {
        fseek(ftmp, outindex[iseq], 0);
        for (iline=0; iline<=leaf; iline++)
          if (!fgets(s, 256, ftmp)) *s= 0;
        if (ftell(ftmp) <= outindex[iseq+1])
          fputs( s, foo);
        }

      for (j=0; j<gPretty.interline; j++)
        fputs( "\n", foo);  /* some want spacer line */
      }
    fclose(ftmp); /* tmp disappears */
    fout= foo;
    }

  if (outform == kASN1)  fprintf( foo, "} }\n");
  if (outform == kPAUP)  fprintf( foo,";\n  end;\n");

  if (outindex != NULL) free(outindex);
  exit_main(0);
}


