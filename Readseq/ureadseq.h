/* File: ureadseq.h
 *
 * Header for module UReadSeq
 */

#ifndef UREADSEQ_H
#define UREADSEQ_H



typedef char  boolean;
#define NEWLINE         '\n'
#define false 0
#define true  1
#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a>b)?a:b
#define skipwhitespace(string)  {while (*string <= ' ' && *string != 0) string++;}

  /* NLM strings */
#define is_upper(c) ('A'<=(c) && (c)<='Z')
#define is_lower(c) ('a'<=(c) && (c)<='z')
#define to_lower(c) ((char)(is_upper(c) ? (c)+' ' : (c)))
#define to_upper(c) ((char)(is_lower(c) ? (c)-' ' : (c)))


  /* readSeq errors */
#define eFileNotFound   -1
#define eNoData         -2
#define eMemFull        -3
#define eItemNotFound   -4
#define eOneFormat      -5
#define eUnequalSize    -6
#define eFileCreate     -7
#define eUnknownFormat  -8
#define eOptionBad      -9
#define eASNerr         -10

  /* magic number for readSeq(whichEntry) to give seq list */
#define kListSequences  -1

  /* sequence types parsed by getseqtype */
#define kOtherSeq   0
#define kDNA        1
#define kRNA        2
#define kNucleic    3
#define kAmino      4

  /* formats known to readSeq */
#define kIG             1
#define kGenBank        2
#define kNBRF           3
#define kEMBL           4
#define kGCG            5
#define kStrider        6
#define kFitch          7
#define kPearson        8
#define kZuker          9
#define kOlsen          10
#define kPhylip2        11
#define kPhylip4        12
#define kPhylip3        kPhylip4
#define kPhylip         kPhylip4
#define kPlain          13  /* keep this at #13 */
#define kPIR            14
#define kMSF            15
#define kASN1           16
#define kPAUP           17
#define kPretty         18

#define kMaxFormat      18
#define kMinFormat      1
#define kNoformat       -1    /* format not tested */
#define kUnknown        0     /* format not determinable */

  /* subsidiary types */
#define kASNseqentry    51
#define kASNseqset      52

#define kPhylipInterleave 61
#define kPhylipSequential 62


typedef struct  {
  boolean isactive, baseonlynum;
  boolean numright, numleft, numtop, numbot;
  boolean nameright, nameleft, nametop;
  boolean noleaves, domatch, degap;
  char  matchchar, gapchar;
  short numline, atseq;
  short namewidth, numwidth;
  short interline, spacer, seqwidth, tab;
  } prettyopts;

#define gPrettyInit(p) { \
  p.isactive=false;\
  p.baseonlynum=true;\
  p.numline= p.atseq= 0;\
  p.numright= p.numleft= p.numtop= p.numbot= false;\
  p.nameright= p.nameleft= p.nametop= false;\
  p.noleaves= p.domatch= p.degap= false;\
  p.matchchar='.';\
  p.gapchar='-';\
  p.namewidth=8;\
  p.numwidth=5;\
  p.interline=1;\
  p.spacer=10;\
  p.seqwidth=50;\
  p.tab=0; }

#ifdef UREADSEQ_G
prettyopts  gPretty;
#else
extern  prettyopts  gPretty;
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern short seqFileFormat(const char *filename, long *skiplines, short *error );
extern short seqFileFormatFp(FILE *fseq, long  *skiplines, short *error );

extern char *listSeqs(const char *filename, const long skiplines,
                       const short format, short *nseq, short *error );

extern char *readSeq(const short whichEntry, const char *filename,
                      const long skiplines, const short format,
                      long *seqlen, short *nseq, short *error, char *seqid );

extern char *readSeqFp(const short whichEntry_, FILE  *fp_,
  const long  skiplines_, const short format_,
        long  *seqlen_,  short *nseq_, short *error_, char *seqid_ );

extern short writeSeq(FILE *outf, const char *seq, const long seqlen,
                       const short outform, const char *seqid );

extern unsigned long CRC32checksum(const char *seq, const long seqlen, unsigned long *checktotal);
extern unsigned long GCGchecksum(const char *seq, const long seqlen, unsigned long *checktotal);
#ifdef SMALLCHECKSUM
#define seqchecksum  GCGchecksum
#else
#define seqchecksum  CRC32checksum
#endif

extern short getseqtype(const char *seq, const long seqlen );
extern char *compressSeq( const char gapc, const char *seq, const long seqlen, long *newlen);

#ifdef NCBI

extern char *listASNSeqs(const char *filename, const long skiplines,
                  const short format, short *nseq, short *error );

extern char *readASNSeq(const short whichEntry, const char *filename,
                const long skiplines, const short format,
                long *seqlen, short *nseq, short *error, char **seqid );
#endif


  /* patches for some missing string.h stuff */
extern int Strcasecmp(const char *a, const char *b);
extern int Strncasecmp(const char *a, const char *b, long maxn);

#ifdef __cplusplus
}
#endif

#endif /*UREADSEQ_H*/

