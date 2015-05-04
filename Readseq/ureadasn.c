/* ureadasn.c
  -- parse, mangle and otherwise rewrite ASN1 file/entries for readseq reading
  -- from NCBI toolkit (ncbi.nlm.nih.gov:/toolkit)
*/

#ifdef NCBI

#include <stdio.h>
#include <ctype.h>
#include <string.h>

/* NCBI toolkit :include: must be on lib path */
#include <ncbi.h>
#include <seqport.h>

#define UREADASN
#include "ureadseq.h"

#pragma segment ureadasn

/* this stuff is hacked up from tofasta.c of ncbitools */
#define   kBaseAny   0
#define   kBaseNucleic 1
#define   kBaseAmino   2

typedef struct tofasta {
    Boolean idonly;
    short   *seqnum;
    short   whichSeq;
    char    **seq, **seqid;
    long    *seqlen;
} FastaDat, PNTR FastaPtr;


void BioseqRawToRaw(BioseqPtr bsp, Boolean idonly,
              short whichSeq, short *seqnum,
              char **seq, char **seqid, long *seqlen)
{
  SeqPortPtr spp;
  SeqIdPtr bestid;
  Uint1 repr, code, residue;
  CharPtr tmp, title;
  long  outlen, outmax;
  char  localid[256], *sp;

  /* !!! this may be called several times for a single sequence
    because SeqEntryExplore looks for parts and joins them...
    assume seq, seqid, seqlen may contain data (or NULL)
  */
  if (bsp == NULL) return;
  repr = Bioseq_repr(bsp);
  if (!(repr == Seq_repr_raw || repr == Seq_repr_const)) return;

  (*seqnum)++;
  if (!(whichSeq == *seqnum || whichSeq == 0)) return;

  bestid = SeqIdFindBest(bsp->id, (Uint1) 0);
  title = BioseqGetTitle(bsp);
  if (idonly) {
    sprintf(localid, " %d)  ", *seqnum);
    tmp= localid + strlen(localid)-1;
    }
  else {
    strcpy(localid," ");
    tmp= localid;
    }
  tmp = SeqIdPrint(bestid, tmp, PRINTID_FASTA_SHORT);
  tmp = StringMove(tmp, " ");
  StringNCpy(tmp, title, 200);
/* fprintf(stderr,"BioseqRawToRaw: localid='%s'\n",localid); */

          /* < seqid is fixed storage */
  /* strcpy( *seqid, localid);  */
          /* < seqid is variable sized */
  outmax= strlen(localid) + 3;
  if (*seqid==NULL) {
    *seqid= (char*) malloc(outmax);
    if (*seqid==NULL) return;
    strcpy(*seqid, localid);
    }
  else {
    outmax += strlen(*seqid) + 2;
    *seqid= (char*) realloc( *seqid, outmax);
    if (*seqid==NULL) return;
    if (!idonly) strcat(*seqid, "; ");
    strcat(*seqid, localid);
    }

  if (idonly) {
    strcat(*seqid,"\n");
    return;
    }

  if (ISA_na(bsp->mol)) code = Seq_code_iupacna;
  else code = Seq_code_iupacaa;
  spp = SeqPortNew(bsp, 0, -1, 0, code);
  SeqPortSeek(spp, 0, SEEK_SET);

  sp= *seq;
  if (sp==NULL) {  /* this is always true now !? */
    outlen= 0;
    outmax= 500;
    sp= (char*) malloc(outmax);
    }
  else {
    outlen= strlen(sp);
    outmax= outlen + 500;
    sp= (char*) realloc( sp, outmax);
    }
  if (sp==NULL) return;

  while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF) {
    if (outlen>=outmax) {
      outmax= outlen + 500;
      sp= (char*) realloc(sp, outmax);
      if (sp==NULL) return;
      }
    sp[outlen++] = residue;
    }
  sp= (char*) realloc(sp, outlen+1);
  if (sp!=NULL) sp[outlen]= '\0';
  *seq= sp;
  *seqlen= outlen;
  SeqPortFree(spp);
  return;
}


static void SeqEntryRawseq(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  FastaPtr tfa;
  BioseqPtr bsp;

  if (!IS_Bioseq(sep)) return;
  bsp = (BioseqPtr)sep->data.ptrvalue;
  tfa = (FastaPtr) data;
  BioseqRawToRaw(bsp, tfa->idonly, tfa->whichSeq, tfa->seqnum,
                  tfa->seq, tfa->seqid, tfa->seqlen);
}

void SeqEntryToRaw(SeqEntryPtr sep, Boolean idonly, short whichSeq, short *seqnum,
                        char **seq, char **seqid, long *seqlen)
{
  FastaDat tfa;

  if (sep == NULL) return;
  tfa.idonly= idonly;
  tfa.seqnum= seqnum;
  tfa.whichSeq= whichSeq;
  tfa.seq   = seq;
  tfa.seqid = seqid;
  tfa.seqlen= seqlen;
  SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryRawseq);
}




char *listASNSeqs(const char *filename, const long skiplines,
                  const short format,   /* note: this is kASNseqentry or kASNseqset */
                  short *nseq, short *error )
{
  AsnIoPtr aip = NULL;
  SeqEntryPtr the_set;
  AsnTypePtr atp, atp2;
  AsnModulePtr amp;
  Boolean inIsBinary= FALSE; /* damn, why can't asn routines test this? */
  char  *seq = NULL;
  char  *seqid = NULL, stemp[256];
  long  seqlen;
  int   i, count;

  *nseq= 0;
  *error= 0;

    /* asn dictionary setups */
/*fprintf(stderr,"listASNSeqs: SeqEntryLoad\n");*/
  if (! SeqEntryLoad()) goto errxit; /*  sequence alphabets (and sequence parse trees) */
  amp = AsnAllModPtr();   /* get pointer to all loaded ASN.1 modules */
  if (amp == NULL) goto errxit;
  atp = AsnFind("Bioseq-set");    /* get the initial type pointers */
  if (atp == NULL) goto errxit;
  atp2 = AsnFind("Bioseq-set.seq-set.E");
  if (atp2 == NULL) goto errxit;

/*fprintf(stderr,"listASNSeqs: AsnIoOpen\n");*/
      /* open the ASN.1 input file in the right mode */
      /* !!!! THIS FAILS when filename has MAC PATH (& other paths?) (:folder:filename) */
  if ((aip = AsnIoOpen(filename, inIsBinary?"rb":"r")) == NULL) goto errxit;
  for (i=0; i<skiplines; i++) fgets( stemp, 255, aip->fp);  /* this may mess up asn routines... */

  if (! ErrSetLog ("stderr"))  goto errxit;
  else ErrSetOpts(ERR_CONTINUE, ERR_LOG_ON);    /*??  log errors instead of die */

  if (format == kASNseqentry) {  /* read one Seq-entry */
/*fprintf(stderr,"listASNSeqs: SeqEntryAsnRead\n");*/
    the_set = SeqEntryAsnRead(aip, NULL);
    SeqEntryToRaw(the_set, true,  0, nseq, &seq, &seqid, &seqlen);
    if (seq) free(seq); seq= NULL;
    SeqEntryFree(the_set);
    }
  else   {                   /* read Seq-entry's from a Bioseq-set */
    count = 0;
/*fprintf(stderr,"listASNSeqs: AsnReadId\n");*/
    while ((atp = AsnReadId(aip, amp, atp)) != NULL) {
      if (atp == atp2)  {  /* top level Seq-entry */
        the_set = SeqEntryAsnRead(aip, atp);
        SeqEntryToRaw(the_set, true, 0, nseq, &seq, &seqid, &seqlen);
        SeqEntryFree(the_set);
        if (seq) free(seq); seq= NULL;
        }
      else
        AsnReadVal(aip, atp, NULL);
      count++;
      }
    }

  AsnIoClose(aip);
  *error= 0;
  return seqid;

errxit:
  AsnIoClose(aip);
  if (seqid) free(seqid);
  *error= eASNerr;
  return NULL;
}


char *readASNSeq(const short whichEntry, const char *filename,
                const long skiplines,
                const short format,     /* note: this is kASNseqentry or kASNseqset */
                long *seqlen, short *nseq,
                short *error, char **seqid )
{
  AsnIoPtr aip = NULL;
  SeqEntryPtr the_set;
  AsnTypePtr atp, atp2;
  AsnModulePtr amp;
  Boolean inIsBinary= FALSE; /* damn, why can't asn routines test this? */
  char  *seq, stemp[200];
  int   i, count;

  *seqlen= 0;
  *nseq= 0;
  *error= 0;
  seq= NULL;

/*fprintf(stderr,"readASNseq: SeqEntryLoad\n");*/
    /* asn dictionary setups */
  if (! SeqEntryLoad()) goto errxit; /*  sequence alphabets (and sequence parse trees) */
  amp = AsnAllModPtr();   /* get pointer to all loaded ASN.1 modules */
  if (amp == NULL) goto errxit;
  atp = AsnFind("Bioseq-set");    /* get the initial type pointers */
  if (atp == NULL) goto errxit;
  atp2 = AsnFind("Bioseq-set.seq-set.E");
  if (atp2 == NULL) goto errxit;

      /* open the ASN.1 input file in the right mode */
/*fprintf(stderr,"readASNseq: AsnIoOpen(%s)\n", filename);*/
  if ((aip = AsnIoOpen(filename, inIsBinary?"rb":"r")) == NULL) goto errxit;
  for (i=0; i<skiplines; i++) fgets( stemp, 255, aip->fp);  /* this may mess up asn routines... */

  if (! ErrSetLog ("stderr"))  goto errxit;
  else ErrSetOpts(ERR_CONTINUE, ERR_LOG_ON);    /*??  log errors instead of die */

  seq= NULL;
  if (format == kASNseqentry) {  /* read one Seq-entry */
/*fprintf(stderr,"readASNseq: SeqEntryAsnRead\n");*/
    the_set = SeqEntryAsnRead(aip, NULL);
    SeqEntryToRaw(the_set, false, whichEntry, nseq, &seq, seqid, seqlen);
    SeqEntryFree(the_set);
    goto goodexit;
    }

  else   {                   /* read Seq-entry's from a Bioseq-set */
    count = 0;
/*fprintf(stderr,"readASNseq: AsnReadId\n");*/
    while ((atp = AsnReadId(aip, amp, atp)) != NULL) {
      if (atp == atp2)  {  /* top level Seq-entry */
        the_set = SeqEntryAsnRead(aip, atp);
        SeqEntryToRaw(the_set, false, whichEntry, nseq, &seq, seqid, seqlen);
        SeqEntryFree(the_set);
        if (*nseq >= whichEntry) goto goodexit;
        }
      else
        AsnReadVal(aip, atp, NULL);
      count++;
      }
    }

goodexit:
  AsnIoClose(aip);
  *error= 0;
  return seq;

errxit:
  AsnIoClose(aip);
  *error= eASNerr;
  if (seq) free(seq);
  return NULL;
}


#endif /*NCBI*/
