/*
			       aln_util.c
	       Helper functions frelated to alignments
*/
/* Last changed Time-stamp: <2006-01-16 11:42:38 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <config.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/aln_util.h"

#define MAX_NUM_NAMES    500
int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]) {
   char *line, name[100]="", *seq;
   int  n, nn=0, num_seq = 0, i;

   if ((line=get_line(clust)) == NULL) {
     fprintf(stderr, "Empty CLUSTAL file\n"); return 0;
   }

   if ((strncmp(line,"CLUSTAL", 7) !=0) && (!strstr(line,"STOCKHOLM"))) {
     fprintf(stderr, "This doesn't look like a CLUSTAL/STOCKHOLM file, sorry\n");
     free(line); return 0;
   }
   free(line);
   line = get_line(clust);

   while (line!=NULL) {
    if(strncmp(line, "//", 2) == 0){
      free(line);
      break;
    }

    if (((n=strlen(line))<4) || isspace((int)line[0])) {
      /* skip non-sequence line */
      free(line); line = get_line(clust);
      nn=0; /* reset seqence number */
      continue;
    }
    /* skip comments */
    if(line[0] == '#'){
      free(line);
      line = get_line(clust);
      continue;
    }

     seq = (char *) vrna_alloc( (n+1)*sizeof(char) );
     sscanf(line,"%99s %s", name, seq);

    for(i=0;i<strlen(seq);i++){
      if(seq[i] == '.') seq[i] = '-'; /* replace '.' gaps by '-' */
      /* comment the next line and think about something more difficult to deal with
         lowercase sequence letters if you really want to */
      seq[i] = toupper(seq[i]);
    }

     if (nn == num_seq) { /* first time */
       names[nn] = strdup(name);
       AlignedSeqs[nn] = strdup(seq);
     }
     else {
       if (strcmp(name, names[nn])!=0) {
	 /* name doesn't match */
	 fprintf(stderr,
		 "Sorry, your file is fucked up (inconsitent seq-names)\n");
	 free(line); free(seq);
	 return 0;
       }
       AlignedSeqs[nn] = (char *)
	 vrna_realloc(AlignedSeqs[nn], strlen(seq)+strlen(AlignedSeqs[nn])+1);
       strcat(AlignedSeqs[nn], seq);
     }
     nn++;
     if (nn>num_seq) num_seq = nn;
     free(seq);
     free(line);
     if (num_seq>=MAX_NUM_NAMES) {
       fprintf(stderr, "Too many sequences in CLUSTAL/STOCKHOLM file");
       return 0;
     }

     line = get_line(clust);
   }

   AlignedSeqs[num_seq] = NULL;
   names[num_seq] = NULL;
   if (num_seq == 0) {
     fprintf(stderr, "No sequences found in CLUSTAL/STOCKHOLM file\n");
     return 0;
   }
   n = strlen(AlignedSeqs[0]);
   for (nn=1; nn<num_seq; nn++) {
     if (strlen(AlignedSeqs[nn])!=n) {
       fprintf(stderr, "Sorry, your file is fucked up.\n"
	       "Unequal lengths!\n\n");
       return 0;
     }
   }

   fprintf(stderr, "%d sequences; length of alignment %d.\n", nn, n);
   return num_seq;
}

char *consensus(const char *AS[]) {
  /* simple consensus sequence (most frequent character) */
  char *string;
  int i,n;
  n = strlen(AS[0]);
  string = (char *) vrna_alloc((n+1)*sizeof(char));
  for (i=0; i<n; i++) {
    int s,c,fm, freq[8] = {0,0,0,0,0,0,0,0};
    for (s=0; AS[s]!=NULL; s++)
      freq[encode_char(AS[s][i])]++;
    for (s=c=fm=0; s<8; s++) /* find the most frequent char */
      if (freq[s]>fm) {c=s, fm=freq[c];}
    if (s>4) s++; /* skip T */
    string[i]=Law_and_Order[c];
  }
  return string;
}

/* IUP nucleotide classes indexed by a bit string of the present bases */
/* A C AC G AG CG ACG U AU CU ACU GU AGU CGU ACGU */
static char IUP[17] = "-ACMGRSVUWYHKDBN";
char *consens_mis(const char*AS[]) {
  /* MIS displays the 'most informative sequence' (Freyhult et al 2004),
     elements in columns with frequency greater than the background
     frequency are projected into iupac notation. Columns where gaps are
     over-represented are in lower case. */

  char *cons;
  int i, s, n, N, c;
  int bgfreq[8] = {0,0,0,0,0,0,0,0};

  n = strlen(AS[0]);
  for (N=0; AS[N]!=NULL; N++);
  cons = (char *) vrna_alloc((n+1)*sizeof(char));

  for (i=0; i<n; i++)
    for (s=0; s<N; s++) {
      c = encode_char(AS[s][i]);
      if (c>4) c=5;
      bgfreq[c]++;
    }

  for (i=0; i<n; i++) {
    int freq[8] = {0,0,0,0,0,0,0,0};
    int code = 0;
    for (s=0; s<N; s++) {
      c = encode_char(AS[s][i]);
      if (c>4) c=5;
      freq[c]++;
    }
    for (c=4; c>0; c--) {
      code <<=1;
      if (freq[c]*n>=bgfreq[c]) code++;
    }
    cons[i] = IUP[code];
    if (freq[0]*n>bgfreq[0])
      cons[i] = tolower(IUP[code]);
  }
  return cons;
}

PUBLIC char *
get_ungapped_sequence(const char *seq){

  char  *tmp_sequence, *b;
  int   i;

  tmp_sequence = strdup(seq);

  b = tmp_sequence;
  i = 0;
  do{
    if((*b=='-')||(*b=='_')||(*b=='~')||(*b=='.')) continue;
    tmp_sequence[i] = *b;
    i++;
  }while(*(++b));

  tmp_sequence = (char *)vrna_realloc(tmp_sequence, (i+1)*sizeof(char));
  tmp_sequence[i] = '\0';

  return tmp_sequence;
}

PUBLIC int
vrna_aln_mpi( char *Alseq[],
              int n_seq,
              int length,
              int *mini){

  int   i, j, k, pairnum = 0, sumident = 0;
  float ident = 0, minimum = 1.;

  for(j=0; j<n_seq-1; j++)
    for(k=j+1; k<n_seq; k++) {
      ident=0;
      for (i=1; i<=length; i++){
        if (Alseq[k][i]==Alseq[j][i]) ident++;
        pairnum++;
      }
      if ((ident/length)<minimum) minimum=ident/(float)length;
      sumident+=ident;
    }
  mini[0]=(int)(minimum*100.);
  if (pairnum>0)   return (int) (sumident*100/pairnum);
  else return 0;

}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC int
get_mpi(char *Alseq[],
        int n_seq,
        int length,
        int *mini){

  return vrna_aln_mpi(Alseq, n_seq, length, mini);
}

PUBLIC void
alloc_sequence_arrays(const char **sequences,
                      short ***S,
                      short ***S5,
                      short ***S3,
                      unsigned short ***a2s,
                      char ***Ss,
                      int circ){

  unsigned int s, n_seq, length;
  if(sequences[0] != NULL){
    length = strlen(sequences[0]);
    for (s=0; sequences[s] != NULL; s++);
    n_seq = s;
    *S    = (short **)          vrna_alloc((n_seq+1) * sizeof(short *));
    *S5   = (short **)          vrna_alloc((n_seq+1) * sizeof(short *));
    *S3   = (short **)          vrna_alloc((n_seq+1) * sizeof(short *));
    *a2s  = (unsigned short **) vrna_alloc((n_seq+1) * sizeof(unsigned short *));
    *Ss   = (char **)           vrna_alloc((n_seq+1) * sizeof(char *));
    for (s=0; s<n_seq; s++) {
      if(strlen(sequences[s]) != length) vrna_message_error("uneqal seqence lengths");
      (*S5)[s]  = (short *)         vrna_alloc((length + 2) * sizeof(short));
      (*S3)[s]  = (short *)         vrna_alloc((length + 2) * sizeof(short));
      (*a2s)[s] = (unsigned short *)vrna_alloc((length + 2) * sizeof(unsigned short));
      (*Ss)[s]  = (char *)          vrna_alloc((length + 2) * sizeof(char));
      (*S)[s]   = (short *)         vrna_alloc((length + 2) * sizeof(short));
      encode_ali_sequence(sequences[s], (*S)[s], (*S5)[s], (*S3)[s], (*Ss)[s], (*a2s)[s], circ);
    }
    (*S5)[n_seq]  = NULL;
    (*S3)[n_seq]  = NULL;
    (*a2s)[n_seq] = NULL;
    (*Ss)[n_seq]  = NULL;
    (*S)[n_seq]   = NULL;
  }
  else vrna_message_error("alloc_sequence_arrays: no sequences in the alignment!");
}

PUBLIC void
free_sequence_arrays( unsigned int n_seq,
                      short ***S,
                      short ***S5,
                      short ***S3,
                      unsigned short ***a2s,
                      char ***Ss){

  unsigned int s;
  for (s=0; s<n_seq; s++) {
    free((*S)[s]);
    free((*S5)[s]);
    free((*S3)[s]);
    free((*a2s)[s]);
    free((*Ss)[s]);
  }
  free(*S);   *S    = NULL;
  free(*S5);  *S5   = NULL;
  free(*S3);  *S3   = NULL;
  free(*a2s); *a2s  = NULL;
  free(*Ss);  *Ss   = NULL;
}

PUBLIC void
encode_ali_sequence(const char *sequence,
                    short *S,
                    short *s5,
                    short *s3,
                    char *ss,
                    unsigned short *as,
                    int circular){

  unsigned int i,l;
  unsigned short p;
  l     = strlen(sequence);
  S[0]  = (short) l;
  s5[0] = s5[1] = 0;

  /* make numerical encoding of sequence */
  for(i=1; i<=l; i++){
    short ctemp;
    ctemp=(short) encode_char(toupper(sequence[i-1]));
    S[i]= ctemp ;
  }

  if (oldAliEn){
    /* use alignment sequences in all energy evaluations */
    ss[0]=sequence[0];
    for(i=1; i<l; i++){
      s5[i] = S[i-1];
      s3[i] = S[i+1];
      ss[i] = sequence[i];
      as[i] = i;
    }
    ss[l]   = sequence[l];
    as[l]   = l;
    s5[l]   = S[l-1];
    s3[l]   = 0;
    S[l+1]  = S[1];
    s5[1]   = 0;
    if (circular) {
      s5[1]   = S[l];
      s3[l]   = S[1];
      ss[l+1] = S[1];
    }
  }
  else{
    if(circular){
      for(i=l; i>0; i--){
        char c5;
        c5 = sequence[i-1];
        if ((c5=='-')||(c5=='_')||(c5=='~')||(c5=='.')) continue;
        s5[1] = S[i];
        break;
      }
      for (i=1; i<=l; i++) {
        char c3;
        c3 = sequence[i-1];
        if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.')) continue;
        s3[l] = S[i];
        break;
      }
    }
    else  s5[1]=s3[l]=0;

    for(i=1,p=0; i<=l; i++){
      char c5;
      c5 = sequence[i-1];
      if ((c5=='-')||(c5=='_')||(c5=='~')||(c5=='.'))
        s5[i+1]=s5[i];
      else { /* no gap */
        ss[p++]=sequence[i-1]; /*start at 0!!*/
        s5[i+1]=S[i];
      }
      as[i]=p;
    }
    for (i=l; i>=1; i--) {
      char c3;
      c3 = sequence[i-1];
      if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.'))
        s3[i-1]=s3[i];
      else
        s3[i-1]=S[i];
    }
  }
}

