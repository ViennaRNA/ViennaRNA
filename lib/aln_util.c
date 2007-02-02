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
#include "utils.h"
#include "fold_vars.h"
#include "pair_mat.h"
/*@unused@*/
static char rcsid[] = "$Id: aln_util.c,v 1.4 2007/02/02 15:18:13 ivo Exp $";

#define MAX_NUM_NAMES    500
int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]) {
   char *line, name[100]="", *seq;
   int  n, nn=0, num_seq = 0;

   if ((line=get_line(clust)) == NULL) {
     fprintf(stderr, "Empty CLUSTAL file\n"); return 0;
   }

   if (strncmp(line,"CLUSTAL", 7) !=0) {
     fprintf(stderr, "This doesn't look like a CLUSTAL file, sorry\n");
     free(line); return 0;
   }
   free(line);
   line = get_line(clust);

   while (line!=NULL) {
     if (((n=strlen(line))<4) || isspace((int)line[0])) {
       /* skip non-sequence line */
       free(line); line = get_line(clust);
       nn=0; /* reset seqence number */
       continue;
     }

     seq = (char *) space( (n+1)*sizeof(char) );
     sscanf(line,"%99s %s", name, seq);
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
	 xrealloc(AlignedSeqs[nn], strlen(seq)+strlen(AlignedSeqs[nn])+1);
       strcat(AlignedSeqs[nn], seq);
     }
     nn++;
     if (nn>num_seq) num_seq = nn;
     free(seq);
     free(line);
     if (num_seq>=MAX_NUM_NAMES) {
       fprintf(stderr, "Too many sequences in CLUSTAL file");
       return 0;
     }

     line = get_line(clust);
   }

   AlignedSeqs[num_seq] = NULL;
   names[num_seq] = NULL;
   if (num_seq == 0) {
     fprintf(stderr, "No sequences found in CLSUATL file\n");
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
  string = (char *) space((n+1)*sizeof(char));
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
  cons = (char *) space((n+1)*sizeof(char));

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
