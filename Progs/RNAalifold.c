/* Last changed Time-stamp: <2001-04-05 09:33:03 ivo> */
/*                
		  Access to alifold Routines

		  c Ivo L Hofacker
		  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
#include "pair_mat.h"
#include "alifold.h"
extern void  read_parameter_file(const char fname[]);
/*@unused@*/
static char rcsid[] = "$Id: RNAalifold.c,v 1.2 2001/04/05 07:33:16 ivo Exp $";

#define PRIVATE static

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);
PRIVATE int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]);
PRIVATE char *consensus(const char *AS[]);
PRIVATE char **annote(const char *structure, const char *AS[]);
PRIVATE void print_pi(pair_info pi, FILE *file);
/*--------------------------------------------------------------------------*/
#define MAX_NUM_NAMES    500
int main(int argc, char *argv[])
{
  char *string;
  char *structure=NULL, *cstruc=NULL;
  char  ffname[20], gfname[20], fname[13]="";
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   n_seq, i, length, sym, r;
  double min_en, sfact=1.07;
  int   pf=0, istty;
  int noconv=0;
  char     *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char     *names[MAX_NUM_NAMES];       /* sequence names */
  FILE     *clust_file = stdin;

  do_backtrack = 1; 
  string=NULL;
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') {
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'p':  pf=1;
	  if (argv[i][2]!='\0')
	    (void) sscanf(argv[i]+2, "%d", &do_backtrack);
	  break;
	case 'n':
	  if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	  if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
	  if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	  if ( strcmp(argv[i], "-nsp") ==0) {
	    if (i==argc-1) usage();
	    ns_bases = argv[++i];
	  }
	  if ( strcmp(argv[i], "-noconv")==0) noconv=1;
	  if ( strcmp(argv[i], "-nc")==0) {
	    r=sscanf(argv[++i], "%lf", &nc_fact);
	    if (!r) usage();
	  }
	  break;
	case '4':
	  tetra_loop=0;
	  break;
	case 'e':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &energy_set);
	  if (!r) usage();
	  break;
	case 'C':
	  fold_constrained=1;
	  break;
	case 'S':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%lf", &sfact);
	  if (!r) usage();
	  break;
	case 'd': dangles=0;
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &dangles);
	    if (r!=1) usage();
	  }
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	case 'c':
	  if ( strcmp(argv[i], "-cv")==0) {
	    r=sscanf(argv[++i], "%lf", &cv_fact);
	    if (!r) usage();
	  }
	  break;
	default: usage();
	}
    }
    else { /* doesn't start with '-' should be filename */ 
      if (i!=argc-1) usage();
      clust_file = fopen(argv[i], "r");
      if (!clust_file) {
	fprintf(stderr, "can't open %s\n", argv[i]);
	usage();
      }
    }
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);
   
  if (ns_bases != NULL) {
    nonstandards = space(33);
    c=ns_bases;
    i=sym=0;
    if (*c=='-') {
      sym=1; c++;
    }
    while (*c!='\0') {
      if (*c!=',') {
	nonstandards[i++]=*c++;
	nonstandards[i++]=*c;
	if ((sym)&&(*c!=*(c-1))) {
	  nonstandards[i++]=*c;
	  nonstandards[i++]=*(c-1);
	}
      }
      c++;
    }
  }
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
  if ((fold_constrained)&&(istty)) {
    printf("Input constraints using the following notation:\n");
    printf("| : paired with another base\n");
    printf(". : no constraint at all\n");
    printf("x : base must not pair\n");
    printf("< : base i is paired with a base j<i\n");
    printf("> : base i is paired with a base j>i\n");
    printf("matching brackets ( ): base i pairs base j\n");
  } 

  if (istty && (clust_file == stdin)) {
    printf("\nInput aligned sequences in clustalw format\n");
    printf("%s%s\n", scale1, scale2);
  }
  
  n_seq = read_clustal(clust_file, AS, names);
  if (n_seq==0)
    nrerror("no sequences found");

  length = (int) strlen(AS[0]);
  structure = (char *) space((unsigned) length+1);
  if (fold_constrained) {
    cstruc = get_line(stdin);
    if (cstruc!=NULL) 
      strncpy(structure, cstruc, length);
    else
      fprintf(stderr, "constraints missing\n");
  }
  
  min_en = alifold(AS, structure);
  string = consensus((const char **) AS);
  printf("%s\n%s", string, structure);
  if (istty)
    printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
  else
    printf(" (%6.2f)\n", min_en);
  
  if (fname[0]!='\0') {
    strcpy(ffname, fname);
    strcat(ffname, "_ss.ps");
    strcpy(gfname, fname);
    strcat(gfname, "_ss.g");
  } else {
    strcpy(ffname, "alirna.ps");
    strcpy(gfname, "alirna.g");
  }
  if (length<2000) {
    char **A, **a;
    a = A = annote(structure, (const char**) AS);
    (void) PS_rna_plot_a(string, structure, ffname, A);
    while (*A) {free(*A); A++;}
    free(a);
  } else 
    fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");
  
  { /* free mfe arrays but preserve base_pair for PS_dot_plot */
    struct bond  *bp;
    bp = base_pair; base_pair = space(16);
    free_arrays();  /* free's base_pair */
    base_pair = bp;
  }
  if (pf) {
    double energy, kT;
    pair_info *pi;
    char * mfe_struc;

    mfe_struc = strdup(structure);
    	 
    kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
    pf_scale = exp(-(sfact*min_en)/kT/length);
    if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
    fflush(stdout);
    
    /* init_alipf_fold(length); */
    
    if (cstruc!=NULL)
      strncpy(structure, cstruc, length+1);
    energy = alipf_fold(AS, structure, &pi);
    
    if (do_backtrack) {
      printf("%s", structure);
      if (!istty) printf(" [%6.2f]\n", energy);
      else printf("\n");
    }
    if ((istty)||(!do_backtrack)) 
      printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
    printf(" frequency of mfe structure in ensemble %g\n",
	   exp((energy-min_en)/kT));
    
    if (do_backtrack) {
      FILE *aliout;
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_dp.ps");
      } else strcpy(ffname, "alidot.ps");
      (void) PS_dot_plot(string, ffname);
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_ali.out");
      } else strcpy(ffname, "alifold.out");
      aliout = fopen(ffname, "w");
      if (!aliout) {
	fprintf(stderr, "can't open %s    skipping output\n", ffname);
      } else {
	short *ptable; int k;
	ptable = make_pair_table(mfe_struc);
	fprintf(aliout, "%d sequence; length of alignment %d\n", n_seq, length);
	fprintf(aliout, "alifold output\n");
	for (k=0; pi[k].i>0; k++) {
	  pi[k].comp = (ptable[pi[k].i] == pi[k].j) ? 1:0;
	  print_pi(pi[k], aliout);
	}
	fprintf(aliout, "%s\n", structure);
	free(ptable);
      }
    }
    free(mfe_struc);
    free(pi);
  }
  if (cstruc!=NULL) free(cstruc);
  free(base_pair);
  (void) fflush(stdout);
  free(string);
  free(structure); 
  return 0;
}

PRIVATE int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]) {
   char *line, workstr[100], *strg;
   int i, n, nn;
   
   line=get_line(clust);
   sscanf(line,"%99s",workstr);
   
   if(strcmp(workstr,"CLUSTAL") !=0) {
      fprintf(stderr,"Sorry. Expecting CLUSTAL file.\n");
      return 0;
   }
   free(line); 
   line = get_line(clust); 
   while ( (line[0]=='\0') || isspace((int)line[0])){
      free(line);
      line = get_line(clust); 
   } 

   nn=0;
   while (!isspace((int)line[0]) ) {
     n = strlen(line);
     strg = (char *) space( (n+1)*sizeof(char) );
     sscanf(line,"%99s %s",workstr,strg);
     names[nn] = (char *) space(strlen(workstr)+1);
     strcpy(names[nn],workstr); 
     AlignedSeqs[nn] =  (char *)space(strlen(strg)+1);
     strcat(AlignedSeqs[nn],strg);
     nn++;
     free(strg);
     free(line);
     line = get_line(clust);
   }

   /* ignore the consensus-line for the moment */
   free(line);
   line = get_line(clust);  

   while(line!=NULL) { 

      while ( (line[0]=='\0') || isspace((int)line[0])){
         free(line);
         line = get_line(clust); 
      } 

      nn=0;
      while (!isspace((int)line[0])) {
         n = strlen(line);
         strg = (char *) space( (n+1)*sizeof(char) );
         sscanf(line,"%s %s",workstr,strg);
         if (strcmp(workstr,names[nn])==0) {
            AlignedSeqs[nn] = (char *)
               realloc(AlignedSeqs[nn], 
                       strlen(strg)+strlen(AlignedSeqs[nn])+1);
            strcat(AlignedSeqs[nn],strg);   
         } else {
            fprintf(stderr, "Sorry, your file is fucked up.\n");
            return 0;
         }
         nn++;
         free(strg);
         free(line);
         line = get_line(clust);
     }
     /* ignore the consensus-line for the moment */
     free(line);
     line = get_line(clust);      
   }   
   
   AlignedSeqs[nn] = NULL;
   n = strlen(AlignedSeqs[0]); 
   for (i=1; i<nn; i++) {
      if (strlen(AlignedSeqs[i])!=n) {
         fprintf(stderr, "Sorry, your file is fucked up.\n"
                 "Unequal lengths!\n\n");
         return 0;
      }
   }

   fprintf(stderr, "%d sequences; length of alignment %d.\n", nn, n);
   return nn;
}
 
void print_pi(pair_info pi, FILE *file) {
  char *pname[8] = {"","CG","GC","GU","UG","AU","UA", "--"};
  int i;
  
  /* numbering starts with 1 in output */
  fprintf(file, "%5d %5d %2d %5.1f%% %7.3f",
	  pi.i, pi.j, pi.bp[0], 100*pi.p, pi.ent);
  for (i=1; i<=7; i++) 
    if (pi.bp[i]) fprintf(file, " %s:%-4d", pname[i], pi.bp[i]);
  /* if ((!pi.sym)&&(pi.j>=0)) printf(" *"); */
  if (!pi.comp) fprintf(file, " +");
  fprintf(file, "\n");
}

PRIVATE char *consensus(const char *AS[]) {
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

PRIVATE char **annote(const char *structure, const char *AS[]) {
  char **A;
  int i,n,s,n_A, i_A=0;
  short *ptable;
  make_pair_matrix();
  n = strlen(AS[0]);
  n_A = 128;
  A = (char **) space(n_A*sizeof(char *));
  ptable = make_pair_table(structure);
  A[i_A++] = strdup("\n");
  for (i=1; i<=n; i++) {
    char ps[256], ci='\0', cj='\0';
    int j, type, pfreq[8] = {0,0,0,0,0,0,0,0}, vi=0, vj=0;
    if ((j=ptable[i])<i) continue;
    for (s=0; AS[s]!=NULL; s++) {
      type = pair[encode_char(AS[s][i-1])][encode_char(AS[s][j-1])];
      pfreq[type]++;
      if (type) {
	if (AS[s][i-1] != ci) { ci = AS[s][i-1]; vi++;}
	if (AS[s][j-1] != cj) { cj = AS[s][j-1]; vj++;}
      }
    }
    if (i_A>=n_A-3) {
      n_A *= 2; A = realloc(A, n_A *sizeof(char *));
      if (A==NULL) nrerror("out of memory in annote");
    }
    if (pfreq[0]>0) {
      sprintf(ps, "%d %d %d gmark", i, j, pfreq[0]);
      A[i_A++] = strdup(ps);
    }
    if (vi>1) {
      sprintf(ps, "%d cmark", i);
      A[i_A++] = strdup(ps);
    }
    if (vj>1) {
      sprintf(ps, "%d cmark", j);
      A[i_A++] = strdup(ps);
    }
  }
  
  A[i_A] = NULL;
  free(ptable);
  return A;
}
/*-------------------------------------------------------------------------*/

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAalifold [-cv float] [-nc float]\n"
	  "        [-p[0]] [-C] [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n" 
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] "
	  "[-noconv]\n");
}
