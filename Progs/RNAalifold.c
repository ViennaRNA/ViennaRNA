/* Last changed Time-stamp: <2005-03-21 13:49:55 ivo> */
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
#include "aln_util.h"
extern void  read_parameter_file(const char fname[]);
/*@unused@*/
static const char rcsid[] = "$Id: RNAalifold.c,v 1.13 2005/03/21 12:50:44 ivo Exp $";

#define PRIVATE static

static const char scale[] = "....,....1....,....2....,....3....,....4"
                            "....,....5....,....6....,....7....,....8";

PRIVATE void /*@exits@*/ usage(void);
PRIVATE char *annote(const char *structure, const char *AS[]);
PRIVATE void print_pi(const pair_info pi, FILE *file);
PRIVATE cpair *make_color_pinfo(const pair_info *pi);
PRIVATE void mark_endgaps(char *seq, char egap);
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
  int   endgaps=0, mis=0;
  double min_en, real_en, sfact=1.07;
  int   pf=0, istty;
  char     *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char     *names[MAX_NUM_NAMES];       /* sequence names */
  FILE     *clust_file = stdin;

  do_backtrack = 1; 
  string=NULL;
  dangles=2;
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
	  if ( strcmp(argv[i], "-nc")==0) {
	    r=sscanf(argv[++i], "%lf", &nc_fact);
	    if (!r) usage();
	  }
	  break;
	case 'm':
	  if ( strcmp(argv[i], "-mis")==0) mis=1;
	  else usage();
	case '4':
	  tetra_loop=0;
	  break;
	case 'e':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &energy_set);
	  if (!r) usage();
	  break;
	case 'E':
	  endgaps=1;
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
      if (clust_file == NULL) {
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
  if (fold_constrained) {
    if (istty) printf("%s\n", scale);
    cstruc = get_line(stdin);
  }

  if (istty && (clust_file == stdin)) {
    printf("\nInput aligned sequences in clustalw format\n");
    if (!fold_constrained) printf("%s\n", scale);
  }
  
  n_seq = read_clustal(clust_file, AS, names);
  if (endgaps)
    for (i=0; i<n_seq; i++) mark_endgaps(AS[i], '~');
  if (clust_file != stdin) fclose(clust_file);
  if (n_seq==0)
    nrerror("no sequences found");

  length = (int) strlen(AS[0]);
  structure = (char *) space((unsigned) length+1);
  if (fold_constrained) {
    if (cstruc!=NULL) 
      strncpy(structure, cstruc, length);
    else
      fprintf(stderr, "constraints missing\n");
  }
  
  min_en = alifold(AS, structure);
  {
    int i; double s=0;
    extern int eos_debug;
    eos_debug=-1; /* shut off warnings about nonstandard pairs */
    for (i=0; AS[i]!=NULL; i++) 
      s += energy_of_struct(AS[i], structure);
    real_en = s/i;
  }
  string = (mis) ?
    consens_mis((const char **) AS) : consensus((const char **) AS);
  printf("%s\n%s", string, structure);
  if (istty)
    printf("\n minimum free energy = %6.2f kcal/mol (%6.2f + %6.2f)\n", 
	   min_en, real_en, min_en - real_en);
  else
    printf(" (%6.2f = %6.2f + %6.2f) \n", min_en, real_en, min_en-real_en );
  
  if (fname[0]!='\0') {
    strcpy(ffname, fname);
    strcat(ffname, "_ss.ps");
    strcpy(gfname, fname);
    strcat(gfname, "_ss.g");
  } else {
    strcpy(ffname, "alirna.ps");
    strcpy(gfname, "alirna.g");
  }
  if (length<=2500) {
    char *A;
    A = annote(structure, (const char**) AS);
    (void) PS_rna_plot_a(string, structure, ffname, NULL, A);
    free(A);
  } else 
    fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");
  
  { /* free mfe arrays but preserve base_pair for PS_dot_plot */
    struct bond  *bp;
    bp = base_pair; base_pair = space(16);
    free_alifold_arrays();  /* free's base_pair */
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
      cpair *cp;
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
	fprintf(aliout, "%d sequence; length of alignment %d\n", 
		n_seq, length);
	fprintf(aliout, "alifold output\n");
	for (k=0; pi[k].i>0; k++) {
	  pi[k].comp = (ptable[pi[k].i] == pi[k].j) ? 1:0;
	  print_pi(pi[k], aliout);
	}
	fprintf(aliout, "%s\n", structure);
	free(ptable);
      }
      fclose(aliout);
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_dp.ps");
      } else strcpy(ffname, "alidot.ps");
      cp = make_color_pinfo(pi);
      (void) PS_color_dot_plot(string, cp, ffname);
      free(cp);
    }
    free(mfe_struc);
    free(pi);
  }
  if (cstruc!=NULL) free(cstruc);
  free(base_pair);
  (void) fflush(stdout);
  free(string);
  free(structure);
  for (i=0; AS[i]; i++) {
    free(AS[i]); free(names[i]);
  }
  return 0;
}

void mark_endgaps(char *seq, char egap) {
  int i,n;
  n = strlen(seq);
  for (i=0; i<n && (seq[i]=='-'); i++) {
    seq[i] = egap;
  }
  for (i=n-1; i>0 && (seq[i]=='-'); i--) {
    seq[i] = egap;
  }
}
 
void print_pi(const pair_info pi, FILE *file) {
  const char *pname[8] = {"","CG","GC","GU","UG","AU","UA", "--"};
  int i;
  
  /* numbering starts with 1 in output */
  fprintf(file, "%5d %5d %2d %5.1f%% %7.3f",
	  pi.i, pi.j, pi.bp[0], 100.*pi.p, pi.ent);
  for (i=1; i<=7; i++) 
    if (pi.bp[i]) fprintf(file, " %s:%-4d", pname[i], pi.bp[i]);
  /* if ((!pi.sym)&&(pi.j>=0)) printf(" *"); */
  if (!pi.comp) fprintf(file, " +");
  fprintf(file, "\n");
}

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

PRIVATE cpair *make_color_pinfo(const pair_info *pi) {
  cpair *cp;
  int i, n;
  for (n=0; pi[n].i>0; n++);
  cp = (cpair *) space(sizeof(cpair)*(n+1));
  for (i=0; i<n; i++) {
    int j, ncomp;
    cp[i].i = pi[i].i;
    cp[i].j = pi[i].j;
    cp[i].p = pi[i].p;
    for (ncomp=0, j=1; j<=6; j++) if (pi[i].bp[j]) ncomp++;
    cp[i].hue = (ncomp-1.0)/6.2;   /* hue<6/6.9 (hue=1 ==  hue=0) */
    cp[i].sat = 1 - MIN2( 1.0, pi[i].bp[0]/2.5);
    cp[i].mfe = pi[i].comp;
  }
  return cp;
}

PRIVATE char *annote(const char *structure, const char *AS[]) {
  char *ps;
  int i, n, s, maxl;
  short *ptable;
  make_pair_matrix();
  n = strlen(AS[0]);
  maxl = 1024;
  ps = (char *) space(maxl);
  ptable = make_pair_table(structure);
  for (i=1; i<=n; i++) {
    char pps[64], ci='\0', cj='\0';
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
    if (maxl - strlen(ps) < 128) { 
      maxl *= 2;
      ps = realloc(ps, maxl);
      if (ps==NULL) nrerror("out of memory in realloc");
    }
    if (pfreq[0]>0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      strcat(ps, pps);
    }
    if (vi>1) {
      snprintf(pps, 64, "%d cmark\n", i);
      strcat(ps, pps);
    }
    if (vj>1) {
      snprintf(pps, 64, "%d cmark\n", j);
      strcat(ps, pps);
    }
  }
  free(ptable);
  return ps;
}

/*-------------------------------------------------------------------------*/

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAalifold [-cv float] [-nc float] [-E] [-mis]\n"
	  "        [-p[0]] [-C] [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n" 
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]\n"
	  );
}
