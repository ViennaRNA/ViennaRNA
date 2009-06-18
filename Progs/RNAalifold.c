/* Last changed Time-stamp: <2009-06-18 14:20:02 ivo> */
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
#include "MEA.h"
extern void  read_parameter_file(const char fname[]);
extern float energy_of_circ_struct(const char *seq, const char *structure);

/*@unused@*/
static const char rcsid[] = "$Id: RNAalifold.c,v 1.23 2009/02/24 14:21:26 ivo Exp $";

#define PRIVATE static

static const char scale[] = "....,....1....,....2....,....3....,....4"
			    "....,....5....,....6....,....7....,....8";

PRIVATE void /*@exits@*/ usage(void);
PRIVATE char **annote(const char *structure, const char *AS[]);
PRIVATE void print_pi(const pair_info pi, FILE *file);
PRIVATE void print_aliout(char **AS, plist *pl, int n_seq, char * mfe, FILE *aliout);
PRIVATE void mark_endgaps(char *seq, char egap);
PRIVATE cpair *make_color_pinfo(char **sequences, plist *pl, int n_seq, bondT *mfe);

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
  int   pf=0, noPS=0, istty;
  char     *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char     *names[MAX_NUM_NAMES];       /* sequence names */
  FILE     *clust_file = stdin;
  int circ=0;
  int doAlnPS=0;
  int doColor=0;
  int n_back=0;
  int eval_energy = 0;
  int doMEA=0;
  double MEAgamma = 1.;


  do_backtrack = 1;
  string=NULL;
  dangles=2;
  oldAliEn=0;
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
	case 'o':
	  if ( strcmp(argv[i], "-old")==0) oldAliEn=1;
	  break;
	case 'm':
	  if ( strcmp(argv[i], "-mis")==0) mis=1;
	  else usage();
	  break;
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
	case 'M':
	  if (strcmp(argv[i], "-MEA")==0) pf=doMEA=1;
	  if (i<argc-1)
	    if (isdigit(argv[i+1][0]) || argv[i+1][0]=='.') {
	      r=sscanf(argv[++i], "%lf", &MEAgamma);
	      if (r!=1) usage();
	    }
	  break;
	case 'c':
	  if ( strcmp(argv[i], "-cv")==0) {
	    r=sscanf(argv[++i], "%lf", &cv_fact);
	    if (!r) usage();
	  } else {
	    if (strcmp(argv[i], "-circ")==0)
	      circ=1;
	    else
	      if (strcmp(argv[i], "-color")==0)
		doColor=1;
	  }
	  break;
	case 'a':
	  if ( strcmp(argv[i], "-aln")==0) {
	    doAlnPS=1;
	  }
	  break;
	case 'R':
	  if (i==argc-1) usage();
	  RibosumFile = argv[++i];
	  ribo=1;
	  break;
	 case 'r':
	  RibosumFile = NULL;
	  ribo=1;
	  break;
	case 's':
	  if (argv[i][2]=='e') eval_energy = 1;
	  else if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r= sscanf(argv[++i], "%d", &n_back);
	  if (!r) usage();
	  do_backtrack=0;
	  pf=1;
	  init_rand();
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

  make_pair_matrix();

  if (circ && noLonelyPairs)
    fprintf(stderr,
	    "warning, depending on the origin of the circular sequence, "
	    "some structures may be missed when using -noLP\n"
	    "Try rotating your sequence a few times\n");
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

  if (circ)
    min_en = circalifold((const char **)AS, structure);
  else
    min_en = alifold(AS, structure);

  if (circ) {
    int i; double s=0;
    extern int eos_debug;
    eos_debug=-1; /* shut off warnings about nonstandard pairs */
    for (i=0; AS[i]!=NULL; i++)
      s +=energy_of_circ_struct(AS[i], structure);
    real_en = s/i;
  } else {
    float CVen;
    real_en = energy_of_alistruct(AS, structure, n_seq, &CVen);
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
  if (!noPS) {
    char **A;
    A = annote(structure, (const char**) AS);
    if (doColor)
      (void) PS_rna_plot_a(string, structure, ffname, A[0], A[1]);
    else
      (void) PS_rna_plot_a(string, structure, ffname, NULL, A[1]);
    free(A[0]); free(A[1]); free(A);
  }
  if (doAlnPS)
    PS_color_aln(structure, "aln.ps", (const char const **) AS, (const char const **) names);

  { /* free mfe arrays but preserve base_pair for PS_dot_plot */
    struct bond  *bp;
    bp = base_pair; base_pair = space(16);
    free_alifold_arrays();  /* free's base_pair */
    base_pair = bp;
  }
  if (pf) {
    float energy, kT;
    plist *pl;
    char * mfe_struc;

    mfe_struc = strdup(structure);

    kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
    pf_scale = exp(-(sfact*min_en)/kT/length);
    if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
    fflush(stdout);

    /* init_alipf_fold(length); */

    if (cstruc!=NULL)
      strncpy(structure, cstruc, length+1);
    energy = (circ) ? alipf_circ_fold(AS, structure, &pl) : alipf_fold(AS, structure, &pl);

    if (n_back>0) {
      /*stochastic sampling*/
      for (i=0; i<n_back; i++) {
	char *s;
	double prob=1.;
	s =alipbacktrack(&prob);
	printf("%s ", s);
	if (eval_energy ) printf("%6g %.2f ",prob, -1*(kT*log(prob)-energy));
	printf("\n");
	free(s);
      }
    }
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
      char *cent;
      double dist;
      float en, CVen;
      cent = centroid_ali(length, &dist, pl);

      en = energy_of_alistruct(AS, cent, n_seq, &CVen);
      printf("%s %6.2f {%6.2f + %6.2f}\n",cent,en-CVen,en,-CVen);
      free(cent);

      if  (doMEA) {
	float mea, mea_en, CVen;
	mea = MEA(pl, structure, MEAgamma);
	mea_en = (circ) ? energy_of_alistruct(AS, structure, n_seq, &CVen) :
	  energy_of_struct(string, structure);
	printf("%s {%6.2f MEA=%.2f}\n", structure, mea_en, mea);
      }

      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_ali.out");
      } else strcpy(ffname, "alifold.out");
      aliout = fopen(ffname, "w");
      if (!aliout) {
	fprintf(stderr, "can't open %s    skipping output\n", ffname);
      } else {
      print_aliout(AS, pl, n_seq, mfe_struc, aliout);
      }
      fclose(aliout);
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_dp.ps");
      } else strcpy(ffname, "alidot.ps");
      cp = make_color_pinfo(AS,pl, n_seq,base_pair);
      (void) PS_color_dot_plot(string, cp, ffname);
      free(cp);
      free(pl);
    }
    free(mfe_struc);
    free_alipf_arrays();
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

PRIVATE void mark_endgaps(char *seq, char egap) {
  int i,n;
  n = strlen(seq);
  for (i=0; i<n && (seq[i]=='-'); i++) {
    seq[i] = egap;
  }
  for (i=n-1; i>0 && (seq[i]=='-'); i--) {
    seq[i] = egap;
  }
}

PRIVATE void print_pi(const pair_info pi, FILE *file) {
  const char *pname[8] = {"","CG","GC","GU","UG","AU","UA", "--"};
  int i;

  /* numbering starts with 1 in output */
  fprintf(file, "%5d %5d %2d %5.1f%% %7.3f",
	  pi.i, pi.j, pi.bp[0], 100.*pi.p, pi.ent);
  for (i=1; i<=7; i++)
    if (pi.bp[i]) fprintf(file, " %s:%-4d", pname[i], pi.bp[i]);
  if (!pi.comp) fprintf(file, " +");
  fprintf(file, "\n");
}

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

/*-------------------------------------------------------------------------*/

PRIVATE char **annote(const char *structure, const char *AS[]) {
  /* produce annotation for colored drawings from PS_rna_plot_a() */
  char *ps, *colorps, **A;
  int i, n, s, pairings, maxl;
  short *ptable;
  char * colorMatrix[6][3] = {
    {"0.0 1", "0.0 0.6",  "0.0 0.2"},  /* red    */
    {"0.16 1","0.16 0.6", "0.16 0.2"}, /* ochre  */
    {"0.32 1","0.32 0.6", "0.32 0.2"}, /* turquoise */
    {"0.48 1","0.48 0.6", "0.48 0.2"}, /* green  */
    {"0.65 1","0.65 0.6", "0.65 0.2"}, /* blue   */
    {"0.81 1","0.81 0.6", "0.81 0.2"}  /* violet */
  };

  n = strlen(AS[0]);
  maxl = 1024;

  A = (char **) space(sizeof(char *)*2);
  ps = (char *) space(maxl);
  colorps = (char *) space(maxl);
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
    for (pairings=0,s=1; s<=7; s++) {
      if (pfreq[s]) pairings++;
    }

    if ((maxl - strlen(ps) < 192) || ((maxl - strlen(colorps)) < 64)) {
      maxl *= 2;
      ps = realloc(ps, maxl);
      colorps = realloc(colorps, maxl);
      if ((ps==NULL) || (colorps == NULL))
	  nrerror("out of memory in realloc");
    }

    if (pfreq[0]<=2 && pairings>0) {
      snprintf(pps, 64, "%d %d %s colorpair\n",
	       i,j, colorMatrix[pairings-1][pfreq[0]]);
      strcat(colorps, pps);
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
  A[0]=colorps;
  A[1]=ps;
  return A;
}

/*-------------------------------------------------------------------------*/

#define PMIN 0.0008
PRIVATE int compare_pair_info(const void *pi1, const void *pi2) {
  pair_info *p1, *p2;
  int  i, nc1, nc2;
  p1 = (pair_info *)pi1;  p2 = (pair_info *)pi2;
  for (nc1=nc2=0, i=1; i<=6; i++) {
    if (p1->bp[i]>0) nc1++;
    if (p2->bp[i]>0) nc2++;
  }
  /* sort mostly by probability, add
     epsilon * comp_mutations/(non-compatible+1) to break ties */
  return (p1->p + 0.01*nc1/(p1->bp[0]+1.)) <
	 (p2->p + 0.01*nc2/(p2->bp[0]+1.)) ? 1 : -1;
}

PRIVATE void print_aliout(char **AS, plist *pl, int n_seq, char * mfe, FILE *aliout) {
  int i, j, k, n, num_p=0, max_p = 64;
  pair_info *pi;
  double *duck, p;
  short *ptable;
  for (n=0; pl[n].i>0; n++);

  max_p = 64; pi = space(max_p*sizeof(pair_info));
  duck =  (double *) space((strlen(mfe)+1)*sizeof(double));
  ptable = make_pair_table(mfe);

  for (k=0; k<n; k++) {
    int s, type;
    p = pl[k].p; i=pl[k].i; j = pl[k].j;
    duck[i] -=  p * log(p);
    duck[j] -=  p * log(p);

    if (p<PMIN) continue;

    pi[num_p].i = i;
    pi[num_p].j = j;
    pi[num_p].p = p;
    pi[num_p].ent =  duck[i]+duck[j]-p*log(p);
    for (type=0; type<8; type++) pi[num_p].bp[type]=0;
    for (s=0; s<n_seq; s++) {
      int a,b;
      a=encode_char(toupper(AS[s][i-1]));
      b=encode_char(toupper(AS[s][j-1]));
      type = pair[a][b];
      if ((AS[s][i-1] == '-')||(AS[s][j-1] == '-')) type = 7;
      if ((AS[s][i-1] == '~')||(AS[s][j-1] == '~')) type = 7;
      pi[num_p].bp[type]++;
      pi[num_p].comp = (ptable[i] == j) ? 1:0;
    }
    num_p++;
    if (num_p>=max_p) {
      max_p *= 2;
      pi = xrealloc(pi, max_p * sizeof(pair_info));
    }
  }
  free(duck);
  pi[num_p].i=0;
  qsort(pi, num_p, sizeof(pair_info), compare_pair_info );

  /* print it */
  fprintf(aliout, "%d sequence; length of alignment %d\n",
	  n_seq, (int) strlen(AS[0]));
  fprintf(aliout, "alifold output\n");

  for (k=0; pi[k].i>0; k++) {
    pi[k].comp = (ptable[pi[k].i] == pi[k].j) ? 1:0;
    print_pi(pi[k], aliout);
  }
  fprintf(aliout, "%s\n", mfe);
  free(ptable);
  free(pi);
}


PRIVATE cpair *make_color_pinfo(char **sequences, plist *pl, int n_seq, bondT *mfe) {
  /* produce info for PS_color_dot_plot */
  cpair *cp;
  int i, n,s, a, b,z,t,j, c;
  int pfreq[7];
  for (n=0; pl[n].i>0; n++);
  c=0;
  cp = (cpair *) space(sizeof(cpair)*(n+1));
  for (i=0; i<n; i++) {
    int ncomp=0;
    if(pl[i].p>PMIN) {
      cp[c].i = pl[i].i;
      cp[c].j = pl[i].j;
      cp[c].p = pl[i].p;
      for (z=0; z<7; z++) pfreq[z]=0;
      for (s=0; s<n_seq; s++) {
	a=encode_char(toupper(sequences[s][cp[c].i-1]));
	b=encode_char(toupper(sequences[s][cp[c].j-1]));
	if ((sequences[s][cp[c].j-1]=='~')||(sequences[s][cp[c].i-1] == '~')) continue;
	pfreq[pair[a][b]]++;
      }
      for (z=1; z<7; z++) {
	if (pfreq[z]>0) {
	  ncomp++;
	}}
      cp[c].hue = (ncomp-1.0)/6.2;   /* hue<6/6.9 (hue=1 ==  hue=0) */
      cp[c].sat = 1 - MIN2( 1.0, (float) (pfreq[0]*2./*pi[i].bp[0]*//(n_seq)));
      c++;
    }
  }
  for (t=1; t<=mfe[0].i; t++) {
    int nofound=1;
      for (j=0; j<c; j++) {
	if ((cp[j].i==mfe[t].i)&&(cp[j].j==mfe[t].j)) {
	  cp[j].mfe=1;
	  nofound=0;
	  break;
	}
      }
      if(nofound) {
	fprintf(stderr,"mfe base pair with very low prob in pf: %d %d\n",mfe[t].i,mfe[t].j);
	cp = (cpair *) realloc(cp,sizeof(cpair)*(c+1));
	cp[c].i = mfe[t].i;
	cp[c].j = mfe[t].j;
	cp[c].p = 0.;
	cp[c].mfe=1;
	c++;
      }
    }
  return cp;
}

/*-------------------------------------------------------------------------*/

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAalifold [-cv float] [-nc float] [-E] [-old] [-r] [-R ribosum]\n"
	  "        [-mis] [-aln] [-color] [-circ] [-s num] [-se num]\n"
	  "        [-p[0]] [-C] [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n"
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]"
	  );
}
