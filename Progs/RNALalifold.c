/* Last changed Time-stamp: <2006-03-02 22:48:15 ivo> */
/*
                  Local version of RNAalifold

                  c Ivo L Hofacker, Stephan Bernhart
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
#include "Lfold.h"
#include "aln_util.h"
#include "read_epars.h"
#include "RNALalifold_cmdl.h"

/*@unused@*/
static const char rcsid[] = "$Id: RNALalifold.c,v 1.1 2007/06/23 09:52:29 ivo Exp $";

/*@exits@*/
PRIVATE void  usage(void);
PRIVATE char  *annote(const char *structure, const char *AS[]);
PRIVATE void  print_pi(const pair_info pi, FILE *file);
PRIVATE cpair *make_color_pinfo(const pair_info *pi);
PRIVATE cpair *make_color_pinfo2(char **sequences, plist *pl, int n_seq);

#define MAX_NUM_NAMES    500

int main(int argc, char *argv[]){
  struct        RNALalifold_args_info args_info;
  char          *string, *structure, *ParamFile, *ns_bases, *c;
  char          ffname[FILENAME_MAX_LENGTH], gfname[FILENAME_MAX_LENGTH], fname[FILENAME_MAX_LENGTH];
  int           n_seq, i, length, sym, r, maxdist, unchangednc, unchangedcv;
  int           mis, pf, istty;
  float         cutoff;
  double        min_en, real_en, sfact;
  char          *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char          *names[MAX_NUM_NAMES];       /* sequence names */
  FILE          *clust_file = stdin;

  string = structure = ParamFile = ns_bases = NULL;
  mis = pf      = 0;
  maxdist       = 70;
  do_backtrack  = unchangednc = unchangedcv = 1;
  dangles       = 2;
  sfact         = 1.07;
  cutoff        = 0.0005;
  ribo          = 0;
  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNALalifold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      warn_user("required dangle model not implemented, falling back to default dangles=2");
    else
      dangles = args_info.dangles_arg;
  }
  /* do not allow weak pairs */
  if(args_info.noLP_given)        noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) no_closingGU = 1;
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)     sfact = args_info.pfScale_arg;
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != -1)
      do_backtrack = args_info.partfunc_arg;
  }
  /* set cfactor */
  if(args_info.cfactor_given){
    cv_fact = args_info.cfactor_arg;
    unchangedcv = 0;
  }
  /* set nfactor */
  if(args_info.nfactor_given){
    nc_fact = args_info.nfactor_arg;
    unchangednc = 0;
  }
  /* set the maximum base pair span */
  if(args_info.span_given)        maxdist = args_info.span_arg;
  /* set the pair probability cutoff */
  if(args_info.cutoff_given)      cutoff  = args_info.cutoff_arg;
  /* calculate most informative sequence */
  if(args_info.mis_given)         mis = 1;
  if(args_info.csv_given)         csv = 1;
  if(args_info.ribosum_file_given){
    RibosumFile = strdup(args_info.ribosum_file_arg);
    ribo = 1;
  }
  if(args_info.ribosum_scoring_given){
    RibosumFile = NULL;
    ribo = 1;
  }

  /* check unnamed options a.k.a. filename of input alignment */
  if(args_info.inputs_num == 1){
    clust_file = fopen(args_info.inputs[0], "r");
    if(clust_file == NULL){
      fprintf(stderr, "can't open %s\n", args_info.inputs[0]);
    }
  }
  else{
    RNALalifold_cmdline_parser_print_help();
    exit(1);
  }

  /* free allocated memory of command line data structure */
  RNALalifold_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  if ((ribo==1)&&(unchangednc)) nc_fact=0.5;
  if ((ribo==1)&&(unchangedcv)) cv_fact=0.6;

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

  if (istty && (clust_file == stdin)) {
    print_tty_input_seq_str("Input aligned sequences in clustalw format");
  }

  n_seq = read_clustal(clust_file, AS, names);
  if (clust_file != stdin) fclose(clust_file);
  if (n_seq==0)
    nrerror("no sequences found");

  length = (int) strlen(AS[0]);
  if (length<maxdist) {
    fprintf(stderr, "Alignment length < window size: setting L=%d\n",length);
    maxdist=length;
  }

  structure = (char *) space((unsigned) length+1);

  /*
  #############################################
  # begin calculations
  #############################################
  */
  update_fold_params();
  if(!pf)
    min_en = aliLfold((const char **) AS, structure, maxdist);
  {
    eos_debug=-1; /* shut off warnings about nonstandard pairs */
    /*   for (i=0; AS[i]!=NULL; i++)
    s += energy_of_struct(AS[i], structure);
    real_en = s/i;*/
  }
  string = (mis) ? consens_mis((const char **) AS) : consensus((const char **) AS);
  printf("%s\n%s\n", string, structure);
  /*  if (istty)
    printf("\n minimum free energy = %6.2f kcal/mol (%6.2f + %6.2f)\n",
           min_en, real_en, min_en - real_en);
  else
    printf(" (%6.2f = %6.2f + %6.2f) \n", min_en, real_en, min_en-real_en );
  */
  strcpy(ffname, "alirna.ps");
  strcpy(gfname, "alirna.g");

  /*  if (length<=2500) {
    char *A;
    A = annote(structure, (const char**) AS);
    (void) PS_rna_plot_a(string, structure, ffname, NULL, A);
    free(A);
  } else
    fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");
  */
  /* {*/ /* free mfe arrays but preserve base_pair for PS_dot_plot */
  /*  struct bond  *bp;
    bp = base_pair; base_pair = space(16);
    free_alifold_arrays();  / * frees base_pair *  /
    base_pair = bp;
  }*/
  if (pf) {
    double energy, kT;
    plist *pl;
    char * mfe_struc;

    mfe_struc = strdup(structure);

    kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
    pf_scale = -1;/*exp(-(sfact*min_en)/kT/length);*/
    if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
    fflush(stdout);

    /* init_alipf_fold(length); */

    /* energy = alipfW_fold(AS, structure, &pl, maxdist, cutoff); */

    if (do_backtrack) {
      printf("%s", structure);
      /*if (!istty) printf(" [%6.2f]\n", energy);
        else */
      printf("\n");
    }
    /*if ((istty)||(!do_backtrack))
      printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
    useless!!*/
    /* printf(" frequency of mfe structure in ensemble %g\n",
       exp((energy-min_en)/kT));*/

    if (do_backtrack) {
      FILE *aliout;
      cpair *cp;
      strcpy(ffname, "alifold.out");
      aliout = fopen(ffname, "w");
      if (!aliout) {
        fprintf(stderr, "can't open %s    skipping output\n", ffname);
      } else {
        fprintf(aliout, "%d sequence; length of alignment %d\n",
                n_seq, length);
        fprintf(aliout, "alifold output\n");

        fprintf(aliout, "%s\n", structure);
      }
      strcpy(ffname, "alidotL.ps");
      cp = make_color_pinfo2(AS,pl,n_seq);
      (void) PS_color_dot_plot_turn(string, cp, ffname, maxdist);
      free(cp);
    }
    free(mfe_struc);
    free(pl);
  }
  free(base_pair);
  (void) fflush(stdout);
  free(string);
  free(structure);
  for (i=0; AS[i]; i++) {
    free(AS[i]); free(names[i]);
  }
  return 0;
}

PRIVATE void print_pi(const pair_info pi, FILE *file) {
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


#if 0
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
#endif
/*-------------------------------------------------------------------------*/

PRIVATE cpair *make_color_pinfo2(char **sequences, plist *pl, int n_seq) {
  cpair *cp;
  int i, n,s, a, b,z;
  int franz[7];
  for (n=0; pl[n].i>0; n++);
  cp = (cpair *) space(sizeof(cpair)*(n+1));
  for (i=0; i<n; i++) {
    int ncomp=0;
    cp[i].i = pl[i].i;
    cp[i].j = pl[i].j;
    cp[i].p = pl[i].p;
    for (z=0; z<7; z++) franz[z]=0;
    for (s=0; s<n_seq; s++) {
      a=encode_char(toupper(sequences[s][cp[i].i-1]));
      b=encode_char(toupper(sequences[s][cp[i].j-1]));
      if ((sequences[s][cp[i].j-1]=='~')||(sequences[s][cp[i].i-1] == '~')) continue;
      franz[BP_pair[a][b]]++;
    }
    for (z=1; z<7; z++) {
      if (franz[z]>0) {
        ncomp++;
      }}
    cp[i].hue = (ncomp-1.0)/6.2;   /* hue<6/6.9 (hue=1 ==  hue=0) */
    cp[i].sat = 1 - MIN2( 1.0, franz[0]/*pi[i].bp[0]*//2.5);
    /*computation of entropy is sth for the ivo*/
    /* cp[i].mfe = pi[i].comp;  don't have that .. yet*/
  }
  return cp;
}
