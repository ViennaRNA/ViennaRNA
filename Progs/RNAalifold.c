/* Last changed Time-stamp: <2009-02-24 14:49:24 ivo> */
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
#include "read_epars.h"
#include "MEA.h"
#include "params.h"
#include "RNAalifold_cmdl.h"

/*@unused@*/
static const char rcsid[] = "$Id: RNAalifold.c,v 1.23 2009/02/24 14:21:26 ivo Exp $";

#define MAX_NUM_NAMES    500

PRIVATE char  **annote(const char *structure, const char *AS[]);
PRIVATE void  print_pi(const pair_info pi, FILE *file);
PRIVATE void  print_aliout(char **AS, plist *pl, int n_seq, char * mfe, FILE *aliout);
PRIVATE void  mark_endgaps(char *seq, char egap);
PRIVATE cpair *make_color_pinfo(char **sequences, plist *pl, int n_seq, plist *mfel);

/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[]){
  struct        RNAalifold_args_info args_info;
  unsigned int  input_type;
  char          ffname[FILENAME_MAX_LENGTH], gfname[FILENAME_MAX_LENGTH], fname[FILENAME_MAX_LENGTH];
  char          *input_string, *string, *structure, *cstruc, *ParamFile, *ns_bases, *c;
  int           n_seq, i, length, sym, r, noPS;
  int           endgaps, mis, circular, doAlnPS, doColor, doMEA, n_back, eval_energy, pf, istty;
  double        min_en, real_en, sfact, MEAgamma, bppmThreshold, betaScale;
  char          *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char          *names[MAX_NUM_NAMES];       /* sequence names */
  FILE          *clust_file = stdin;
  pf_paramT     *pf_parameters;
  model_detailsT  md;

  fname[0] = ffname[0] = gfname[0] = '\0';
  string = structure = cstruc = ParamFile = ns_bases = NULL;
  pf_parameters = NULL;
  endgaps = mis = pf = circular = doAlnPS = doColor = n_back = eval_energy = oldAliEn = doMEA = ribo = noPS = 0;
  do_backtrack  = 1;
  dangles       = 2;
  sfact         = 1.07;
  bppmThreshold = 1e-6;
  MEAgamma      = 1.0;
  betaScale     = 1.;

  set_model_details(&md);

  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(RNAalifold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given)  fold_constrained=1;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     md.special_hp = tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)     md.dangles = dangles = args_info.dangles_arg;
  /* do not allow weak pairs */
  if(args_info.noLP_given)        md.noLP = noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        md.noGU = noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) md.noGUclosure = no_closingGU = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)     sfact = args_info.pfScale_arg;
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)        circular=1;
  /* do not produce postscript output */
  if(args_info.noPS_given)        noPS = 1;
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != -1)
      do_backtrack = args_info.partfunc_arg;
  }
  /* MEA (maximum expected accuracy) settings */
  if(args_info.MEA_given){
    pf = doMEA = 1;
    if(args_info.MEA_arg != -1)
      MEAgamma = args_info.MEA_arg;
  }
  if(args_info.betaScale_given)   betaScale = args_info.betaScale_arg;
  /* set the bppm threshold for the dotplot */
  if(args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0.,args_info.bppmThreshold_arg));
  /* set cfactor */
  if(args_info.cfactor_given)     cv_fact = args_info.cfactor_arg;
  /* set nfactor */
  if(args_info.nfactor_given)     nc_fact = args_info.nfactor_arg;
  if(args_info.endgaps_given)     endgaps = 1;
  if(args_info.mis_given)         mis = 1;
  if(args_info.color_given)       doColor=1;
  if(args_info.aln_given)         doAlnPS=1;
  if(args_info.old_given)         oldAliEn = 1;
  if(args_info.stochBT_given){
    n_back = args_info.stochBT_arg;
    do_backtrack = 0;
    pf = 1;
    init_rand();
  }
  if(args_info.stochBT_en_given){
    n_back = args_info.stochBT_en_arg;
    do_backtrack = 0;
    pf = 1;
    eval_energy = 1;
    init_rand();
  }
  if(args_info.ribosum_file_given){
    RibosumFile = strdup(args_info.ribosum_file_arg);
    ribo = 1;
  }
  if(args_info.ribosum_scoring_given){
    RibosumFile = NULL;
    ribo = 1;
  }
  /* alignment file name given as unnamed option? */
  if(args_info.inputs_num == 1){
    clust_file = fopen(args_info.inputs[0], "r");
    if (clust_file == NULL) {
      fprintf(stderr, "can't open %s\n", args_info.inputs[0]);
    }
  }

  /* free allocated memory of command line data structure */
  RNAalifold_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  make_pair_matrix();

  if (circular && noLonelyPairs)
    warn_user("depending on the origin of the circular sequence, "
            "some structures may be missed when using -noLP\n"
            "Try rotating your sequence a few times\n");

  if (ParamFile != NULL) read_parameter_file(ParamFile);

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

  /*
  ########################################################
  # handle user input from 'stdin' if necessary
  ########################################################
  */
  if(fold_constrained){
    if(istty){
      print_tty_constraint_full();
      print_tty_input_seq_str("");
    }
    input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
    if(input_type & VRNA_INPUT_QUIT){ return 0;}
    else if((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)){
      cstruc = strdup(input_string);
      free(input_string);
    }
    else warn_user("constraints missing");
  }

  if (istty && (clust_file == stdin))
    print_tty_input_seq_str("Input aligned sequences in clustalw format");

  n_seq = read_clustal(clust_file, AS, names);
  if (n_seq==0) nrerror("no sequences found");

  if (clust_file != stdin) fclose(clust_file);
  /*
  ########################################################
  # done with 'stdin' handling, now init everything properly
  ########################################################
  */

  length    = (int)   strlen(AS[0]);
  structure = (char *)space((unsigned) length+1);

  if(fold_constrained && cstruc != NULL)
    strncpy(structure, cstruc, length);

  if (endgaps)
    for (i=0; i<n_seq; i++) mark_endgaps(AS[i], '~');

  /*
  ########################################################
  # begin actual calculations
  ########################################################
  */

  if (circular) {
    int     i;
    double  s = 0;
    min_en    = circalifold((const char **)AS, structure);
    for (i=0; AS[i]!=NULL; i++)
      s += energy_of_circ_structure(AS[i], structure, -1);
    real_en = s/i;
  } else {
    float *ens  = (float *)space(2*sizeof(float));
    min_en      = alifold((const char **)AS, structure);
    energy_of_alistruct((const char **)AS, structure, n_seq, ens);
    real_en     = ens[0];
    free(ens);
  }

  string = (mis) ? consens_mis((const char **) AS) : consensus((const char **) AS);
  printf("%s\n%s", string, structure);

  if (istty)
    printf("\n minimum free energy = %6.2f kcal/mol (%6.2f + %6.2f)\n",
           min_en, real_en, min_en - real_en);
  else
    printf(" (%6.2f = %6.2f + %6.2f) \n", min_en, real_en, min_en-real_en );

  strcpy(ffname, "alirna.ps");
  strcpy(gfname, "alirna.g");

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

  /* free mfe arrays */
  free_alifold_arrays();

  if (pf) {
    float energy, kT;
    char * mfe_struc;

    mfe_struc = strdup(structure);

    kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
    pf_scale = exp(-(sfact*min_en)/kT/length);
    if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
    fflush(stdout);

    if (cstruc!=NULL) strncpy(structure, cstruc, length+1);

    pf_parameters = get_boltzmann_factors_ali(n_seq, temperature, betaScale, md, pf_scale);
    energy = alipf_fold_par((const char **)AS, structure, NULL, pf_parameters, do_backtrack, fold_constrained, circular);

    if (n_back>0) {
      /*stochastic sampling*/
      for (i=0; i<n_back; i++) {
        char *s;
        double prob=1.;
        s = alipbacktrack(&prob);
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
      FLT_OR_DBL *probs = export_ali_bppm();
      plist *pl, *mfel;

      assign_plist_from_pr(&pl, probs, length, bppmThreshold);
      assign_plist_from_db(&mfel, mfe_struc, 0.95*0.95);

      if (!circular){
        float *ens;
        cent = get_centroid_struct_pr(length, &dist, probs);
        ens=(float *)space(2*sizeof(float));
        energy_of_alistruct((const char **)AS, cent, n_seq, ens);
        /*cent_en = energy_of_struct(string, cent);*//*ali*/
        printf("%s %6.2f {%6.2f + %6.2f}\n",cent,ens[0]-ens[1],ens[0],(-1)*ens[1]);
        free(cent);
        free(ens);
      }
      if(doMEA){
        float mea, *ens;
        plist *pl2;
        assign_plist_from_pr(&pl2, probs, length, 1e-4/(1+MEAgamma));
        mea = MEA(pl2, structure, MEAgamma);
        ens = (float *)space(2*sizeof(float));
        if(circular)
          energy_of_alistruct((const char **)AS, structure, n_seq, ens);
        else
          ens[0] = energy_of_structure(string, structure, 0);
        printf("%s {%6.2f MEA=%.2f}\n", structure, ens[0], mea);
        free(ens);
        free(pl2);
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
      cp = make_color_pinfo(AS,pl, n_seq, mfel);
      (void) PS_color_dot_plot(string, cp, ffname);
      free(cp);
      free(pl);
    }
    free(mfe_struc);
    free_alipf_arrays();
    free(pf_parameters);
  }
  if (cstruc!=NULL) free(cstruc);
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


PRIVATE cpair *make_color_pinfo(char **sequences, plist *pl, int n_seq, plist *mfel) {
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
  for (t=0; mfel[t].i > 0; t++) {
    int nofound=1;
      for (j=0; j<c; j++) {
        if ((cp[j].i==mfel[t].i)&&(cp[j].j==mfel[t].j)) {
          cp[j].mfe=1;
          nofound=0;
          break;
        }
      }
      if(nofound) {
        fprintf(stderr,"mfe base pair with very low prob in pf: %d %d\n",mfel[t].i,mfel[t].j);
        cp = (cpair *) realloc(cp,sizeof(cpair)*(c+1));
        cp[c].i = mfel[t].i;
        cp[c].j = mfel[t].j;
        cp[c].p = 0.;
        cp[c].mfe=1;
        c++;
      }
    }
  return cp;
}
