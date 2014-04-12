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
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "RNAalifold_cmdl.h"

/*@unused@*/
static const char rcsid[] = "$Id: RNAalifold.c,v 1.23 2009/02/24 14:21:26 ivo Exp $";

#define MAX_NUM_NAMES    500

PRIVATE char  **annote(const char *structure, const char *AS[]);
PRIVATE void  print_pi(const pair_info pi, FILE *file);
PRIVATE void  print_aliout(vrna_fold_compound *vc, plist *pl, double threshold, char * mfe, FILE *aliout);
PRIVATE void  mark_endgaps(char *seq, char egap);
PRIVATE cpair *make_color_pinfo(char **sequences, plist *pl, double threshold, int n_seq, plist *mfel);

PRIVATE void
add_shape_constraints(vrna_fold_compound *vc,
                      const char *shape_method,
                      const char **shape_files,
                      const int *shape_file_association,
                      int verbose,
                      unsigned int constraint_type){

  float p1, p2;
  char method;

  if(!parse_soft_constraints_shape_method(shape_method, &method, &p1, &p2)){
    warn_user("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if(verbose){
    fprintf(stderr, "Using SHAPE method '%c'", method);
    if(method != 'W'){
      if(method == 'C')
        fprintf(stderr, " with parameter p1=%f", p1);
      else
        fprintf(stderr, " with parameters p1=%f and p2=%f", p1, p2);
    }
    fputc('\n', stderr);
  }

  if(method == 'M'){
    vrna_sc_add_mathews_ali(vc, shape_files, shape_file_association, p1, p2, constraint_type);
    return;
  }
}

/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[]){
  struct        RNAalifold_args_info args_info;
  unsigned int  input_type;
  char          ffname[FILENAME_MAX_LENGTH], gfname[FILENAME_MAX_LENGTH], fname[FILENAME_MAX_LENGTH];
  char          *input_string, *string, *structure, *cstruc, *ParamFile, *ns_bases, *c;
  int           s, n_seq, i, length, sym, noPS, with_shapes, verbose;
  int           endgaps, mis, circular, doAlnPS, doColor, doMEA, n_back, eval_energy, pf, istty;
  double        min_en, real_en, sfact, MEAgamma, bppmThreshold, betaScale;
  char          *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char          *names[MAX_NUM_NAMES];       /* sequence names */
  char          **shape_files, *shape_method;
  int           *shape_file_association;
  FILE          *clust_file = stdin;
  pf_paramT     *pf_parameters;
  model_detailsT  md;

  fname[0] = ffname[0] = gfname[0] = '\0';
  string = structure = cstruc = ParamFile = ns_bases = NULL;
  pf_parameters = NULL;
  endgaps = mis = pf = circular = doAlnPS = doColor = n_back = eval_energy = oldAliEn = doMEA = ribo = noPS = 0;
  do_backtrack  = 1;
  dangles       = 2;
  gquad         = 0;
  sfact         = 1.07;
  bppmThreshold = 1e-6;
  MEAgamma      = 1.0;
  betaScale     = 1.;
  shape_files   = NULL;
  shape_file_association = 0;
  shape_method  = NULL;
  with_shapes   = 0;
  max_bp_span   = -1;
  verbose       = 0;

  set_model_details(&md);

  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(RNAalifold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given)
    fold_constrained = 1;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
      warn_user("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }
  /* do not allow weak pairs */
  if(args_info.noLP_given)
    md.noLP = noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)
    md.noGU = noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;
  /* gquadruplex support */
  if(args_info.gquad_given)
    md.gquad = gquad = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  /* set energy model */
  if(args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)
    sfact = args_info.pfScale_arg;
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)
    md.circ = circular = 1;
  /* do not produce postscript output */
  if(args_info.noPS_given)
    noPS = 1;
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != -1)
      md.compute_bpp = do_backtrack = args_info.partfunc_arg;
  }
  /* MEA (maximum expected accuracy) settings */
  if(args_info.MEA_given){
    pf = doMEA = 1;
    if(args_info.MEA_arg != -1)
      MEAgamma = args_info.MEA_arg;
  }
  if(args_info.betaScale_given)
    betaScale = args_info.betaScale_arg;
  /* set the bppm threshold for the dotplot */
  if(args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0.,args_info.bppmThreshold_arg));
  /* set cfactor */
  if(args_info.cfactor_given)
    md.cv_fact = cv_fact = args_info.cfactor_arg;
  /* set nfactor */
  if(args_info.nfactor_given)
    md.nc_fact = nc_fact = args_info.nfactor_arg;
  if(args_info.endgaps_given)
    endgaps = 1;
  if(args_info.mis_given)
    mis = 1;
  if(args_info.color_given)
    doColor=1;
  if(args_info.aln_given)
    doAlnPS=1;
  if(args_info.old_given)
    md.oldAliEn = oldAliEn = 1;
  if(args_info.stochBT_given){
    n_back = args_info.stochBT_arg;
    md.compute_bpp = do_backtrack = 0;
    pf = 1;
    init_rand();
  }
  if(args_info.stochBT_en_given){
    n_back = args_info.stochBT_en_arg;
    md.compute_bpp = do_backtrack = 0;
    pf = 1;
    eval_energy = 1;
    init_rand();
  }
  if(args_info.ribosum_file_given){
    RibosumFile = strdup(args_info.ribosum_file_arg);
    md.ribo = ribo = 1;
  }
  if(args_info.ribosum_scoring_given){
    RibosumFile = NULL;
    md.ribo = ribo = 1;
  }
  if(args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  if(args_info.maxBPspan_given){
    md.max_bp_span = max_bp_span = args_info.maxBPspan_arg;
  }

  if(args_info.verbose_given){
    verbose = 1;
  }

  /* SHAPE reactivity data */
  if(args_info.shape_given){
    if(verbose)
      fprintf(stderr, "SHAPE reactivity data correction activated\n");

    with_shapes             = 1;
    shape_files             = (char **)space(sizeof(char*) * (args_info.shape_given + 1));
    shape_file_association  = (int *)space(sizeof(int*) * (args_info.shape_given + 1));

    /* find longest string in argument list */
    unsigned int longest_string = 0;
    for(s = 0; s < args_info.shape_given; s++)
      if(strlen(args_info.shape_arg[s]) > longest_string)
        longest_string = strlen(args_info.shape_arg[s]);

    char *tmp_string  = (char *)space(sizeof(char) * (longest_string + 1));
    int   tmp_number  = 0;

    for(s = 0; s < args_info.shape_given; s++){
      /* check whether we have int=string style that specifies a SHAPE file for a certain sequence number in the alignment */
      if(sscanf(args_info.shape_arg[s], "%d=%s", &tmp_number, tmp_string) == 2){
        shape_files[s]            = strdup(tmp_string);
        shape_file_association[s] = tmp_number - 1;
      } else {
        shape_files[s] = strdup(args_info.shape_arg[s]);
        shape_file_association[s] = s;
      }
      if(verbose)
        fprintf(stderr, "using SHAPE reactivity data provided in file %s for sequence %d\n", shape_files[s], shape_file_association[s]+1);
    }
    
    shape_file_association[s] = -1;

    free(tmp_string);
  }

  if(args_info.shapeMethod_given){
    shape_method = strdup(args_info.shapeMethod_arg);
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
  if(circular && gquad){
    nrerror("G-Quadruplex support is currently not available for circular RNA structures");
  }

  make_pair_matrix(); /* for make_color_pinfo */

  if (circular && noLonelyPairs)
    warn_user("depending on the origin of the circular sequence, "
              "some structures may be missed when using --noLP\n"
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
    print_tty_input_seq_str("Input aligned sequences in clustalw or stockholm format\n(enter a line starting with \"//\" to indicate the end of your input)");

  n_seq = read_clustal(clust_file, AS, names);
  if (n_seq==0) nrerror("no sequences found");

  if(with_shapes){
    
    if(s != n_seq)
      warn_user("number of sequences in alignment does not match number of provided SHAPE reactivity data files! ");

    shape_files             = (char **)xrealloc(shape_files, (n_seq + 1) * sizeof(char *));
    shape_file_association  = (int *)xrealloc(shape_file_association, (n_seq + 1) * sizeof(int));

  }

  if (clust_file != stdin) fclose(clust_file);
  /*
  ########################################################
  # done with 'stdin' handling, now init everything properly
  ########################################################
  */

  length    = (int)   strlen(AS[0]);
  structure = (char *)space((unsigned)length + 1);

  if(fold_constrained && cstruc != NULL)
    strncpy(structure, cstruc, length);

  if (endgaps)
    for (i=0; i<n_seq; i++) mark_endgaps(AS[i], '~');

  /*
  ########################################################
  # begin actual calculations
  ########################################################
  */
  unsigned int options = VRNA_CONSTRAINT_SOFT_MFE;

  if(pf)
    options |= VRNA_CONSTRAINT_SOFT_PF;

  vrna_fold_compound *vc = vrna_get_fold_compound_ali((const char **)AS, &md, VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));

  if(fold_constrained){
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_PIPE
                          | VRNA_CONSTRAINT_DOT
                          | VRNA_CONSTRAINT_X
                          | VRNA_CONSTRAINT_ANG_BRACK
                          | VRNA_CONSTRAINT_RND_BRACK;

    vrna_hc_add(vc, (const char *)structure, constraint_options);
  }

  if(with_shapes)
    add_shape_constraints(vc, \
                          shape_method, \
                          (const char **)shape_files, \
                          shape_file_association, \
                          verbose, \
                          VRNA_CONSTRAINT_SOFT_MFE | ((pf) ? VRNA_CONSTRAINT_SOFT_PF : 0));

  min_en = vrna_ali_fold(vc, structure);

  if(md.circ){
    int     i;
    double  s = 0;
    for (i=0; AS[i]!=NULL; i++)
      s += vrna_eval_structure(AS[i], structure, vc->params);
    real_en   = s/i;
  } else {

    float *ens    = (float *)space(2*sizeof(float));

    if(md.gquad)
      energy_of_ali_gquad_structure((const char **)AS, structure, n_seq, ens);
    else
      energy_of_alistruct((const char **)AS, structure, n_seq, ens);

    real_en       = ens[0];

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

    if(md.gquad){
      if (doColor)
        (void) PS_rna_plot_a_gquad(string, structure, ffname, A[0], A[1]);
      else
        (void) PS_rna_plot_a_gquad(string, structure, ffname, NULL, A[1]);
    } else {
      if (doColor)
        (void) PS_rna_plot_a(string, structure, ffname, A[0], A[1]);
      else
        (void) PS_rna_plot_a(string, structure, ffname, NULL, A[1]);
    }
    free(A[0]); free(A[1]); free(A);
  }
  if (doAlnPS)
    PS_color_aln(structure, "aln.ps", (const char **) AS, (const char **) names);

  /* free mfe arrays */
  destroy_mfe_matrices(vc->matrices);
  vc->matrices = NULL;

  if (pf) {
    float energy, kT;
    char * mfe_struc;

    mfe_struc = strdup(structure);

    kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
    pf_scale = exp(-(sfact*min_en)/kT/length);
    if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
    fflush(stdout);

    if (cstruc!=NULL) strncpy(structure, cstruc, length+1);

    /* rescale energy parameters according to above calculated pf_scale */
    pf_parameters = get_boltzmann_factors_ali(n_seq, temperature, betaScale, md, pf_scale);

    /* change energy parameters in vc */
    vrna_update_pf_params(vc, pf_parameters);

    energy = vrna_ali_pf_fold(vc, structure, NULL);

    if (n_back>0) {
      /*stochastic sampling*/
      for (i=0; i<n_back; i++) {
        char *s;
        double prob=1.;
        s = vrna_ali_pbacktrack(vc, &prob);
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
      plist *pl, *mfel;

      pl    = vrna_get_plist_from_pr(vc, bppmThreshold);
      mfel  = vrna_get_plist_from_db(mfe_struc, 0.95*0.95);

      if (!circular){
        float *ens;
        cent = vrna_get_centroid_struct(vc, &dist);
        ens=(float *)space(2*sizeof(float));
        energy_of_alistruct((const char **)AS, cent, n_seq, ens);

        printf("%s %6.2f {%6.2f + %6.2f}\n",cent,ens[0]-ens[1],ens[0],(-1)*ens[1]);
        free(cent);
        free(ens);
      }
      if(doMEA){
        float mea, *ens;
        plist *pl2;
        pl2 = vrna_get_plist_from_pr(vc, 1e-4/(1+MEAgamma));
        mea = MEA(pl2, structure, MEAgamma);
        ens = (float *)space(2*sizeof(float));
        if(circular)
          energy_of_alistruct((const char **)AS, structure, n_seq, ens);
        else
          ens[0] = vrna_eval_structure(string, structure, vc->params);
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
        print_aliout(vc, pl, bppmThreshold, mfe_struc, aliout);
      }
      fclose(aliout);
      if (fname[0]!='\0') {
        strcpy(ffname, fname);
        strcat(ffname, "_dp.ps");
      } else strcpy(ffname, "alidot.ps");
      cp = make_color_pinfo(AS,pl, bppmThreshold, n_seq, mfel);
      (void) PS_color_dot_plot(string, cp, ffname);
      free(cp);
      free(pl);
      free(mfel);
    }
    free(mfe_struc);
    destroy_fold_compound(vc);
    free(pf_parameters);
  }
  if (cstruc!=NULL) free(cstruc);
  (void) fflush(stdout);
  if(shape_files)
    free(shape_files);
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
  ptable = vrna_pt_get(structure);
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

PRIVATE void
print_aliout( vrna_fold_compound *vc,
              plist *pl,
              double threshold,
              char *mfe,
              FILE *aliout){

  int k;
  pair_info *pi;
  char  **AS    = vc->sequences;
  int   n_seq   = vc->n_seq;

  pi = vrna_ali_get_pair_info(vc, (const char *)mfe, threshold);

  /* print it */
  fprintf(aliout, "%d sequence; length of alignment %d\n",
          n_seq, (int) strlen(AS[0]));
  fprintf(aliout, "alifold output\n");

  for (k=0; pi[k].i>0; k++)
    print_pi(pi[k], aliout);

  fprintf(aliout, "%s\n", mfe);
  free(pi);
}


PRIVATE cpair *make_color_pinfo(char **sequences, plist *pl, double threshold, int n_seq, plist *mfel) {
  /* produce info for PS_color_dot_plot */
  cpair *cp;
  int i, n,s, a, b,z,t,j, c;
  int pfreq[7];
  for (n=0; pl[n].i>0; n++);
  c=0;
  cp = (cpair *) space(sizeof(cpair)*(n+1));
  for (i=0; i<n; i++) {
    int ncomp=0;
    if(pl[i].p>threshold) {
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
      cp[c].sat = 1 - MIN2( 1.0, (float) (pfreq[0]*2. /*pi[i].bp[0]*/ /(n_seq)));
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
