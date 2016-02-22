/* Last changed Time-stamp: <2008-12-03 16:38:01 ivo> */
/*
                Ineractive access to suboptimal folding

                           c Ivo L Hofacker
                          Vienna RNA package
*/
#include <config.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/subopt.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/file_formats.h"
#include "RNAsubopt_cmdl.h"

/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAsubopt.c,v 1.20 2008/12/03 16:55:44 ivo Exp $";

static void
add_shape_constraints(vrna_fold_compound_t *vc,
                      const char *shape_method,
                      const char *shape_conversion,
                      const char *shape_file,
                      int verbose,
                      unsigned int constraint_type){

  float p1, p2;
  char method;
  char *sequence;
  double *values;
  int i, length = vc->length;

  if(!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)){
    vrna_message_warning("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if(verbose){
    fprintf(stderr, "Using SHAPE method '%c'", method);
    if(method != 'W'){
      if(method == 'Z')
        fprintf(stderr, " with parameter p1=%f", p1);
      else
        fprintf(stderr, " with parameters p1=%f and p2=%f", p1, p2);
    }
    fputc('\n', stderr);
  }

  sequence = vrna_alloc(sizeof(char) * (length + 1));
  values = vrna_alloc(sizeof(double) * (length + 1));
  vrna_file_SHAPE_read(shape_file, length, method == 'W' ? 0 : -1, sequence, values);

  if(method == 'D'){
    (void)vrna_sc_add_SHAPE_deigan(vc, (const double *)values, p1, p2, constraint_type);
  }
  else if(method == 'Z'){
    (void)vrna_sc_add_SHAPE_zarringhalam(vc, (const double *)values, p1, 0.5, shape_conversion, constraint_type);
  } else {
    assert(method == 'W');
    FLT_OR_DBL *v = vrna_alloc(sizeof(FLT_OR_DBL) * (length + 1));
    for(i = 0; i < length; i++)
      v[i] = values[i];

    vrna_sc_add_up(vc, v, constraint_type);

    free(v);
  }

  free(values);
  free(sequence);
}

PRIVATE void putoutzuker(vrna_subopt_solution_t* zukersolution);

int main(int argc, char *argv[]){
  struct          RNAsubopt_args_info args_info;
  unsigned int    input_type;
  unsigned int    rec_type, read_opt;
  char            fname[FILENAME_MAX_LENGTH], *c, *input_string, *rec_sequence, *rec_id, **rec_rest, *orig_sequence;
  char            *constraints_file, *cstruc, *structure, *ParamFile, *ns_bases, *shape_file, *shape_method, *shape_conversion;
  int             i, length, l, cl, sym, istty;
  double          deltaf, deltap, betaScale;
  int             delta, n_back, noconv, circular, dos, zuker, gquad, with_shapes, verbose, max_bp_span, enforceConstraints;
  vrna_md_t       md;

  do_backtrack  = 1;
  dangles       = 2;
  betaScale     = 1.;
  gquad         = 0;
  delta         = 100;
  deltap = n_back = noconv = circular = dos = zuker = 0;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  input_string  = c = cstruc = structure = ParamFile = ns_bases = NULL;
  constraints_file = NULL;
  shape_file    = NULL;
  shape_method  = NULL;
  with_shapes   = 0;
  verbose       = 0;
  max_bp_span   = -1;
  enforceConstraints  = 0;

  set_model_details(&md);

  /* switch on unique multibranch loop decomposition */
  md.uniq_ML = 1;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAsubopt_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        md.temperature = temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given){
    fold_constrained=1;
    if(args_info.constraint_arg[0] != '\0')
      constraints_file = strdup(args_info.constraint_arg);
  }
  /* enforce base pairs given in constraint string rather than weak enforce */
  if(args_info.enforceConstraint_given)   enforceConstraints = 1;
  /* SHAPE reactivity data */
  if(args_info.shape_given){
    with_shapes = 1;
    shape_file = strdup(args_info.shape_arg);
  }

  shape_method = strdup(args_info.shapeMethod_arg);
  shape_conversion = strdup(args_info.shapeConversion_arg);
  if(args_info.verbose_given){
    verbose = 1;
  }
  if(args_info.maxBPspan_given){
    md.max_bp_span = max_bp_span = args_info.maxBPspan_arg;
  }
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     md.special_hp = tetra_loop = 0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }
  /* do not allow weak pairs */
  if(args_info.noLP_given)        md.noLP = noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        md.noGU = noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) md.noGUclosure = no_closingGU = 1;
  /* gquadruplex support */
  if(args_info.gquad_given)       md.gquad = gquad = 1;
  /* enforce canonical base pairs in any case? */
  if(args_info.canonicalBPonly_given) md.canonicalBPonly = canonicalBPonly = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)      noconv = 1;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* energy range */
  if(args_info.deltaEnergy_given) delta = (int) (0.1+args_info.deltaEnergy_arg*100);
  /* energy range after post evaluation */
  if(args_info.deltaEnergyPost_given) deltap = args_info.deltaEnergyPost_arg;
  /* sorted output */
  if(args_info.sorted_given)      subopt_sorted = 1;
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)        md.circ = circular = 1;
  /* stochastic backtracking */
  if(args_info.stochBT_given){
    n_back = args_info.stochBT_arg;
    vrna_init_rand();
    md.compute_bpp = 0;
  }
  if(args_info.betaScale_given)   md.betaScale = betaScale = args_info.betaScale_arg;
  /* density of states */
  if(args_info.dos_given){
    dos = 1;
    print_energy = -999999;
  }
  /* logarithmic multiloop energies */
  if(args_info.logML_given) md.logML = logML = 1;
  /* zuker subopts */
  if(args_info.zuker_given) zuker = 1;

  if(zuker){
    if(circular){
      vrna_message_warning("Sorry, zuker subopts not yet implemented for circfold");
      RNAsubopt_cmdline_parser_print_help();
      exit(1);
    }
    else if(n_back>0){
      vrna_message_warning("Can't do zuker subopts and stochastic subopts at the same time");
      RNAsubopt_cmdline_parser_print_help();
      exit(1);
    }
    else if(gquad){
      vrna_message_warning("G-quadruplex support for Zuker subopts not implemented yet");
      RNAsubopt_cmdline_parser_print_help();
      exit(1);
    }
  }

  if(gquad && (n_back > 0)){
    vrna_message_warning("G-quadruplex support for stochastic backtracking not implemented yet");
    RNAsubopt_cmdline_parser_print_help();
    exit(1);
  }

  /* free allocated memory of command line data structure */
  RNAsubopt_cmdline_parser_free(&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */

  if (ParamFile != NULL) read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards = vrna_alloc(33);
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

  /* print user help if we get input from tty */
  if(istty){
    if(!zuker)
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
    if(fold_constrained){
      vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK);
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n");
    }
    else vrna_message_input_seq_simple();
  }

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if(istty)             read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  if(!fold_constrained) read_opt |= VRNA_INPUT_NO_REST;

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    if(rec_id){
      /* if(!istty) printf("%s\n", rec_id); */
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) vrna_seq_toRNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound_t *vc = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE | (circ ? 0 : VRNA_OPTION_HYBRID) | ((n_back > 0) ? VRNA_OPTION_PF : 0));
    length    = vc->length;

    structure = (char *) vrna_alloc((char) length+1);

    /* parse the rest of the current dataset to obtain a structure constraint */
    if(fold_constrained){
      if(constraints_file){
        vrna_constraints_add(vc, constraints_file, VRNA_OPTION_MFE | ((n_back > 0) ? VRNA_OPTION_PF : 0));
      } else {
        cstruc = NULL;
        int cp = -1;
        unsigned int coptions = (rec_id) ? VRNA_OPTION_MULTILINE : 0;
        cstruc = vrna_extract_record_rest_structure((const char **)rec_rest, 0, coptions);
        cstruc = vrna_cut_point_remove(cstruc, &cp);
        if(vc->cutpoint != cp){
          fprintf(stderr,"cut_point = %d cut = %d\n", vc->cutpoint, cp);
          vrna_message_error("Sequence and Structure have different cut points.");
        }
        cl = (cstruc) ? (int)strlen(cstruc) : 0;

        if(cl == 0)           vrna_message_warning("structure constraint is missing");
        else if(cl < length)  vrna_message_warning("structure constraint is shorter than sequence");
        else if(cl > length)  vrna_message_error("structure constraint is too long");

        if(cstruc){
          strncpy(structure, cstruc, sizeof(char)*(cl+1));

          /* convert pseudo-dot-bracket to actual hard constraints */
          unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;

          if(enforceConstraints)
            constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;
          vrna_constraints_add(vc, (const char *)structure, constraint_options);
        }
      }
    }

    if(with_shapes)
      add_shape_constraints(vc, shape_method, shape_conversion, shape_file, verbose, VRNA_OPTION_MFE | ((n_back > 0) ? VRNA_OPTION_PF : 0));

    if(istty){
      if (vc->cutpoint == -1)
        printf("length = %d\n", length);
      else
        printf("length1 = %d\nlength2 = %d\n", vc->cutpoint - 1, length - vc->cutpoint + 1);
    }

    /*
    ########################################################
    # begin actual computations
    ########################################################
    */

    if((logML != 0 || dangles==1 || dangles==3) && dos == 0)
      if(deltap<=0) deltap = delta/100. + 0.001;
    if (deltap>0)
      print_energy = deltap;

    /* stochastic backtracking */
    if(n_back>0){
      double mfe, kT;
      char *ss;
      if (fname[0] != '\0') printf(">%s\n", fname);

      ss = (char *) vrna_alloc(strlen(rec_sequence)+1);
      strncpy(ss, structure, length);
      mfe = vrna_mfe(vc, ss);
      /* rescale Boltzmann factors according to predicted MFE */
      vrna_exp_params_rescale(vc, &mfe);
      /* ignore return value, we are not interested in the free energy */
      (void)vrna_pf(vc, ss);

      free(ss);
      for (i=0; i<n_back; i++) {
        char *s;
        s = vrna_pbacktrack(vc);
        printf("%s\n", s);
        free(s);
      }
    }
    /* normal subopt */
    else if(!zuker){
      /* first lines of output (suitable  for sort +1n) */
      if (fname[0] != '\0')
        printf("> %s [%d]\n", fname, delta);

      vrna_subopt(vc, delta, subopt_sorted, stdout);

      if (dos) {
        int i;
        for (i=0; i<= MAXDOS && i<=delta/10; i++) {
          printf("%4d %6d\n", i, density_of_states[i]);
        }
      }
    }
    /* Zuker suboptimals */
    else{
      vrna_subopt_solution_t *zr;

      if (vc->cutpoint != -1) {
        vrna_message_error("Sorry, zuker subopts not yet implemented for cofold\n");
      }

      int i;
      if (fname[0] != '\0') printf(">%s\n%s\n", fname, rec_sequence);

      zr = vrna_subopt_zuker(vc);

      putoutzuker(zr);
      (void)fflush(stdout);
      for (i=0; zr[i].structure; i++) {
        free(zr[i].structure);
      }
      free(zr);
    }
    
    (void)fflush(stdout);

    /* clean up */
    vrna_fold_compound_free(vc);

    if(cstruc) free(cstruc);
    if(rec_id) free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    free(structure);

    /* free the rest of current dataset */
    if(rec_rest){
      for(i=0;rec_rest[i];i++) free(rec_rest[i]);
      free(rec_rest);
    }
    rec_id = rec_sequence = orig_sequence = structure = cstruc = NULL;
    rec_rest = NULL;

    /* print user help for the next round if we get input from tty */
    if(istty){
      if(!zuker)
        printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      if(fold_constrained){
        vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK);
        vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n");
      }
      else vrna_message_input_seq_simple();
    }
  }
  return EXIT_SUCCESS;
}

PRIVATE void putoutzuker(vrna_subopt_solution_t* zukersolution) {
  int i;
  printf("%s [%.2f]\n",zukersolution[0].structure,zukersolution[0].energy/100.);
  for(i=1; zukersolution[i].structure; i++) {
    printf("%s [%.2f]\n", zukersolution[i].structure, zukersolution[i].energy/100.);
  }
  return;
}
