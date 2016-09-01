/*
                  Ineractive Access to folding Routines

                  c Ivo L Hofacker
                  Vienna RNA package
*/

/** \file
 *  \brief RNAfold program source code
 *
 *  This code provides an interface for MFE and Partition function folding
 *  of single linear or circular RNA molecules.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
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
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/centroid.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/ligand.h"
#include "ViennaRNA/file_formats.h"
#include "RNAfold_cmdl.h"

#include "ViennaRNA/color_output.inc"

/*--------------------------------------------------------------------------*/

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
    if(method != 'W'){
      char *msg = NULL;
      if(method == 'Z')
        asprintf( &msg,
                  "Using SHAPE method '%c' with parameter p1=%f",
                  method, p1);
      else
        asprintf( &msg,
                  "Using SHAPE method '%c' with parameters p1=%f and p2=%f",
                  method, p1, p2);
      vrna_message_info(stderr, msg);
      free(msg);
    }
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

static void
add_ligand_motif( vrna_fold_compound_t *vc,
                  char *motifstring,
                  int verbose,
                  unsigned int options){

  int r, l, n, i, error;
  char *seq, *str, *ptr;
  float energy;

  l = strlen(motifstring);
  seq = vrna_alloc(sizeof(char) * (l + 1));
  str = vrna_alloc(sizeof(char) * (l + 1));

  error = 1;

  if(motifstring){
    error = 0;
    /* parse sequence */
    for(r = 0, ptr = motifstring; *ptr != '\0'; ptr++){
      if(*ptr == ',')
        break;
      seq[r++] = *ptr;
      toupper(seq[r-1]);
    }
    seq[r] = '\0';
    seq = vrna_realloc(seq, sizeof(char) * (strlen(seq) + 1));

    for(ptr++, r = 0; *ptr != '\0'; ptr++){
      if(*ptr == ',')
        break;
      str[r++] = *ptr;
    }
    str[r] = '\0';
    str = vrna_realloc(str, sizeof(char) * (strlen(seq) + 1));

    ptr++;
    if(!(sscanf(ptr, "%f", &energy) == 1)){
      vrna_message_warning("Energy contribution in ligand motif missing!");
      error = 1;
    }
    if(strlen(seq) != strlen(str)){
      vrna_message_warning("Sequence and structure length in ligand motif have unequal lengths!");
      error = 1;
    }
    if(strlen(seq) == 0){
      vrna_message_warning("Sequence length in ligand motif is zero!");
      error = 1;
    }

    if(!error && verbose){
      char *msg = NULL;
      asprintf( &msg,
                "Read ligand motif: %s, %s, %f\n",
                seq, str, energy);
      vrna_message_info(stderr, msg);
      free(msg);
    }
  }

  if(error || (!vrna_sc_add_hi_motif(vc, seq, str, energy, options))){
    vrna_message_warning("Malformatted ligand motif! Skipping stabilizing motif.");
  }

  free(seq);
  free(str);
}

int main(int argc, char *argv[]){
  FILE            *input, *output;
  struct          RNAfold_args_info args_info;
  char            *buf, *rec_sequence, *rec_id, **rec_rest, *structure, *cstruc, *orig_sequence;
  char            *constraints_file, *shape_file, *shape_method, *shape_conversion;
  char            fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH], *ParamFile;
  char            *ns_bases, *c;
  int             i, length, l, cl, sym, istty, pf, noPS, noconv, do_bpp, enforceConstraints,
                  batch, auto_id, id_digits;
  long int        seq_number;
  unsigned int    rec_type, read_opt;
  double          energy, min_en, kT, sfact;
  int             doMEA, circular, lucky, with_shapes, verbose, istty_in, istty_out;
  double          MEAgamma, bppmThreshold, betaScale;
  char            *infile, *outfile, *ligandMotif, *id_prefix;

  vrna_md_t         md;

  rec_type      = read_opt = 0;
  rec_id        = buf = rec_sequence = structure = cstruc = orig_sequence = NULL;
  rec_rest      = NULL;
  ParamFile     = NULL;
  ns_bases      = NULL;
  do_bpp        = do_backtrack  = 1;  /* set local (do_bpp) and global (do_backtrack) default for bpp computation */
  pf            = 0;
  sfact         = 1.07;
  noPS          = 0;
  noconv        = 0;
  circular      = 0;
  gquad         = 0;
  cl            = l = length = 0;
  dangles       = 2;
  MEAgamma      = 1.;
  bppmThreshold = 1e-5;
  lucky         = 0;
  doMEA         = 0;
  betaScale     = 1.;
  shape_file    = NULL;
  shape_method  = NULL;
  with_shapes   = 0;
  verbose       = 0;
  max_bp_span   = -1;
  constraints_file = NULL;
  enforceConstraints  = 0;
  batch               = 0;
  seq_number          = 1;
  id_prefix           = NULL;
  auto_id             = 0;
  id_digits           = 4;
  outfile       = NULL;
  infile        = NULL;
  input         = NULL;
  output        = NULL;
  ligandMotif   = NULL;

  /* apply default model details */
  set_model_details(&md);


  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given){
    fold_constrained=1;
    if(args_info.constraint_arg[0] != '\0')
      constraints_file = strdup(args_info.constraint_arg);
  }
  /* enforce base pairs given in constraint string rather than weak enforce */
  if(args_info.enforceConstraint_given)
    enforceConstraints = 1;

  /* do batch jobs despite constraints read from input file */
  if(args_info.batch_given)
    batch = 1;

  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     md.special_hp = tetra_loop=0;
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
  if(args_info.canonicalBPonly_given)       md.canonicalBPonly = canonicalBPonly = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)      noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given) md.energy_set = energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)     md.sfact = sfact = args_info.pfScale_arg;
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)        md.circ = circular = 1;
  /* always look on the bright side of life */
  if(args_info.ImFeelingLucky_given)  md.uniq_ML = lucky = pf = st_back = 1;
  /* set the bppm threshold for the dotplot */
  if(args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0.,args_info.bppmThreshold_arg));
  if(args_info.betaScale_given)   md.betaScale = betaScale = args_info.betaScale_arg;
  /* do not produce postscript output */
  if(args_info.noPS_given)        noPS=1;
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != 1)
      do_bpp = md.compute_bpp = do_backtrack = args_info.partfunc_arg;
  }
  /* MEA (maximum expected accuracy) settings */
  if(args_info.MEA_given){
    pf = doMEA = 1;
    if(args_info.MEA_arg != -1)
      MEAgamma = args_info.MEA_arg;
  }
  if(args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;
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
  if(args_info.outfile_given){
    outfile = strdup(args_info.outfile_arg);
  }
  if(args_info.infile_given){
    infile = strdup(args_info.infile_arg);
  }
  if(args_info.motif_given){
    ligandMotif = strdup(args_info.motif_arg);
  }

  if(args_info.auto_id_given){
    auto_id = 1;
  }

  if(args_info.id_prefix_given){
    id_prefix   = strdup(args_info.id_prefix_arg);
    auto_id  = 1;
  } else {
    id_prefix = strdup("sequence");
  }

  /* set width of alignment number in the output */
  if(args_info.id_digits_given){
    if((args_info.id_digits_arg > 0) && (args_info.id_digits_arg < 19))
      id_digits = args_info.id_digits_arg;
    else
      vrna_message_warning("ID number digits out of allowed range! Using defaults...");
  }

  /* set first alignment number in the output */
  if(args_info.id_start_given){
    if((args_info.id_start_arg >= 0) && (args_info.id_start_arg <= LONG_MAX)){
      seq_number  = args_info.id_start_arg;
      auto_id  = 1;
    } else
      vrna_message_warning("ID number start out of allowed range! Using defaults...");
  }


  /* free allocated memory of command line data structure */
  RNAfold_cmdline_parser_free (&args_info);


  /*
  #############################################
  # begin initializing
  #############################################
  */
  if(infile){
    input = fopen((const char *)infile, "r");
    if(!input)
      vrna_message_error("Could not read input file");
  }

  if(circular && gquad){
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (circular && noLonelyPairs)
    vrna_message_warning("depending on the origin of the circular sequence, some structures may be missed when using -noLP\nTry rotating your sequence a few times");

  if (ns_bases != NULL) {
    vrna_md_set_nonstandards(&md, ns_bases);
  }

  istty_in  = isatty(fileno(stdin));
  istty_out = isatty(fileno(stdout)) && (!outfile);
  istty = (!infile) && isatty(fileno(stdout)) && isatty(fileno(stdin));

  /* print user help if we get input from tty */
  if(istty){
    if(fold_constrained){
      vrna_message_constraint_options_all();
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint");
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
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, input, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    char *prefix      = NULL;
    char *v_file_name = NULL;
    char *SEQ_ID      = NULL;

    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    if(rec_id){
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    /* construct the sequence ID */
    if(outfile){
      if((fname[0] != '\0') && (!auto_id)){ /* append ID as read from file to the user prefix */
        asprintf( &SEQ_ID,
                  "%s_%s",
                  outfile,
                  fname);
      } else if(auto_id){
        asprintf(&SEQ_ID, "%s_%s_%0*ld", outfile, id_prefix, id_digits, seq_number);
      } else {
        asprintf(&SEQ_ID, "%s_%0*ld", outfile, id_digits, seq_number);
      }
    } else if((fname[0] != '\0') && (!auto_id)){ /* we've read an ID from file, so we use it */
      SEQ_ID = strdup(fname);
    } else if(auto_id){ /* we have nuffin', Jon Snow (...so we simply generate an ID) */
      asprintf(&SEQ_ID, "%s_%0*ld", id_prefix, id_digits, seq_number);
    }

    if(outfile){
      /* prepare the file prefix */
      if(fname[0] != '\0'){
        prefix = (char *)vrna_alloc(sizeof(char) * (strlen(fname) + strlen(outfile) + 1));
        strcpy(prefix, outfile);
        strcat(prefix, "_");
        strcat(prefix, fname);
      } else {
        prefix = (char *)vrna_alloc(sizeof(char) * (strlen(outfile) + 1));
        strcpy(prefix, outfile);
      }
    }

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) vrna_seq_toRNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound_t *vc = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));

    length    = vc->length;

    structure = (char *) vrna_alloc(sizeof(char) * (length+1));

    /* parse the rest of the current dataset to obtain a structure constraint */
    if(fold_constrained){
      if(constraints_file){
        /** [Adding hard constraints from file] */
        vrna_constraints_add(vc, constraints_file, VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));
        /** [Adding hard constraints from file] */
      } else {
        cstruc = NULL;
        unsigned int coptions = (rec_id) ? VRNA_OPTION_MULTILINE : 0;
        cstruc = vrna_extract_record_rest_structure((const char **)rec_rest, 0, coptions);
        cl = (cstruc) ? (int)strlen(cstruc) : 0;

        if(cl == 0)           vrna_message_warning("structure constraint is missing");
        else if(cl < length)  vrna_message_warning("structure constraint is shorter than sequence");
        else if(cl > length)  vrna_message_error("structure constraint is too long");
        if(cstruc){
          strncpy(structure, cstruc, sizeof(char)*(cl+1));

          /** [Adding hard constraints from pseudo dot-bracket] */
          unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;
          if(enforceConstraints)
            constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;
          vrna_constraints_add(vc, (const char *)structure, constraint_options);
          /** [Adding hard constraints from pseudo dot-bracket] */
        }
      }
    }

    if(with_shapes)
      add_shape_constraints(vc, shape_method, shape_conversion, shape_file, verbose, VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));

    if(ligandMotif)
      add_ligand_motif(vc, ligandMotif, verbose, VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));

    if(istty){
      char *msg = NULL;
      asprintf( &msg,
                "length = %d\n",
                length);
      vrna_message_info(stdout, msg);
      free(msg);
    }

    if(outfile){
      v_file_name = (char *)vrna_alloc(sizeof(char) * (strlen(prefix) + 8));
      strcpy(v_file_name, prefix);
      strcat(v_file_name, ".fold");

      if(infile && !strcmp(infile, v_file_name))
        vrna_message_error("Input and output file names are identical");

      output = fopen((const char *)v_file_name, "a");
      if(!output)
        vrna_message_error("Failed to open file for writing");
    } else {
      output = stdout;
    }

    /*
    ########################################################
    # begin actual computations
    ########################################################
    */
    min_en = (double)vrna_mfe(vc, structure);

    /* check whether the constraint allows for any solution */
    if(fold_constrained && constraints_file){
      if(min_en == (double)(INF/100.)){
        char *msg = NULL;
        asprintf( &msg,
                  "Supplied structure constraints create empty solution set for sequence:\n%s",
                  orig_sequence);
        vrna_message_error(msg);
        free(msg);
        exit(EXIT_FAILURE);
      }
    }

    if(output){
      print_fasta_header(output, SEQ_ID);
      fprintf(output, "%s\n", orig_sequence);
    }

    if(!lucky){
      if(output){
        char *msg = NULL;
        if(istty)
          asprintf( &msg,
                    "\n minimum free energy = %6.2f kcal/mol",
                    min_en);
        else
          asprintf( &msg,
                    " (%6.2f)",
                    min_en);
        print_structure(output, structure, msg);
        free(msg);
        (void) fflush(output);

      }

      if(fname[0] != '\0'){
        strcpy(ffname, fname);
        strcat(ffname, "_ss.ps");
      } else strcpy(ffname, "rna.ps");

      if(!noPS){
        char *filename_plot = NULL;
        if(SEQ_ID)
          asprintf( &filename_plot,
                    "%s_ss.ps",
                    SEQ_ID);
        else
          filename_plot = strdup("rna.ps");

        if(ligandMotif){
          int a,b,c,d, cnt;
          char *annote;
          cnt     = 1;
          a       = 1;
          annote  = vrna_alloc(sizeof(char) * cnt * 64);
          while(vrna_sc_detect_hi_motif(vc, structure, &a, &b, &c, &d)){
            char segment[64];
            if(c != 0){
              if(verbose){
                char *msg = NULL;
                asprintf( &msg,
                          "specified motif detected in MFE structure: (%d,%d) (%d,%d)",
                          a, b, c, d);
                vrna_message_info(stdout, msg);
                free(msg);
              }
              sprintf(segment, " %d %d %d %d 1. 0 0 BFmark", a, b, c, d);
              a = c;
            }
            else{
              if(verbose){
                char *msg = NULL;
                asprintf( &msg,
                          "specified motif detected in MFE structure: (%d,%d)",
                          a, b);
                vrna_message_info(stdout, msg);
                free(msg);
              }
              sprintf(segment, " %d %d 1. 0 0 Fomark", a, b);
              a = b;
            }
            strcat(annote, segment);
            annote = vrna_realloc(annote, sizeof(char) * ++cnt * 64);
          }
          (void) vrna_file_PS_rnaplot_a(orig_sequence, structure, filename_plot, annote, NULL, &md);
          free(annote);
        } else {
          (void) vrna_file_PS_rnaplot_a(orig_sequence, structure, filename_plot, NULL, NULL, &md);
        }
        free(filename_plot);
      }
    }

    if (length>2000)
      vrna_mx_mfe_free(vc);

    if (pf) {
      char *pf_struc = (char *) vrna_alloc((unsigned) length+1);
      if (vc->params->model_details.dangles==1) {
          vc->params->model_details.dangles=2;   /* recompute with dangles as in pf_fold() */
          min_en = vrna_eval_structure(vc, structure);
          vc->params->model_details.dangles=1;
      }

      vrna_exp_params_rescale(vc, &min_en);

      kT = vc->exp_params->kT/1000.;

      if (length>2000){
        char *msg = NULL;
        asprintf(&msg, "scaling factor %f", vc->exp_params->pf_scale);
        vrna_message_info(stderr, msg);
        free(msg);
      }

      fflush(stdout);

      if (cstruc!=NULL) strncpy(pf_struc, cstruc, length+1);

      energy = (double)vrna_pf(vc, pf_struc);

      /* in case we abort because of floating point errors */
      if (length>1600){
        char *msg = NULL;
        asprintf( &msg,
                  "free energy = %8.2f",
                  energy);
        vrna_message_info(stderr, msg);
        free(msg);
      }
      if(lucky){
        vrna_init_rand();
        char *filename_plot = NULL;
        char *s = vrna_pbacktrack(vc);
        min_en = vrna_eval_structure(vc, (const char *)s);
        if(output){
          char *energy_string = NULL;
          if(istty_in)
            asprintf( &energy_string,
                      "\n free energy = %6.2f kcal/mol",
                      min_en);
          else
            asprintf(&energy_string, " (%6.2f)", min_en);
          print_structure(output, s, energy_string);
          free(energy_string);
          (void) fflush(output);
        }

        if(SEQ_ID)
          asprintf( &filename_plot,
                    "%s_ss.ps",
                    SEQ_ID);
        else
          filename_plot = strdup("rna.ps");

        if (!noPS)
          (void) vrna_file_PS_rnaplot(orig_sequence, s, filename_plot, &md);
        free(s);
        free(filename_plot);
      }
      else{
        if (do_bpp) {
          if(output){
            char *msg = NULL;
            if(istty_in)
              asprintf( &msg,
                        "\n free energy of ensemble = %6.2f kcal/mol",
                        energy);
            else
              asprintf( &msg, " [%6.2f]", energy);

            print_structure(output, pf_struc, msg);
            free(msg);
          }
        } else {
          char *msg = NULL;
          asprintf(&msg, " free energy of ensemble = %6.2f kcal/mol", energy);
          print_structure(output, NULL, msg);
          free(msg);
        }

        if (do_bpp) {
          plist *pl1,*pl2;
          char *cent;
          double dist, cent_en;

          pl1     = vrna_plist_from_probs(vc, bppmThreshold);
          pl2     = vrna_plist(structure, 0.95*0.95);

          if(ligandMotif){
            /* append motif positions to the plists of base pair probabilities */
            vrna_plist_t *motifs, *ptr;
            int a,b,c,d, cnt, size, add;
            cnt = 0;
            a   = 1;
            add = 10;
            /* get size of pl1 */
            for(size = 0, ptr = pl1; ptr->i; size++, ptr++);

            /* increase length of pl1 */
            pl1 = vrna_realloc(pl1, sizeof(vrna_plist_t) * (size + add + 1));

            while(vrna_sc_get_hi_motif(vc, &a, &b, &c, &d)){
              if(c == 0){ /* hairpin motif */
                pl1[size + cnt].i = a;
                pl1[size + cnt].j = b;
                pl1[size + cnt].p = 0.95*0.95;
                pl1[size + cnt].type = 2;
                cnt++;
                if(cnt == add){
                  add += 10;
                  /* increase length of pl1 */
                  pl1 = vrna_realloc(pl1, sizeof(vrna_plist_t) * (size + add + 1));
                }
              } else { /* interior loop motif */
                pl1[size + cnt].i = a;
                pl1[size + cnt].j = b;
                pl1[size + cnt].p = 0.95*0.95;
                pl1[size + cnt].type = 3;
                cnt++;
                pl1[size + cnt].i = c;
                pl1[size + cnt].j = d;
                pl1[size + cnt].p = 0.95*0.95;
                pl1[size + cnt].type = 3;
                cnt++;
                if(cnt == add){
                  add += 10;
                  /* increase length of pl1 */
                  pl1 = vrna_realloc(pl1, sizeof(vrna_plist_t) * (size + add + 1));
                }
              }

              a = b;
            }
            /* resize pl1 to actual needs */
            pl1 = vrna_realloc(pl1, sizeof(vrna_plist_t) * (size + cnt + 1));
            pl1[size + cnt].i = 0;
            pl1[size + cnt].j = 0;

            /* now scan for the motif in MFE structure again */
            add = 10;
            a   = 1;
            cnt = 0;
            /* get size of pl2 */
            for(size = 0, ptr = pl2; ptr->i; size++, ptr++);

            /* increase length of pl2 */
            pl2 = vrna_realloc(pl2, sizeof(vrna_plist_t) * (size + add + 1));
            while(vrna_sc_detect_hi_motif(vc, structure, &a, &b, &c, &d)){
              if(c == 0){ /* hairpin motif */
                pl2[size + cnt].i = a;
                pl2[size + cnt].j = b;
                pl2[size + cnt].p = 0.95*0.95;
                pl2[size + cnt].type = 2;
                cnt++;
                if(cnt == add){
                  add += 10;
                  /* increase length of pl1 */
                  pl2 = vrna_realloc(pl2, sizeof(vrna_plist_t) * (size + add + 1));
                }
              } else { /* interior loop motif */
                pl2[size + cnt].i = a;
                pl2[size + cnt].j = b;
                pl2[size + cnt].p = 0.95*0.95;
                pl2[size + cnt].type = 3;
                cnt++;
                pl2[size + cnt].i = c;
                pl2[size + cnt].j = d;
                pl2[size + cnt].p = 0.95*0.95;
                pl2[size + cnt].type = 3;
                cnt++;
                if(cnt == add){
                  add += 10;
                  /* increase length of pl1 */
                  pl2 = vrna_realloc(pl2, sizeof(vrna_plist_t) * (size + add + 1));
                }
              }
              a = b;
            }
            /* resize pl1 to actual needs */
            pl2 = vrna_realloc(pl2, sizeof(vrna_plist_t) * (size + cnt + 1));
            pl2[size + cnt].i = 0;
            pl2[size + cnt].j = 0;
          }

          char *filename_dotplot = NULL;
          if(SEQ_ID)
            asprintf( &filename_dotplot,
                      "%s_dp.ps",
                      SEQ_ID);
          else
            filename_dotplot = strdup("dot.ps");

          (void) PS_dot_plot_list(orig_sequence, filename_dotplot, pl1, pl2, "");
          free(filename_dotplot);

          cent    = vrna_centroid(vc, &dist);
          cent_en = vrna_eval_structure(vc, (const char *)cent);
          if(output){
            char *msg = NULL;
            asprintf( &msg, " {%6.2f d=%.2f}", cent_en, dist);
            print_structure(output, cent, msg);
            free(msg);
          }

          if(ligandMotif){
            int a,b,c,d;
            a = 1;
            while(vrna_sc_detect_hi_motif(vc, structure, &a, &b, &c, &d)){
              if(c != 0){
                if(verbose){
                  char *msg = NULL;
                  asprintf( &msg,
                            "specified motif detected in centroid structure: (%d,%d) (%d,%d)",
                            a, b, c, d);
                  vrna_message_info(stdout, msg);
                  free(msg);
                }
                a = c;
              }
              else{
                if(verbose){
                  char *msg = NULL;
                  asprintf( &msg,
                            "specified motif detected in centroid structure: (%d,%d)",
                            a, b);
                  vrna_message_info(stdout, msg);
                  free(msg);
                }
                a = b;
              }
            }
          }

          free(cent);

          free(pl2);
          if (do_bpp==2) {
            char *filename_stackplot = NULL;
            if(SEQ_ID)
              asprintf( &filename_stackplot,
                        "%s_dp2.ps",
                        SEQ_ID);
            else
              filename_stackplot = strdup("dot2.ps");

            pl2 = vrna_stack_prob(vc, 1e-5);

            PS_dot_plot_list(orig_sequence, filename_stackplot, pl1, pl2,
                             "Probabilities for stacked pairs (i,j)(i+1,j-1)");
            free(pl2);
            free(filename_stackplot);
          }
          free(pl1);
          free(pf_struc);
          if(doMEA){
            float mea, mea_en;
            /*  this is a hack since vrna_plist_from_probs() always resolves g-quad pairs,
                while MEA_seq() still expects unresolved gquads */
            int gq = vc->exp_params->model_details.gquad;
            vc->exp_params->model_details.gquad = 0;
            plist *pl = vrna_plist_from_probs(vc, 1e-4/(1+MEAgamma));
            vc->exp_params->model_details.gquad = gq;

            if(gquad){
              mea = MEA_seq(pl, rec_sequence, structure, MEAgamma, vc->exp_params);
            } else {
              mea = MEA(pl, structure, MEAgamma);
            }
            mea_en = vrna_eval_structure(vc, (const char *)structure);
            if(output){
              char *msg = NULL;
              asprintf( &msg, " {%6.2f MEA=%.2f}", mea_en, mea);
              print_structure(output, structure, msg);
              free(msg);
            }
            if(ligandMotif){
              int a,b,c,d;
              a = 1;
              while(vrna_sc_detect_hi_motif(vc, structure, &a, &b, &c, &d)){
                if(c != 0){
                  if(verbose){
                    char *msg = NULL;
                    asprintf( &msg,
                              "specified motif detected in MEA structure: (%d,%d) (%d,%d)",
                              a, b, c, d);
                    vrna_message_info(stdout, msg);
                    free(msg);
                  }
                  a = c;
                }
                else{
                  if(verbose){
                    char *msg = NULL;
                    asprintf( &msg,
                              "specified motif detected in MEA structure: (%d,%d)",
                              a, b);
                    vrna_message_info(stdout, msg);
                    free(msg);
                  }
                  a = b;
                }
              }
            }

            free(pl);
          }
        }
        if(output){
          char *msg = NULL;
          if(do_bpp){
            asprintf( &msg,
                      " frequency of mfe structure in ensemble %g"
                      "; ensemble diversity %-6.2f",
                      exp((energy-min_en)/kT),
                      vrna_mean_bp_distance(vc));
          } else {
            asprintf( &msg,
                      " frequency of mfe structure in ensemble %g;",
                      exp((energy-min_en)/kT));
          }
          print_structure(output, NULL, msg);
          free(msg);
        }
      }
    }
    if(output)
      (void) fflush(output);
    if(outfile && output){
      fclose(output);
      output = NULL;
    }

    /* clean up */
    vrna_fold_compound_free(vc);
    free(cstruc);
    free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    free(structure);

    /* free the rest of current dataset */
    if(rec_rest){
      for(i=0;rec_rest[i];i++) free(rec_rest[i]);
      free(rec_rest);
    }
    rec_id = rec_sequence = structure = cstruc = NULL;
    rec_rest = NULL;

    if(with_shapes || (constraints_file && (!batch)))
      break;

    free(SEQ_ID);

    if(seq_number == LONG_MAX){
      vrna_message_warning("Sequence ID number overflow, beginning with 1 (again)!");
      seq_number = 1;
    } else
      seq_number++;

    /* print user help for the next round if we get input from tty */
    if(istty){
      if(fold_constrained){
        vrna_message_constraint_options_all();
        vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint");
      }
      else vrna_message_input_seq_simple();
    }
  }
  
  if(input)
    fclose(input);

  free(constraints_file);
  free(ligandMotif);
  free(shape_method);
  free(shape_conversion);
  free(id_prefix);

  return EXIT_SUCCESS;
}
