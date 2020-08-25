/*
 *                Ineractive Access to folding Routines
 *
 *                c Ivo L Hofacker
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/commands.h"
#include "ViennaRNA/constraints/SHAPE.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/loops/external.h"

#include "RNALfold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#include "ViennaRNA/color_output.inc"

typedef struct {
  FILE  *output;
  int   dangle_model;
} hit_data;


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data);


#endif

PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data);

struct local_struct {
  short               *pt;
  unsigned int        start;
  unsigned int        shift;
  unsigned int        length;
  int                 energy;
  unsigned char       valid;
  struct local_struct *next_entry;
};


int
main(int  argc,
     char *argv[])
{
  FILE                        *input, *output;
  struct  RNALfold_args_info  args_info;
  char                        *ParamFile, *ns_bases, *rec_sequence, *rec_id, **rec_rest,
                              *command_file, *orig_sequence, *infile, *outfile, *filename_delim,
                              *shape_file, *shape_method, *shape_conversion;
  unsigned int                rec_type, read_opt;
  int                         length, istty, noconv, maxdist, zsc, tofile, filename_full,
                              with_shapes, verbose, backtrack;
  double                      min_en, min_z;
  long                        file_pos;
  vrna_md_t                   md;
  vrna_cmd_t                  commands;
  dataset_id                  id_control;

  ParamFile     = ns_bases = NULL;
  do_backtrack  = 1;
  noconv        = 0;
  backtrack     = 0;
  dangles       = 2;
  maxdist       = 150;
  zsc           = 0;
  min_z         = -2.0;
  gquad         = 0;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  outfile       = NULL;
  infile        = NULL;
  input         = NULL;
  output        = NULL;
  tofile        = 0;
  filename_full = 0;
  command_file  = NULL;
  commands      = NULL;
  file_pos      = -1;

  /* apply default model details */
  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNALfold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* parse options for ID manipulation */
  ggo_get_id_control(args_info, id_control, "Sequence", "sequence", "_", 4, 1);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      vrna_message_warning(
        "required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }

  /* do not allow weak pairs */
  if (args_info.backtrack_given)
    backtrack = tofile = 1;

  /* do not allow weak pairs */
  if (args_info.noLP_given)
    md.noLP = noLonelyPairs = 1;

  /* do not allow wobble pairs (GU) */
  if (args_info.noGU_given)
    md.noGU = noGU = 1;

  /* do not allow weak closing pairs (AU,GU) */
  if (args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    noconv = 1;

  /* set energy model */
  if (args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;

  /* take another energy parameter set */
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* set the maximum base pair span */
  if (args_info.span_given)
    maxdist = args_info.span_arg;

  if (args_info.zscore_given) {
#ifdef VRNA_WITH_SVM
    zsc = 1;
    if (args_info.zscore_arg != -2)
      min_z = args_info.zscore_arg;

#else
    vrna_message_error("\'z\' option is available only if compiled with SVM support!");
#endif
  }

  /* gquadruplex support */
  if (args_info.gquad_given)
    md.gquad = gquad = 1;

  if (args_info.verbose_given)
    verbose = 1;

  /* SHAPE reactivity data */
  ggo_get_SHAPE(args_info, with_shapes, shape_file, shape_method, shape_conversion);

  if (args_info.outfile_given) {
    tofile = 1;
    if (args_info.outfile_arg)
      outfile = strdup(args_info.outfile_arg);
  }

  if (args_info.infile_given)
    infile = strdup(args_info.infile_arg);

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(id_control))
    filename_delim = strdup(get_id_delim(id_control));
  else
    filename_delim = NULL;

  if ((filename_delim) && isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

  /* full filename from FASTA header support */
  if (args_info.filename_full_given)
    filename_full = 1;

  if (args_info.commands_given)
    command_file = strdup(args_info.commands_arg);

  /* check for errorneous parameter options */
  if (maxdist <= 0) {
    RNALfold_cmdline_parser_print_help();
    exit(EXIT_FAILURE);
  }

  /* free allocated memory of command line data structure */
  RNALfold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */

  md.max_bp_span = md.window_size = maxdist;

  if (infile) {
    input = fopen((const char *)infile, "r");
    if (!input)
      vrna_message_error("Could not read input file");
  } else {
    input = stdin;
  }

  if (ParamFile != NULL) {
    if (!strcmp(ParamFile, "DNA"))
        vrna_params_load_DNA_Mathews2004();
    else
      vrna_params_load(ParamFile, VRNA_PARAMETER_FORMAT_DEFAULT);
  }

  if (command_file != NULL)
    commands = vrna_file_commands_read(command_file, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

  if (ns_bases != NULL)
    vrna_md_set_nonstandards(&md, ns_bases);

  istty     = (!infile) && isatty(fileno(stdout)) && isatty(fileno(stdin));
  read_opt  |= VRNA_INPUT_NO_REST;
  if (istty) {
    vrna_message_input_seq_simple();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  while (
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, input, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    char  *SEQ_ID       = NULL;
    char  *v_file_name  = NULL;
    char  *tmp_string   = NULL;
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    if (rec_id) /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));

    /* construct the sequence ID */
    set_next_id(&rec_id, id_control);
    SEQ_ID = fileprefix_from_id(rec_id, id_control, filename_full);

    if (tofile) {
      /* prepare the file name */
      if (outfile)
        v_file_name = vrna_strdup_printf("%s", outfile);
      else
        v_file_name = (SEQ_ID) ?
                      vrna_strdup_printf("%s.lfold", SEQ_ID) :
                      vrna_strdup_printf("RNALfold_output.lfold");

      tmp_string = vrna_filename_sanitize(v_file_name, filename_delim);
      free(v_file_name);
      v_file_name = tmp_string;

      if (infile && !strcmp(infile, v_file_name))
        vrna_message_error("Input and output file names are identical");

      output = fopen((const char *)v_file_name, "a");

      if (!output)
        vrna_message_error("Failed to open file for writing");

      file_pos = ftell(output);
    } else {
      output = stdout;
    }

    if (!istty)
      print_fasta_header(output, rec_id);

    length = (int)strlen(rec_sequence);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(rec_sequence);

    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    if (!tofile && istty)
      vrna_message_info(output, "length = %d", length);

    /*
     ########################################################
     # done with 'stdin' handling
     # begin actual computations
     ########################################################
     */

    vrna_fold_compound_t *vc = vrna_fold_compound((const char *)rec_sequence,
                                                  &md,
                                                  VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

    if (commands)
      vrna_commands_apply(vc, commands, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

    if (with_shapes) {
      vrna_constraints_add_SHAPE(vc,
                                 shape_file,
                                 shape_method,
                                 shape_conversion,
                                 verbose,
                                 VRNA_OPTION_WINDOW);
    }

    hit_data data;
    data.output       = output;
    data.dangle_model = md.dangles;

#ifdef VRNA_WITH_SVM
    min_en =
      (zsc) ? vrna_mfe_window_zscore_cb(vc, min_z, &default_callback_z,
                                        (void *)&data) : vrna_mfe_window_cb(vc, &default_callback,
                                                                            (void *)&data);
#else
    min_en = vrna_mfe_window_cb(vc, &default_callback, (void *)&data);
#endif
    fprintf(output, "%s\n", orig_sequence);

    char *msg = NULL;
    if (!tofile && istty)
      msg = vrna_strdup_printf(" minimum free energy = %6.2f kcal/mol", min_en);
    else
      msg = vrna_strdup_printf(" (%6.2f)", min_en);

    print_structure(output, NULL, msg);
    free(msg);

    if (output)
      (void)fflush(output);

    if (tofile && output) {
      fclose(output);
      output = NULL;
    }

    if (backtrack) {
      FILE *f = fopen((const char *)v_file_name, "r");
      if (f) {
        if (fseek(f, file_pos, SEEK_SET) != -1) {
          size_t num_lines, mem_lines;
          num_lines = 0;
          mem_lines = 1024;

          long *lines = (long *)vrna_alloc(sizeof(long) * mem_lines);

          lines[num_lines++] = ftell(f);

          do {
            /* increase memory if necessary */
            if (num_lines == mem_lines) {
              mem_lines *= 1.4;
              lines = (long *)vrna_realloc(lines, sizeof(long) * mem_lines);
            }

            /* seek to next newline char */
            do {
              char c = fgetc(f);
              if ((feof(f)) || (c == '\n'))
                break;
            } while (1);
            
            /* stop at end of file */
            if (feof(f))
              break;

            lines[num_lines++] = ftell(f);
          } while(1);

          if (num_lines > 0) {
            num_lines--;
            char                *mfe_structure = (char *)vrna_alloc(sizeof(char) * (vc->length + 1));
            struct local_struct *ss = (struct local_struct *)vrna_alloc(sizeof(struct local_struct) * (maxdist + 1));
            size_t              num_ss = 0;

            for (size_t i = num_lines - 1; ; i--) {

              printf("line %u at %ld\n", i, lines[i]);
              fseek(f, lines[i], SEEK_SET);
              char *l = vrna_read_line(f);
              long int  start = 0;
              float     en = INF / 100.;
              char      *structure = (char *)vrna_alloc(sizeof(char) * (strlen(l) + 1));
              if (sscanf(l, "%[.()] %*c %f %*c %ld", structure, &en, &start) == 3) {
                printf("s: %s, en: %6.2f, start: %ld\n", structure, en, start);
                ss[num_ss].pt     = vrna_ptable(structure);
                ss[num_ss].start  = start;
                ss[num_ss].shift  = 0;
                ss[num_ss].length = strlen(structure);
                ss[num_ss].energy = (int)(en * 100.);
                ss[num_ss].valid  = 1;
                num_ss++;
              } else {
                printf("%s\n", l);
              }
              free(structure);
              free(l);

              if (i == 0)
                break;
            }

            /* start backtracing */
            memset(mfe_structure, (int)'.', vc->length);

            int *f3 = vc->matrices->f3_local;

            /* The last structure is always part of the full length MFE */
            size_t s = 0;
            for (size_t l = 1; l <= ss[s].length; l++)
              if (ss[s].pt[l] > l) {
                mfe_structure[ss[s].start + l - 1] = '(';
                mfe_structure[ss[s].start + ss[s].pt[l] - 1] = ')';
              }

            /* truncate other structures overlapping with the last one */
            size_t min_s = 1;
            for (unsigned int l = ss[s].start + 1; l < ss[s].start + ss[s].length ; l++) {
              printf("l=%u\n", l);
              for (size_t sss = min_s; sss < num_ss; sss++) {
                /* stop correction if structure starts later */
                if (ss[sss].start > l)
                  break;

                /* only correct structure if it is not entirely subsumed */
                if (ss[sss].valid) {
                  unsigned int i, j, i_local, j_local;
                  i       = ss[sss].start;
                  j       = i + ss[sss].length - 1;
                  i_local = l - i + 1;
                  j_local = ss[sss].pt[i_local];
                  char *seq_local = (char *)vrna_alloc(sizeof(char) * (j - i + 2));
                  memcpy(seq_local, vc->sequence + i - 1, (sizeof(char) * (j - i + 1)));

                  printf("truncating and updating structure %u, interval [%d:%d]\n%s\n%s (%6.2f) shift:%d\n", sss, i, j, seq_local, vrna_db_from_ptable(ss[sss].pt), (float)ss[sss].energy / 100., ss[sss].shift);
                  if (j <= l) {
                    ss[sss].valid = 0;
                  } else if (j_local != 0) {
                    /* we will remove the base pair (i_local, j_local) */
                    /* compute energy difference of the removal step */
                    /* This requires further work to acknowledge dangle model!!!! */
                    int ediff = vrna_eval_move_pt_simple(seq_local, ss[sss].pt, -i_local, -j_local);
                    /* remove the pair */
                    ss[sss].pt[i_local] = ss[sss].pt[j_local] = 0;
                    /* update energy for the remaining substructure */
                    ss[sss].energy += ediff;
                  } else if ((j > l) &&
                             (ss[sss].pt[i_local + 1] != 0)) {
                    /* i is unpaired, we will remove it's 5' dangle contribution to base pair (i + 1, p) */
                    unsigned int  type;
                    int           ediff, d5, d3;
                    
                    j_local = ss[sss].pt[i_local + 1];
                    
                    switch (dangles) {
                      case 2:
                        d5    = vc->sequence_encoding2[i + i_local - 1];
                        d3    = (j_local + 1 <= i_local + ss[sss].length - 1) ? vc->sequence_encoding2[i + j_local + 1 - 1] : -1;
                        type  = vrna_get_ptype_md(vc->sequence_encoding2[i + i_local - 1 + 1],
                                                  vc->sequence_encoding2[i + j_local - 1],
                                                  &(vc->params->model_details));
                        printf("pair (%d, %d) energy: %d vs. %d (%d, %d, type: %d)\n", i_local + 1, j_local, vrna_E_ext_stem(type, -1, d3, vc->params), vrna_E_ext_stem(type, d5, d3, vc->params), d5, d3, type);
                        ediff = vrna_E_ext_stem(type, -1, d3, vc->params) -
                                vrna_E_ext_stem(type, d5, d3, vc->params);
                        break;

                      case 0:
                        ediff = 0;
                        break;

                      default:
                        ediff = 0;
                        break;
                    }
                    /* update energy for the remaining substructure */
                    ss[sss].energy += ediff;
                  } else if (j > l) {
                    /* nothing to do here, since we simply remove an unpaired nucleotide that contributes 0 */
                  }
                  ss[sss].shift++;
                  printf("to\n%s (%6.2f) shift:%d\n", vrna_db_from_ptable(ss[sss].pt), (float)ss[sss].energy / 100., ss[sss].shift);
                }
              }
            }
            
            printf("%s\n", mfe_structure);
            for (unsigned int i = 1; i <= vc->length; i++)
              printf("f3[%d] = %d\n", i, f3[i]);
          }
        }
        fclose(f);

      }
    }

    /* clean up */
    vrna_fold_compound_free(vc);
    free(rec_id);
    free(SEQ_ID);
    free(rec_sequence);
    free(orig_sequence);
    rec_id    = rec_sequence = orig_sequence = NULL;
    rec_rest  = NULL;

    free(v_file_name);

    if (with_shapes)
      break;

    /* print user help for the next round if we get input from tty */

    if (istty)
      vrna_message_input_seq_simple();
  }

  if (infile && input)
    fclose(input);

  free(filename_delim);
  free(command_file);
  vrna_commands_free(commands);

  free_id_data(id_control);

  return EXIT_SUCCESS;
}


PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data)
{
  FILE  *output       = ((hit_data *)data)->output;
  int   dangle_model  = ((hit_data *)data)->dangle_model;
  char  *struct_d2    = NULL;
  char  *msg          = NULL;

  if ((dangle_model == 2) && (start > 1)) {
    msg       = vrna_strdup_printf(" (%6.2f) %4d", en, start - 1);
    struct_d2 = vrna_strdup_printf(".%s", structure);
    print_structure(output, struct_d2, msg);
    free(struct_d2);
  } else {
    msg = vrna_strdup_printf(" (%6.2f) %4d", en, start);
    print_structure(output, structure, msg);
  }

  free(msg);
}


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data)
{
  FILE  *output       = ((hit_data *)data)->output;
  int   dangle_model  = ((hit_data *)data)->dangle_model;
  char  *struct_d2    = NULL;
  char  *msg          = NULL;

  if ((dangle_model == 2) && (start > 1)) {
    msg       = vrna_strdup_printf(" (%6.2f) %4d z= %.3f", en, start - 1, zscore);
    struct_d2 = vrna_strdup_printf(".%s", structure);
    print_structure(output, struct_d2, msg);
    free(struct_d2);
  } else {
    msg = vrna_strdup_printf(" (%6.2f) %4d z= %.3f", en, start, zscore);
    print_structure(output, structure, msg);
  }

  free(msg);
}


#endif
