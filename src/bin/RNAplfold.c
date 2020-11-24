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
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/SHAPE.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/commands.h"
#include "RNAplfold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#include "ViennaRNA/color_output.inc"

#ifndef isnan
#define isnan(x) \
  (sizeof(x) == sizeof(long double) ? isnan_ld(x) \
   : sizeof(x) == sizeof(double) ? isnan_d(x) \
   : isnan_f(x))

/* Use volatile tmp variable to prevent optimization by compiler. */
static inline int
isnan_f(float x)
{
  volatile float tmp = x;

  return tmp != x;
}


static inline int
isnan_d(double x)
{
  volatile double tmp = x;

  return tmp != x;
}


static inline int
isnan_ld(long double x)
{
  volatile long double tmp = x;

  return tmp != x;
}


#endif /* ifndef isnan */

typedef struct {
  float     cutoff;
  FILE      *pUfp;
  FILE      *spup;
  vrna_ep_t *plist;
  int       plist_cnt;
  int       plexoutput;
  int       simply_putout;
  int       openenergies;
  double    **pup;
  int       ulength;
  int       n;
  double    kT;
} plfold_data;

int unpaired;

PRIVATE void
putoutphakim_u(vrna_fold_compound_t *fc,
               double               **pU,
               int                  length,
               int                  ulength,
               FILE                 *fp);


PRIVATE void
plfold_callback(FLT_OR_DBL    *pr,
                int           pr_size,
                int           i,
                int           max,
                unsigned int  type,
                void          *data);


PRIVATE void
prepare_up_file(plfold_data *data);


PRIVATE void
print_up_open(FILE    *fp,
              int     i,
              double  *pr,
              int     pr_size,
              int     ulength,
              double  kT);


PRIVATE void
print_up(FILE   *fp,
         int    i,
         double *pr,
         int    pr_size,
         int    ulength);


PRIVATE void
print_pu_bin(vrna_fold_compound_t *fc,
             plfold_data          *data,
             int                  ulength);


/*--------------------------------------------------------------------------*/
int
main(int  argc,
     char *argv[])
{
  FILE                        *pUfp;
  struct RNAplfold_args_info  args_info;
  char                        *structure, *ParamFile, *ns_bases, *rec_sequence, *rec_id,
                              **rec_rest, *orig_sequence, *filename_delim, *command_file,
                              *shape_file, *shape_method, *shape_conversion;
  unsigned int                rec_type, read_opt;
  int                         length, istty, winsize, pairdist, tempwin, temppair, tempunpaired,
                              noconv, i, plexoutput, simply_putout, openenergies, binaries,
                              filename_full, with_shapes, verbose;
  float                       cutoff;
  vrna_exp_param_t            *pf_parameters;
  vrna_md_t                   md;
  vrna_cmd_t                  commands;
  dataset_id                  id_control;

  pUfp          = NULL;
  dangles       = 2;
  cutoff        = 0.01;
  winsize       = 70;
  pairdist      = 0;
  unpaired      = 0;
  simply_putout = plexoutput = openenergies = noconv = 0;
  binaries      = 0;
  tempwin       = temppair = tempunpaired = 0;
  structure     = ParamFile = ns_bases = NULL;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  pf_parameters = NULL;
  filename_full = 0;
  command_file  = NULL;
  commands      = NULL;
  verbose       = 0;

  set_model_details(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAplfold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  if (args_info.verbose_given)
    verbose = 1;

  /* SHAPE reactivity data */
  ggo_get_SHAPE(args_info, with_shapes, shape_file, shape_method, shape_conversion);

  /* parse options for ID manipulation */
  ggo_get_id_control(args_info, id_control, "Sequence", "sequence", "_", 4, 1);

  ggo_get_md_part(args_info, md);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
      vrna_message_warning(
        "required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }

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
    pairdist = args_info.span_arg;

  /* set the pair probability cutoff */
  if (args_info.cutoff_given)
    cutoff = args_info.cutoff_arg;

  /* set the windowsize */
  if (args_info.winsize_given)
    winsize = args_info.winsize_arg;

  /* set the length of unstructured region */
  if (args_info.ulength_given)
    unpaired = args_info.ulength_arg;

  /* compute opening energies */
  if (args_info.opening_energies_given)
    openenergies = 1;

  /* print output on the fly */
  if (args_info.print_onthefly_given)
    simply_putout = 1;

  /* turn on RNAplex output */
  if (args_info.plex_output_given)
    plexoutput = 1;

  /* turn on binary output*/
  if (args_info.binaries_given)
    binaries = 1;

  /* check for errorneous parameter options */
  if ((pairdist < 0) || (cutoff < 0.) || (unpaired < 0) || (winsize < 0)) {
    RNAplfold_cmdline_parser_print_help();
    exit(EXIT_FAILURE);
  }

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

  /* free allocated memory of command line data structure */
  RNAplfold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (ParamFile != NULL) {
    if (!strcmp(ParamFile, "DNA"))
        vrna_params_load_DNA_Mathews2004();
    else
      vrna_params_load(ParamFile, VRNA_PARAMETER_FORMAT_DEFAULT);
  }

  if (ns_bases != NULL)
    vrna_md_set_nonstandards(&md, ns_bases);

  if (command_file != NULL)
    commands = vrna_file_commands_read(command_file, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

  /* check parameter options again and reset to reasonable values if needed */
  if (openenergies && !unpaired)
    unpaired = 31;

  if (pairdist == 0)
    pairdist = winsize;

  if (pairdist > winsize) {
    vrna_message_warning("pairdist (-L %d) should be <= winsize (-W %d);"
                         "Setting pairdist=winsize",
                         pairdist, winsize);
    pairdist = winsize;
  }

  if (dangles % 2) {
    vrna_message_warning("using default dangles = 2");
    md.dangles = dangles = 2;
  }

  istty     = isatty(fileno(stdout)) && isatty(fileno(stdin));
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
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    char *SEQ_ID = NULL;
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

    length    = (int)strlen(rec_sequence);
    structure = (char *)vrna_alloc(sizeof(char) * (length + 1));

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(rec_sequence);

    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    if (istty)
      vrna_message_info(stdout, "length = %d", length);

    /*
     ########################################################
     # done with 'stdin' handling
     ########################################################
     */

    if (length > 1000000) {
      if (!simply_putout && !unpaired) {
        vrna_message_warning("Switched to simple output mode!!!");
        simply_putout = 1;
      }
    }

    if ((simply_putout) && (plexoutput)) {
      vrna_message_warning("plexoutput not available in simple output mode!\n"
                           "Switching back to full mode instead!");
      simply_putout = 0;
    }

    if ((simply_putout) && (binaries)) {
      vrna_message_warning("binary output not available in simple output mode!\n"
                           "Switching back to full mode instead!");
      simply_putout = 0;
    }

    /* restore winsize if altered before */
    if (tempwin != 0) {
      winsize = tempwin;
      tempwin = 0;
    }

    /* restore pairdist if altered before */
    if (temppair != 0) {
      pairdist  = temppair;
      temppair  = 0;
    }

    /* restore ulength if altered before */
    if (tempunpaired != 0) {
      unpaired      = tempunpaired;
      tempunpaired  = 0;
    }

    /* adjust winsize, pairdist and ulength if necessary */
    if (length < winsize) {
      vrna_message_warning("window size %d larger than sequence length %d", winsize, length);
      tempwin = winsize;
      winsize = length;
      if (pairdist > winsize) {
        temppair  = pairdist;
        pairdist  = winsize;
      }

      if (unpaired > winsize) {
        tempunpaired  = unpaired;
        unpaired      = winsize;
      }
    }

    /*
     ########################################################
     # begin actual computations
     ########################################################
     */

    if (length > 0) {
      /* construct output file names */
      char *fname1, *fname2, *fname3, *fname4, *ffname, *tmp_string;

      if (!SEQ_ID)
        SEQ_ID = strdup("plfold");

      fname1  = vrna_strdup_printf("%s%slunp", SEQ_ID, filename_delim);
      fname2  = vrna_strdup_printf("%s%sbasepairs", SEQ_ID, filename_delim);
      fname3  = vrna_strdup_printf("%s%suplex", SEQ_ID, filename_delim);
      fname4  = (binaries) ?
                vrna_strdup_printf("%s%sopenen%sbin",
                                   SEQ_ID,
                                   filename_delim,
                                   filename_delim) :
                vrna_strdup_printf("%s%sopenen",
                                   SEQ_ID,
                                   filename_delim);
      ffname = vrna_strdup_printf("%s%sdp.ps", SEQ_ID, filename_delim);

      /* sanitize filenames */
      tmp_string = vrna_filename_sanitize(fname1, filename_delim);
      free(fname1);
      fname1      = tmp_string;
      tmp_string  = vrna_filename_sanitize(fname2, filename_delim);
      free(fname2);
      fname2      = tmp_string;
      tmp_string  = vrna_filename_sanitize(fname3, filename_delim);
      free(fname3);
      fname3      = tmp_string;
      tmp_string  = vrna_filename_sanitize(fname4, filename_delim);
      free(fname4);
      fname4      = tmp_string;
      tmp_string  = vrna_filename_sanitize(ffname, filename_delim);
      free(ffname);
      ffname = tmp_string;


      md.compute_bpp  = 1;
      md.window_size  = winsize;
      md.max_bp_span  = pairdist;

      vrna_fold_compound_t *fc = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_WINDOW);

      if (with_shapes) {
        vrna_constraints_add_SHAPE(fc,
                                   shape_file,
                                   shape_method,
                                   shape_conversion,
                                   verbose,
                                   VRNA_OPTION_DEFAULT | VRNA_OPTION_WINDOW);
      }

      if (commands)
        vrna_commands_apply(fc, commands, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

      pf_parameters = vrna_exp_params(&md);

      /* prepare data structure for callback */
      plfold_data data;

      data.cutoff         = cutoff;
      data.spup           = (simply_putout) ? fopen(fname2, "w") : NULL;
      data.plexoutput     = plexoutput;
      data.simply_putout  = simply_putout;
      data.openenergies   = openenergies;
      data.plist          = NULL;
      data.plist_cnt      = 0;
      data.ulength        = unpaired;
      data.n              = length;
      data.kT             = pf_parameters->kT;

      if (unpaired > 0) {
        if (simply_putout) {
          data.pup  = NULL;
          data.pUfp = fopen(openenergies ? fname4 : fname1, "w");
          prepare_up_file(&data);
        } else {
          /* if we don't print on-the-fly we store unpaired probabilities for later */
          data.pup        = (double **)vrna_alloc(MAX2(unpaired, length + 1) * sizeof(double *));
          data.pup[0]     = (double *)vrna_alloc(sizeof(double));   /*I only need entry 0*/
          data.pup[0][0]  = unpaired;
          data.pUfp       = NULL;
        }
      } else {
        data.pup  = NULL;
        data.pUfp = NULL;
      }

      /* prepare option flags */
      unsigned int plfold_opt = 0;

      /* always compute base pair probabilities */
      plfold_opt |= VRNA_PROBS_WINDOW_BPP;

      if (unpaired > 0)
        plfold_opt |= VRNA_PROBS_WINDOW_UP;

      /* perform recursions */
      int r = vrna_probs_window(fc, unpaired, plfold_opt, &plfold_callback, (void *)&data);

      if (!r) {
        vrna_message_warning("Something bad happened while processing the input! "
                             "Aborting now...");
        goto rnaplfold_exit;
      }

      if (!simply_putout) {
        /* create dot plot output */
        PS_dot_plot_turn(orig_sequence, data.plist, ffname, pairdist);

        /* print unpaired probabilities */
        if (unpaired > 0) {
          if (plexoutput) {
            pUfp = fopen(fname3, "w");
            putoutphakim_u(fc, data.pup, length, unpaired, pUfp);
            fclose(pUfp);
          }

          /* print unpaired probabilities to file */

          data.pUfp = fopen(openenergies ? fname4 : fname1, "w");
          if (binaries) {
            print_pu_bin(fc, &data, unpaired);
          } else {
            prepare_up_file(&data);
            if (openenergies) {
              for (i = 1; i <= length; i++)
                print_up_open(data.pUfp,
                              i,
                              data.pup[i],
                              (i > unpaired) ? unpaired : i,
                              unpaired,
                              data.kT / 1000.);
            } else {
              for (i = 1; i <= length; i++)
                print_up(data.pUfp, i, data.pup[i], (i > unpaired) ? unpaired : i, unpaired);
            }
          }

          fclose(data.pUfp);
          data.pUfp = NULL;

          for (i = 0; i <= length; i++)
            free(data.pup[i]);
          free(data.pup);
        }
      }

      vrna_fold_compound_free(fc);

      free(pf_parameters);

      /* clean up data */
      if (data.pUfp)
        fclose(data.pUfp);

      if (data.spup)
        fclose(data.spup);

      free(data.plist);


      free(fname1);
      free(fname2);
      free(fname3);
      free(fname4);
      free(ffname);
    }

    (void)fflush(stdout);

    /* clean up */
    free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    free(structure);
    free(SEQ_ID);

    if (with_shapes)
      break;

    rec_id    = rec_sequence = orig_sequence = NULL;
    rec_rest  = NULL;

    /* print user help for the next round if we get input from tty */
    if (istty)
      vrna_message_input_seq_simple();
  }

rnaplfold_exit:

  free(filename_delim);
  free(command_file);
  free(shape_method);
  free(shape_conversion);
  vrna_commands_free(commands);

  free_id_data(id_control);

  return EXIT_SUCCESS;
}


PRIVATE void
print_pu_bin(vrna_fold_compound_t *fc,
             plfold_data          *data,
             int                  ulength)
{
  unsigned int  length;
  int           i, k, *p;
  double        kT = fc->exp_params->kT / 1000.0;

  length = fc->length;

  p = (int *)vrna_alloc(sizeof(int));

  /* write first line */
  p[0] = ulength; /* u length */
  fwrite(p, sizeof(int), 1, data->pUfp);
  p[0] = length;  /* seq length */
  fwrite(p, sizeof(int), 1, data->pUfp);
  for (i = 3; i <= (length + 20); i++) {
    /* all the other lines are set to 1000000 because we are at ulength=0 */
    p[0] = 1000000;
    fwrite(p, sizeof(int), 1, data->pUfp);
  }

  /* write data */
  for (i = 1; i <= ulength; i++) {
    for (k = 1; k <= 11; k++) {
      /* write first ten entries to 1000000 */
      p[0] = 1000000;
      fwrite(p, sizeof(int), 1, data->pUfp);
    }
    for (k = 1; k <= length; k++) {
      /* write data now */
      if (i > k) {
        p[0] = 1000000;         /* check if u > pos */
        fwrite(p, sizeof(int), 1, data->pUfp);
        continue;
      } else {
        p[0] = (int)rint(100 * (-log(data->pup[k][i]) * kT));
        fwrite(p, sizeof(int), 1, data->pUfp);
      }
    }
    for (k = 1; k <= 9; k++) {
      /* finish by writing the last 10 entries */
      p[0] = 1000000;
      fwrite(p, sizeof(int), 1, data->pUfp);
    }
  }
  free(p);
}


PRIVATE void
prepare_up_file(plfold_data *data)
{
  int i;

  if (data->openenergies)
    fprintf(data->pUfp, "#opening energies\n #i$\tl=");
  else
    fprintf(data->pUfp, "#unpaired probabilities\n #i$\tl=");

  for (i = 1; i <= data->ulength; i++)
    fprintf(data->pUfp, "%d\t", i);

  fprintf(data->pUfp, "\n");
}


PRIVATE void
plfold_callback(FLT_OR_DBL    *pr,
                int           pr_size,
                int           i,
                int           max,
                unsigned int  type,
                void          *data)
{
  int         cnt;
  plfold_data *d;

  d = (plfold_data *)data;

  if (type & VRNA_PROBS_WINDOW_BPP) {
    if (!d->simply_putout) {
      /* store pair probabilities in plist */
      d->plist = (vrna_ep_t *)vrna_realloc(d->plist,
                                           sizeof(vrna_ep_t) * (d->plist_cnt + pr_size + 1));

      for (cnt = i + 1; cnt <= pr_size; cnt++) {
        if (pr[cnt] >= d->cutoff) {
          d->plist[d->plist_cnt].i    = i;
          d->plist[d->plist_cnt].j    = cnt;
          d->plist[d->plist_cnt].p    = pr[cnt];
          d->plist[d->plist_cnt].type = VRNA_PLIST_TYPE_BASEPAIR;
          (d->plist_cnt)++;
        }
      }

      /* resize list to actual size */
      d->plist = (vrna_ep_t *)vrna_realloc(d->plist, sizeof(vrna_ep_t) * (d->plist_cnt + 1));

      /* add end-marker to last element */
      d->plist[d->plist_cnt].i    = 0;
      d->plist[d->plist_cnt].j    = 0;
      d->plist[d->plist_cnt].p    = 0.;
      d->plist[d->plist_cnt].type = VRNA_PLIST_TYPE_BASEPAIR;
    } else {
      /* print pair probabilities to output file handle */
      for (cnt = i + 1; cnt <= pr_size; cnt++)
        if (pr[cnt] >= d->cutoff)
          fprintf(d->spup, "%d  %d  %g\n", i, cnt, pr[cnt]);
    }
  }

  /* limit output to full unpaired probabilities */
  if ((type & VRNA_PROBS_WINDOW_UP) && ((type & VRNA_ANY_LOOP) == VRNA_ANY_LOOP)) {
    if (!d->simply_putout) {
      /* store unpaired probabilities in an array */

      /* first allocate some memory */
      d->pup[i]     = (double *)vrna_realloc(d->pup[i], sizeof(double) * (max + 1));
      d->pup[i][0]  = 0.;
      /* copy over unpaired probabilities */
      for (cnt = 1; cnt <= pr_size; cnt++)
        d->pup[i][cnt] = pr[cnt];
      for (cnt = pr_size + 1; cnt <= max; cnt++)
        d->pup[i][cnt] = 0.;
    } else {
      /* print unpaired probabilities to output file handle */
      if (d->openenergies)
        print_up_open(d->pUfp, i, pr, pr_size, max, d->kT / 1000.);
      else
        print_up(d->pUfp, i, pr, pr_size, max);
    }
  }
}


PRIVATE void
print_up_open(FILE    *fp,
              int     i,
              double  *pr,
              int     pr_size,
              int     ulength,
              double  kT)
{
  int cnt;

  fprintf(fp, "%d\t", i);
  for (cnt = 1; cnt < pr_size; cnt++) {
    if (isnan(pr[cnt]) || (pr[cnt] == 0.))
      fprintf(fp, "NA\t");
    else
      fprintf(fp, "%.7g\t", -log(pr[cnt]) * kT);
  }

  if (isnan(pr[pr_size]) || (pr[pr_size] == 0.))
    fprintf(fp, "NA");
  else
    fprintf(fp, "%.7g", -log(pr[pr_size]) * kT);

  for (cnt = pr_size + 1; cnt <= ulength; cnt++)
    fprintf(fp, "\tNA");
  fprintf(fp, "\n");
}


PRIVATE void
print_up(FILE   *fp,
         int    i,
         double *pr,
         int    pr_size,
         int    ulength)
{
  int cnt;

  fprintf(fp, "%d\t", i);
  for (cnt = 1; cnt < pr_size; cnt++) {
    if (isnan(pr[cnt]))
      fprintf(fp, "NA\t");
    else
      fprintf(fp, "%.7g\t", pr[cnt]);
  }

  if (isnan(pr[pr_size]))
    fprintf(fp, "NA");
  else
    fprintf(fp, "%.7g", pr[pr_size]);

  for (cnt = pr_size + 1; cnt <= ulength; cnt++)
    fprintf(fp, "\tNA");
  fprintf(fp, "\n");
}


PRIVATE void
putoutphakim_u(vrna_fold_compound_t *fc,
               double               **pU,
               int                  length,
               int                  ulength,
               FILE                 *fp)
{
  /*put out Fopen in dekacalories per mol, and F(cond,open) also in dekacal*/
  int   k;

  float RT = fc->exp_params->kT;
  float p0;
  float pdep;
  int   f0;
  int   fdep;

  fprintf(fp,
          "#energy necessary to unpair as well as to unpair if i-1 is unpaired also, if i+1 is unpaired also in dekacal/mol\n");
  for (k = 1; k <= length; k++) {
    fprintf(fp, "%d\t", k);
    p0  = pU[k][1];
    f0  = (int)-RT *log(p0) / 10;


    fprintf(fp, "%d\t", f0);
    if (k > 1) {
      pdep  = pU[k][2] / pU[k - 1][1];
      fdep  = (int)-RT *log(pdep) / 10;


      fprintf(fp, "%d\t", fdep);
    } else {
      fprintf(fp, "-0\t");
    }

    if (k < length) {
      pdep  = pU[k + 1][2] / pU[k + 1][1];
      fdep  = (int)-RT *log(pdep) / 10;


      fprintf(fp, "%d\t", fdep);
    } else {
      fprintf(fp, "-0\t");
    }

    fprintf(fp, "\n");
  }

  fflush(fp);
}
