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
#include "ViennaRNA/utils.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/energy_const.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/file_formats.h"
#include "RNAplfold_cmdl.h"

#include "ViennaRNA/color_output.inc"

int unpaired;
PRIVATE void putoutphakim_u(double  **pU,
                            int     length,
                            int     ulength,
                            FILE    *fp);


/*--------------------------------------------------------------------------*/
int
main(int  argc,
     char *argv[])
{
  FILE                        *pUfp, *spup;
  struct RNAplfold_args_info  args_info;
  char                        *fname, *structure, *ParamFile, *ns_bases, *rec_sequence, *rec_id,
                              **rec_rest, *orig_sequence;
  unsigned int                rec_type, read_opt;
  int                         length, istty, winsize, pairdist, tempwin, temppair, tempunpaired, noconv,
                              plexoutput, simply_putout, openenergies, binaries;
  float                       cutoff;
  double                      **pup, betaScale;
  plist                       *pl, *dpp;
  vrna_exp_param_t            *pf_parameters;
  vrna_md_t                   md;

  pUfp          = NULL;
  spup          = NULL;
  pup           = NULL; /*prob of being unpaired, lengthwise*/
  dpp           = NULL;
  dangles       = 2;
  cutoff        = 0.01;
  winsize       = 70;
  pairdist      = 0;
  unpaired      = 0;
  betaScale     = 1.;
  simply_putout = plexoutput = openenergies = noconv = 0;
  binaries      = 0;
  tempwin       = temppair = tempunpaired = 0;
  structure     = ParamFile = ns_bases = NULL;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  pf_parameters = NULL;
  fname         = NULL;

  set_model_details(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAplfold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
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

  if (args_info.betaScale_given)
    md.betaScale = betaScale = args_info.betaScale_arg;

  /* check for errorneous parameter options */
  if ((pairdist < 0) || (cutoff < 0.) || (unpaired < 0) || (winsize < 0)) {
    RNAplfold_cmdline_parser_print_help();
    exit(EXIT_FAILURE);
  }

  /* free allocated memory of command line data structure */
  RNAplfold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL)
    vrna_md_set_nonstandards(&md, ns_bases);

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
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    if (rec_id) {
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
      if (!istty)
        print_fasta_header(stdout, rec_id);

      fname = strdup(rec_id);
      (void)sscanf(rec_id, "%s", fname); /* truncate to first word */
      fname = (char *)vrna_realloc(fname, sizeof(char) * (strlen(fname) + 1));
    }

    length    = (int)strlen(rec_sequence);
    structure = (char *)vrna_alloc((unsigned)length + 1);

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

    if (unpaired && simply_putout) {
      vrna_message_warning("Output simplification not possible if unpaired is switched on");
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

      if ((!fname) || (fname[0] == '\0')) {
        free(fname);
        fname = strdup("plfold");
      }

      fname1  = vrna_strdup_printf("%s_lunp", fname);
      fname2  = vrna_strdup_printf("%s_basepairs", fname);
      fname3  = vrna_strdup_printf("%s_uplex", fname);
      fname4  = (binaries) ? vrna_strdup_printf("%s_openen_bin", fname) : vrna_strdup_printf("%s_openen", fname);
      ffname  = vrna_strdup_printf("%s_dp.ps", fname);

      /* sanitize filenames */
      tmp_string = vrna_filename_sanitize(fname1, "_");
      free(fname1);
      fname1      = tmp_string;
      tmp_string  = vrna_filename_sanitize(fname2, "_");
      free(fname2);
      fname2      = tmp_string;
      tmp_string  = vrna_filename_sanitize(fname3, "_");
      free(fname3);
      fname3      = tmp_string;
      tmp_string  = vrna_filename_sanitize(fname4, "_");
      free(fname4);
      fname4      = tmp_string;
      tmp_string  = vrna_filename_sanitize(ffname, "_");
      free(ffname);
      ffname = tmp_string;

      pf_parameters = vrna_exp_params(&md);

      if (unpaired > 0) {
        pup       = (double **)vrna_alloc((length + 1) * sizeof(double *));
        pup[0]    = (double *)vrna_alloc(sizeof(double));   /*I only need entry 0*/
        pup[0][0] = unpaired;
      }

      pUfp = spup = NULL;
      if (simply_putout) {
        spup  = fopen(fname2, "w");
        pUfp  = (unpaired > 0) ? fopen(fname1, "w") : NULL;

        pl = pfl_fold_par(rec_sequence, winsize, pairdist, cutoff, pup, &dpp, pUfp, spup, pf_parameters);

        if (pUfp != NULL)
          fclose(pUfp);

        if (spup != NULL)
          fclose(spup);
      } else {
        pl = pfl_fold_par(rec_sequence, winsize, pairdist, cutoff, pup, &dpp, pUfp, spup, pf_parameters);
        PS_dot_plot_turn(orig_sequence, pl, ffname, pairdist);
        if (unpaired > 0) {
          if (plexoutput) {
            pUfp = fopen(fname3, "w");
            putoutphakim_u(pup, length, unpaired, pUfp);
            fclose(pUfp);
          }

          pUfp = fopen(openenergies ? fname4 : fname1, "w");
          if (binaries)
            putoutpU_prob_bin_par(pup, length, unpaired, pUfp, openenergies, pf_parameters);
          else
            putoutpU_prob_par(pup, length, unpaired, pUfp, openenergies, pf_parameters);

          fclose(pUfp);
        }
      }

      free(pl);

      if (unpaired > 0) {
        free(pup[0]);
        free(pup);
      }

      free(pf_parameters);
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
    free(fname);
    fname     = NULL;
    rec_id    = rec_sequence = orig_sequence = NULL;
    rec_rest  = NULL;
    /* print user help for the next round if we get input from tty */

    if (istty)
      vrna_message_input_seq_simple();
  }
  return EXIT_SUCCESS;
}


PRIVATE void
putoutphakim_u(double **pU,
               int    length,
               int    ulength,
               FILE   *fp)
{
  /*put out Fopen in dekacalories per mol, and F(cond,open) also in dekacal*/
  int   k;

  float RT = (temperature + K0) * GASCONST;
  float p0;
  float pdep;
  int   f0;
  int   fdep;

  fprintf(fp, "#energy necessary to unpair as well as to unpair if i-1 is unpaired also, if i+1 is unpaired also in dekacal/mol\n");
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
