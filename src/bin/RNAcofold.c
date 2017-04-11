/*
 *                c Ivo L Hofacker, Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/commands.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func_co.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "RNAcofold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helper.h"

#include "ViennaRNA/color_output.inc"

PRIVATE vrna_dimer_pf_t do_partfunc(char              *string,
                                    int               length,
                                    int               Switch,
                                    plist             **tpr,
                                    plist             **mf,
                                    vrna_exp_param_t  *parameters);


PRIVATE double *read_concentrations(FILE *fp);


PRIVATE void            do_concentrations(double            FEAB,
                                          double            FEAA,
                                          double            FEBB,
                                          double            FEA,
                                          double            FEB,
                                          double            *startconces,
                                          vrna_exp_param_t  *parameters);


PRIVATE double bppmThreshold;

/*--------------------------------------------------------------------------*/

int
main(int  argc,
     char *argv[])
{
  struct        RNAcofold_args_info args_info;
  char                              *constraints_file, *structure, *cstruc, *rec_sequence, *orig_sequence,
                                    *rec_id, **rec_rest, *Concfile, *id_prefix,
                                    *command_file, *id_delim, *filename_delim, *tmp_string;
  unsigned int                      rec_type, read_opt;
  int                               i, length, cl, pf, istty, noconv, noPS, enforceConstraints,
                                    doT, doC, cofi, auto_id, id_digits, istty_in, istty_out, batch,
                                    filename_full;
  long int                          seq_number;
  double                            min_en, kT, *ConcAandB;
  plist                             *prAB, *prAA, *prBB, *prA, *prB, *mfAB, *mfAA, *mfBB, *mfA, *mfB;
  vrna_md_t                         md;
  vrna_cmd_t                        *commands;

  /*
   #############################################
   # init variables and parameter options
   #############################################
   */
  bppmThreshold = 1e-5;
  noconv        = 0;
  noPS          = 0;
  do_backtrack  = 1;
  pf            = 0;
  doT           = 0;        /* compute dimer free energies etc. */
  doC           = 0;        /* toggle to compute concentrations */
  cofi          = 0;        /* toggle concentrations stdin / file */
  Concfile      = NULL;
  structure     = NULL;
  cstruc        = NULL;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  auto_id       = 0;
  command_file  = NULL;
  commands      = NULL;
  filename_full = 0;

  set_model_details(&md);
  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (RNAcofold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, md);
  ggo_get_md_fold(args_info, md);
  ggo_get_md_part(args_info, md);

  /* check dangle model */
  if ((md.dangles < 0) || (md.dangles > 3)) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    md.dangles = dangles = 2;
  }

  ggo_get_ID_manipulation(args_info,
                          auto_id,
                          id_prefix, "sequence",
                          id_delim, "_",
                          id_digits, 4,
                          seq_number, 1);

  ggo_get_constraints_settings(args_info,
                               fold_constrained,
                               constraints_file,
                               enforceConstraints,
                               batch);

  /* enforce canonical base pairs in any case? */
  if (args_info.canonicalBPonly_given)
    md.canonicalBPonly = canonicalBPonly = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    noconv = 1;

  /*  */
  if (args_info.noPS_given)
    noPS = 1;

  /* concentrations from stdin */
  if (args_info.concentrations_given)
    doC = doT = pf = 1;

  /* set the bppm threshold for the dotplot */
  if (args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0., args_info.bppmThreshold_arg));

  /* concentrations in file */
  if (args_info.concfile_given) {
    Concfile  = strdup(args_info.concfile_arg);
    doC       = cofi = doT = pf = 1;
  }

  /* partition function settings */
  if (args_info.partfunc_given) {
    pf = 1;
    if (args_info.partfunc_arg != -1)
      md.compute_bpp = do_backtrack = args_info.partfunc_arg;
  }

  if (args_info.all_pf_given) {
    doT = pf = 1;
    if (args_info.all_pf_arg != 1)
      md.compute_bpp = do_backtrack = args_info.all_pf_arg;
    else
      md.compute_bpp = do_backtrack = 1;
  }

  if (args_info.commands_given)
    command_file = strdup(args_info.commands_arg);

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    filename_delim = strdup(args_info.filename_delim_arg);
  else
    filename_delim = strdup(id_delim);

  if (isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

  /* full filename from FASTA header support */
  if (args_info.filename_full_given)
    filename_full = 1;

  /* free allocated memory of command line data structure */
  RNAcofold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (pf && gquad)
    vrna_message_error("G-Quadruplex support is currently not available for partition function computations");

  if (command_file != NULL)
    commands = vrna_file_commands_read(command_file, VRNA_CMD_PARSE_DEFAULTS);

  istty_in  = isatty(fileno(stdin));
  istty_out = isatty(fileno(stdout));
  istty     = isatty(fileno(stdout)) && isatty(fileno(stdin));

  /* print user help if we get input from tty */
  if (istty) {
    if (fold_constrained) {
      vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK);
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n"
                             "Use '&' to connect 2 sequences that shall form a complex.");
    } else {
      vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.");
    }
  }

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty)
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;

  if (!fold_constrained)
    read_opt |= VRNA_INPUT_NO_REST;

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
    char  *SEQ_ID         = NULL;
    int   maybe_multiline = 0;

    if (rec_id) {
      maybe_multiline = 1;
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
    }

    /* construct the sequence ID */
    ID_generate(SEQ_ID, rec_id, auto_id, id_prefix, id_delim, id_digits, seq_number, filename_full);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(rec_sequence);

    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound_t *vc = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_HYBRID | ((pf) ? VRNA_OPTION_PF : 0));
    length    = vc->length;
    structure = (char *)vrna_alloc((unsigned)length + 1);

    /* parse the rest of the current dataset to obtain a structure constraint */
    if (fold_constrained) {
      if (constraints_file) {
        vrna_constraints_add(vc, constraints_file, VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));
      } else {
        cstruc  = NULL;
        cstruc  = NULL;
        int           cp        = -1;
        unsigned int  coptions  = (maybe_multiline) ? VRNA_OPTION_MULTILINE : 0;
        cstruc  = vrna_extract_record_rest_structure((const char **)rec_rest, 0, coptions);
        cstruc  = vrna_cut_point_remove(cstruc, &cp);
        if (vc->cutpoint != cp) {
          vrna_message_error("Sequence and Structure have different cut points.\n"
                             "sequence: %d, structure: %d",
                             vc->cutpoint, cp);
        }

        cl = (cstruc) ? (int)strlen(cstruc) : 0;

        if (cl == 0)
          vrna_message_warning("Structure constraint is missing");
        else if (cl < length)
          vrna_message_warning("Structure constraint is shorter than sequence");
        else if (cl > length)
          vrna_message_error("Structure constraint is too long");

        if (cstruc) {
          strncpy(structure, cstruc, sizeof(char) * (cl + 1));

          unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;
          if (enforceConstraints)
            constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

          vrna_constraints_add(vc, (const char *)structure, constraint_options);
        }
      }
    }

    if (commands)
      vrna_commands_apply(vc, commands, VRNA_CMD_PARSE_DEFAULTS);

    if (istty) {
      if (cut_point == -1)
        vrna_message_info(stdout, "length = %d", length);
      else
        vrna_message_info(stdout, "length1 = %d\nlength2 = %d", cut_point - 1, length - cut_point + 1);
    }

    if (doC) {
      FILE *fp;
      if (cofi) {
        /* read from file */
        fp = fopen(Concfile, "r");
        if (fp == NULL)
          vrna_message_error("could not open concentration file %s", Concfile);

        ConcAandB = read_concentrations(fp);
        fclose(fp);
      } else {
        printf("Please enter concentrations [mol/l]\n format: ConcA ConcB\n return to end\n");
        ConcAandB = read_concentrations(stdin);
      }
    }

    /*
     ########################################################
     # begin actual computations
     ########################################################
     */

    /* compute mfe of AB dimer */
    min_en  = vrna_mfe_dimer(vc, structure);
    mfAB    = vrna_plist(structure, 0.95);

    /* check whether the constraint allows for any solution */
    if (fold_constrained && constraints_file) {
      if (min_en == (double)(INF / 100.)) {
        vrna_message_error("Supplied structure constraints create empty solution set for sequence:\n%s",
                           orig_sequence);
        exit(EXIT_FAILURE);
      }
    }

    {
      char *pstring, *pstruct, *msg = NULL;
      pstring = strdup(orig_sequence);
      pstruct = vrna_cut_point_insert(structure, vc->cutpoint);

      print_fasta_header(stdout, rec_id);
      fprintf(stdout, "%s\n", orig_sequence);

      if (istty)
        msg = vrna_strdup_printf("\n minimum free energy = %6.2f kcal/mol", min_en);
      else
        msg = vrna_strdup_printf(" (%6.2f)", min_en);

      print_structure(stdout, pstruct, msg);
      (void)fflush(stdout);

      if (!noPS) {
        char *filename_plot = NULL, *annot = NULL;
        if (SEQ_ID) {
          filename_plot = vrna_strdup_printf("%s%sss.ps", SEQ_ID, id_delim);
          tmp_string    = vrna_filename_sanitize(filename_plot, filename_delim);
          free(filename_plot);
          filename_plot = tmp_string;
        } else {
          filename_plot = strdup("rna.ps");
        }

        if (vc->cutpoint >= 0) {
          annot = vrna_strdup_printf("1 %d 9  0 0.9 0.2 omark\n"
                                     "%d %d 9  1 0.1 0.2 omark\n",
                                     vc->cutpoint - 1,
                                     vc->cutpoint + 1,
                                     length + 1);
        }

        if (filename_plot)
          (void)vrna_file_PS_rnaplot_a(pstring, pstruct, filename_plot, annot, NULL, &md);

        free(filename_plot);
        free(annot);
      }

      free(pstring);
      free(pstruct);
      free(msg);
    }

    if (length > 2000)
      vrna_mx_mfe_free(vc);

    /* compute partition function */
    if (pf) {
      vrna_dimer_pf_t AB, AA, BB;

      prAB = NULL;
      prAA = NULL;
      prBB = NULL;
      prA  = NULL;
      prB  = NULL;

      if (md.dangles == 1) {
        vc->params->model_details.dangles = dangles = 2;   /* recompute with dangles as in pf_fold() */
        min_en                            = vrna_eval_structure(vc, structure);
        vc->params->model_details.dangles = dangles = 1;
      }

      vrna_exp_params_rescale(vc, &min_en);
      kT = vc->exp_params->kT / 1000.;

      if (length > 2000)
        vrna_message_info(stderr, "scaling factor %f", vc->exp_params->pf_scale);

      /* do we need to add hard constraints? */
      if (cstruc != NULL)
        strncpy(structure, cstruc, length + 1);

      /* compute partition function */
      AB = vrna_pf_dimer(vc, structure);

      if (do_backtrack) {
        char *costruc, *msg = NULL;
        costruc = vrna_cut_point_insert(structure, vc->cutpoint);
        if (istty_in)
          msg = vrna_strdup_printf("\n free energy of ensemble = %6.2f kcal/mol", AB.FAB);
        else
          msg = vrna_strdup_printf(" [%6.2f]", AB.FAB);

        print_structure(stdout, costruc, msg);
        free(msg);
        free(costruc);
        prAB = vrna_plist_from_probs(vc, bppmThreshold);
      } else {
        char *msg = vrna_strdup_printf(" free energy of ensemble = %6.2f kcal/mol", AB.FAB);
        print_structure(stdout, NULL, msg);
        free(msg);
      }

      {
        char *msg = vrna_strdup_printf(" frequency of mfe structure in ensemble %g"
                                       "; delta G binding=%6.2f",
                                       exp((AB.FAB - min_en) / kT),
                                       AB.FcAB - AB.FA - AB.FB);
        print_structure(stdout, NULL, msg);
        free(msg);
      }

      /* if (doQ) make_probsum(length,fname); */ /*compute prob of base paired*/
      /* free_co_arrays(); */
      if (doT) {
        /* cofold of all dimers, monomers */
        int   Blength, Alength;
        char  *Astring, *Bstring, *orig_Astring, *orig_Bstring;
        char  *Newstring;
        if (vc->cutpoint <= 0) {
          vrna_message_warning("Sorry, i cannot do that with only one molecule, please give me two or leave it");
          free(mfAB);
          free(prAB);
          continue;
        }

        if (dangles == 1)
          dangles = 2;

        Alength = vc->cutpoint - 1;                                 /* length of first molecule */
        Blength = length - vc->cutpoint + 1;                        /* length of 2nd molecule   */

        Astring = (char *)vrna_alloc(sizeof(char) * (Alength + 1)); /*Sequence of first molecule*/
        Bstring = (char *)vrna_alloc(sizeof(char) * (Blength + 1)); /*Sequence of second molecule*/
        strncat(Astring, rec_sequence, Alength);
        strncat(Bstring, rec_sequence + Alength + 1, Blength);

        orig_Astring  = (char *)vrna_alloc(sizeof(char) * (Alength + 1)); /*Sequence of first molecule*/
        orig_Bstring  = (char *)vrna_alloc(sizeof(char) * (Blength + 1)); /*Sequence of second molecule*/
        strncat(orig_Astring, orig_sequence, Alength);
        strncat(orig_Bstring, orig_sequence + Alength + 1, Blength);

        /* compute AA dimer */
        AA = do_partfunc(Astring, Alength, 2, (do_backtrack) ? &prAA : NULL, &mfAA, vc->exp_params);
        /* compute BB dimer */
        BB = do_partfunc(Bstring, Blength, 2, (do_backtrack) ? &prBB : NULL, &mfBB, vc->exp_params);
        /*free_co_pf_arrays();*/

        /* compute A monomer */
        do_partfunc(Astring, Alength, 1, (do_backtrack) ? &prA : NULL, &mfA, vc->exp_params);

        /* compute B monomer */
        do_partfunc(Bstring, Blength, 1, (do_backtrack) ? &prB : NULL, &mfB, vc->exp_params);

        if (do_backtrack) {
          vrna_pf_dimer_probs(AB.F0AB, AB.FA, AB.FB, prAB, prA, prB, Alength, vc->exp_params);
          vrna_pf_dimer_probs(AA.F0AB, AA.FA, AA.FA, prAA, prA, prA, Alength, vc->exp_params);
          vrna_pf_dimer_probs(BB.F0AB, BB.FA, BB.FA, prBB, prA, prB, Blength, vc->exp_params);
        }

        print_comment(stdout, "Free Energies:");
        char *thead = NULL, *tline = NULL;
        thead = strdup("AB\t\tAA\t\tBB\t\tA\t\tB");
        tline = vrna_strdup_printf("%.6f\t%6f\t%6f\t%6f\t%6f",
                                   AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB);
        print_table(stdout, thead, tline);
        free(thead);
        free(tline);

        if (doC) {
          do_concentrations(AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB, ConcAandB, vc->exp_params);
          free(ConcAandB);/*freeen*/
        }

        /*output of the 5 dot plots*/

        if (do_backtrack) {
          char  *comment    = NULL;
          char  *fname_dot  = NULL;
          char *filename_dot = NULL;

          if (SEQ_ID) {
            filename_dot  = vrna_strdup_printf("%s%sdp5.ps", SEQ_ID, id_delim);
            tmp_string    = vrna_filename_sanitize(filename_dot, filename_delim);
            free(filename_dot);
            filename_dot = tmp_string;
          } else {
            filename_dot = strdup("dot5.ps");
          }

          /*AB dot_plot*/

          /*write Free Energy into comment*/
          comment = vrna_strdup_printf("\n%%Heterodimer AB FreeEnergy= %.9f\n", AB.FcAB);
          /*reset cut_point*/
          cut_point = Alength + 1;

          /*write New name*/
          fname_dot = vrna_strdup_printf("AB%s", filename_dot);
          int   cp      = -1;
          char  *coseq  = vrna_cut_point_remove(orig_sequence, &cp);
          (void)PS_dot_plot_list(coseq, fname_dot, prAB, mfAB, comment);
          free(coseq);
          free(fname_dot);
          free(comment);
          fname_dot = NULL;
          comment   = NULL;

          /*AA dot_plot*/
          comment = vrna_strdup_printf("\n%%Homodimer AA FreeEnergy= %.9f\n", AA.FcAB);
          /*write New name*/
          fname_dot = vrna_strdup_printf("AA%s", filename_dot);
          /*write AA sequence*/
          Newstring = (char *)vrna_alloc((2 * Alength + 1) * sizeof(char));
          strcpy(Newstring, orig_Astring);
          strcat(Newstring, orig_Astring);
          (void)PS_dot_plot_list(Newstring, fname_dot, prAA, mfAA, comment);
          free(Newstring);
          free(fname_dot);
          free(comment);
          fname_dot = NULL;
          comment   = NULL;

          /*BB dot_plot*/
          comment = vrna_strdup_printf("\n%%Homodimer BB FreeEnergy= %.9f\n", BB.FcAB);
          /*write New name*/
          fname_dot = vrna_strdup_printf("BB%s", filename_dot);
          /*write BB sequence*/
          Newstring = (char *)vrna_alloc((2 * Blength + 1) * sizeof(char));
          strcpy(Newstring, orig_Bstring);
          strcat(Newstring, orig_Bstring);
          /*reset cut_point*/
          cut_point = Blength + 1;
          (void)PS_dot_plot_list(Newstring, fname_dot, prBB, mfBB, comment);
          free(Newstring);
          free(fname_dot);
          free(comment);
          fname_dot = NULL;
          comment   = NULL;

          /*A dot plot*/
          /*reset cut_point*/
          cut_point = -1;
          comment   = vrna_strdup_printf("\n%%Monomer A FreeEnergy= %.9f\n", AB.FA);
          /*write New name*/
          fname_dot = vrna_strdup_printf("A%s", filename_dot);
          /*write BB sequence*/
          (void)PS_dot_plot_list(orig_Astring, fname_dot, prA, mfA, comment);
          free(fname_dot);
          free(comment);
          fname_dot = NULL;
          comment   = NULL;

          /*B monomer dot plot*/
          comment = vrna_strdup_printf("\n%%Monomer B FreeEnergy= %.9f\n", AB.FB);
          /*write New name*/
          fname_dot = vrna_strdup_printf("B%s", filename_dot);
          /*write BB sequence*/
          (void)PS_dot_plot_list(orig_Bstring, fname_dot, prB, mfB, comment);
          free(fname_dot);
          free(comment);
          free(filename_dot);
        }

        free(Astring);
        free(Bstring);
        free(orig_Astring);
        free(orig_Bstring);
        free(prAB);
        free(prAA);
        free(prBB);
        free(prA);
        free(prB);
        free(mfAB);
        free(mfAA);
        free(mfBB);
        free(mfA);
        free(mfB);
      } /*end if(doT)*/
    }   /*end if(pf)*/

    if (do_backtrack) {
      if (!doT) {
        if (pf) {
          int   cp;
          char  *seq          = vrna_cut_point_remove(rec_sequence, &cp);
          char  *filename_dot = NULL;
          if (SEQ_ID) {
            filename_dot  = vrna_strdup_printf("%s%sdp.ps", SEQ_ID, id_delim);
            tmp_string    = vrna_filename_sanitize(filename_dot, filename_delim);
            free(filename_dot);
            filename_dot = tmp_string;
          } else {
            filename_dot = strdup("dot.ps");
          }

          if (filename_dot)
            (void)vrna_plot_dp_PS_list(seq, cp, filename_dot, prAB, mfAB, "doof");

          free(filename_dot);
          free(prAB);
          free(seq);
        }

        free(mfAB);
      }
    }

    if (!doT)
      vrna_mx_pf_free(vc);

    (void)fflush(stdout);

    /* clean up */
    free(cstruc);
    free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    free(structure);
    /* free the rest of current dataset */
    if (rec_rest) {
      for (i = 0; rec_rest[i]; i++)
        free(rec_rest[i]);
      free(rec_rest);
    }

    rec_id    = rec_sequence = orig_sequence = structure = cstruc = NULL;
    rec_rest  = NULL;
    vrna_fold_compound_free(vc);

    free(SEQ_ID);

    if (constraints_file && (!batch))
      break;

    ID_number_increase(seq_number, "Sequence");

    /* print user help for the next round if we get input from tty */
    if (istty) {
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      if (fold_constrained) {
        vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK | VRNA_CONSTRAINT_DB_RND_BRACK);
        vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n");
      } else {
        vrna_message_input_seq_simple();
      }
    }
  }

  free(command_file);
  vrna_commands_free(commands);

  free(id_prefix);
  free(id_delim);
  free(filename_delim);
  free(Concfile);

  return EXIT_SUCCESS;
}


PRIVATE vrna_dimer_pf_t
do_partfunc(char              *string,
            int               length,
            int               Switch,
            plist             **tpr,
            plist             **mfpl,
            vrna_exp_param_t  *parameters)
{
  /*compute mfe and partition function of dimer or monomer*/
  char                  *Newstring;
  char                  *tempstruc;
  double                min_en;
  double                sfact = 1.07;
  double                kT;
  vrna_exp_param_t      *par;
  vrna_dimer_pf_t       X;
  vrna_fold_compound_t  *vc;

  kT = parameters->kT / 1000.;
  switch (Switch) {
    case 1:   /* monomer */
      tempstruc = (char *)vrna_alloc((unsigned)length + 1);
      //parameters->model_details.min_loop_size = TURN; /* we need min_loop_size of 0 to correct for Q_AB */
      vc      = vrna_fold_compound(string, &(parameters->model_details), VRNA_OPTION_MFE | VRNA_OPTION_PF);
      min_en  = vrna_mfe(vc, tempstruc);
      *mfpl   = vrna_plist(tempstruc, 0.95);
      vrna_mx_mfe_free(vc);

      par           = vrna_exp_params_copy(parameters);
      par->pf_scale = exp(-(sfact * min_en) / kT / (length));
      vrna_exp_params_subst(vc, par);
      X = vrna_pf_dimer(vc, tempstruc);
      if (tpr)
        *tpr = vrna_plist_from_probs(vc, bppmThreshold);

      vrna_fold_compound_free(vc);
      free(tempstruc);
      free(par);
      parameters->model_details.min_loop_size = 0;
      break;

    case 2:   /* dimer */
      tempstruc = (char *)vrna_alloc((unsigned)length * 2 + 2);
      Newstring = (char *)vrna_alloc(sizeof(char) * (length * 2 + 2));
      strcat(Newstring, string);
      strcat(Newstring, "&");
      strcat(Newstring, string);
      parameters->model_details.min_loop_size = 0;
      vc                                      = vrna_fold_compound(Newstring, &(parameters->model_details), VRNA_OPTION_MFE | VRNA_OPTION_PF | VRNA_OPTION_HYBRID);
      min_en                                  = vrna_mfe_dimer(vc, tempstruc);
      *mfpl                                   = vrna_plist(tempstruc, 0.95);
      vrna_mx_mfe_free(vc);

      par           = vrna_exp_params_copy(parameters);
      par->pf_scale = exp(-(sfact * min_en) / kT / (2 * length));
      vrna_exp_params_subst(vc, par);
      X = vrna_pf_dimer(vc, tempstruc);
      if (tpr)
        *tpr = vrna_plist_from_probs(vc, bppmThreshold);

      vrna_fold_compound_free(vc);

      free(Newstring);
      free(tempstruc);
      free(par);
      break;

    default:
      printf("Error in get_partfunc\n, computing neither mono- nor dimere!\n");
      exit(42);
  }

  return X;
}


PRIVATE void
do_concentrations(double            FEAB,
                  double            FEAA,
                  double            FEBB,
                  double            FEA,
                  double            FEB,
                  double            *startconc,
                  vrna_exp_param_t  *parameters)
{
  /* compute and print concentrations out of free energies, calls get_concentrations */
  vrna_dimer_conc_t *result;
  int               i, n;

  result = vrna_pf_dimer_concentrations(FEAB, FEAA, FEBB, FEA, FEB, startconc, parameters);

  print_table(stdout, "Initial concentrations\t\trelative Equilibrium concentrations", NULL);
  print_table(stdout, "A\t\tB\t\tAB\t\tAA\t\tBB\t\tA\t\tB", NULL);
  for (n = 0; (startconc[2 * n] > 0) || (startconc[2 * n + 1] > 0); n++); /* count */
  for (i = 0; i < n; i++) {
    double  tot     = result[i].Ac_start + result[i].Bc_start;
    char    *tline  = vrna_strdup_printf("%-10g\t%-10g\t%.5f \t%.5f \t%.5f \t%.5f \t%.5f",
                                         result[i].Ac_start,
                                         result[i].Bc_start,
                                         result[i].ABc / tot,
                                         result[i].AAc / tot,
                                         result[i].BBc / tot,
                                         result[i].Ac / tot,
                                         result[i].Bc / tot);
    print_table(stdout, NULL, tline);
    free(tline);
  }
  free(result);
}


PRIVATE double *
read_concentrations(FILE *fp)
{
  /* reads concentrations, returns list of double, -1. marks end */
  char    *line;
  double  *startc;
  int     i = 0, n = 2;

  startc = (double *)vrna_alloc((2 * n + 1) * sizeof(double));

  while ((line = vrna_read_line(fp)) != NULL) {
    int c;
    if (i + 4 >= 2 * n) {
      n       *= 2;
      startc  = (double *)vrna_realloc(startc, (2 * n + 1) * sizeof(double));
    }

    c = sscanf(line, "%lf %lf", &startc[i], &startc[i + 1]);
    free(line);
    if (c < 2)
      break;

    i += 2;
  }
  startc[i] = startc[i + 1] = 0;
  return startc;
}
