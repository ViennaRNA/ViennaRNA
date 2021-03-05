/*
 *           Compute pseudoknotted structure of an RNA
 *
 *                         c Ivo L Hofacker
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/units.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/PKplex.h"
#include "RNAPKplex_cmdl.h"

int
PlexHit_cmp(const void  *c1,
            const void  *c2);


int
PlexHit_cmp_energy(const void *c1,
                   const void *c2);

int
PlexHit_cmp_active_energy(const void *c1,
                          const void *c2);


/*--------------------------------------------------------------------------*/

int
main(int  argc,
     char *argv[])
{
  FILE                    *pUfp, *spup;
  struct PKplex_args_info args_info;
  char                    *id_s1, *s1, *orig_s1, *ParamFile, *ns_bases, *plexstring,
                          *constraint, fname[FILENAME_MAX_LENGTH], *annotation, **rest;
  unsigned int            options;
  int                     istty, i, j, noconv, length, pairdist, current, unpaired, winsize, verbose;
  float                   cutoff, constrainedEnergy;
  double                  **pup, subopts, pk_penalty;
  plist                   *pl, *dpp;
  vrna_md_t               md;

  options     = 0;
  pup         = NULL;           /*prob of being unpaired, lengthwise*/
  pUfp        = NULL;
  dpp         = NULL;
  spup        = NULL;
  subopts     = 0.0;
  dangles     = 2;
  winsize     = 70;
  cutoff      = 0.01;
  pairdist    = 0;
  unpaired    = 0;
  noconv      = 0;
  pk_penalty  = 8.10;
  ParamFile   = ns_bases = NULL;
  s1          = id_s1 = orig_s1 = NULL;
  annotation  = NULL;
  verbose     = 0;

  vrna_md_set_default(&md);

  /*
   #############################################
   # check command line parameters
   #############################################
   */
  if (PKplex_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

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

  /* take another energy parameter set */
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* set the pair probability cutoff */
  if (args_info.cutoff_given)
    cutoff = args_info.cutoff_arg;

  /* turn on verbose output (mainly for debugging) */
  if (args_info.verbose_given)
    verbose = 1;

  /* set energy cutoff */
  if (args_info.energyCutoff_given)
    pk_penalty = args_info.energyCutoff_arg;

  /* show suboptimal structures which are better than given value difference */
  if (args_info.subopts_given)
    subopts = args_info.subopts_arg;

  /* free allocated memory of command line data structure */
  PKplex_cmdline_parser_free(&args_info);

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

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
  if (istty)
    options |= VRNA_INPUT_NOSKIP_BLANK_LINES;

  options |= VRNA_INPUT_NO_REST;
  if (istty)
    vrna_message_input_seq_simple();

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  while (!(vrna_file_fasta_read_record(&id_s1, &s1, &rest, NULL,
                                       options) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    /*
     ########################################################
     # handle user input from 'stdin'
     ########################################################
     */
    if (id_s1) {
      if (!istty)
        printf("%s\n", id_s1);

      (void)sscanf(id_s1, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    } else {
      strcpy(fname, "PKplex");
    }

    strcat(fname, ".ps");

    length    = strlen(s1);
    winsize   = pairdist = length;
    unpaired  = MIN2(30, length - 3);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(s1);

    /* store case-unmodified sequence */
    orig_s1 = strdup(s1);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(s1);

    printf("%s\n", orig_s1);
    if (verbose)
      printf("length = %d\n", length);

    /*
     ########################################################
     # do PLfold computations
     ########################################################
     */
    if (length >= 5) {
      vrna_pk_plex_t *hits, *hit_ptr;

      if (verbose)
        printf("Winsize = %d\nPairdist = %d\nUnpaired = %d\n", winsize, pairdist, unpaired);

      /*
       ########################################################
       # do Plex computations
       ########################################################
       */
      int **access = vrna_pk_plex_accessibility(s1, unpaired, cutoff);

      if (verbose)
        printf("EnergyCutoff = %f\n", pk_penalty);

      vrna_fold_compound_t  *fc;
      vrna_pk_plex_opt_t    pk_plex_options;

      fc = vrna_fold_compound(s1, &md, VRNA_OPTION_DEFAULT);


      pk_plex_options = vrna_pk_plex_opt((unsigned int)vrna_convert_kcal_to_dcal(subopts),
                                         MIN2(12, length - 3),
                                         (unsigned int)vrna_convert_kcal_to_dcal(pk_penalty));
      hits = vrna_pk_plex(fc,
                          (const int **)access,
                          pk_plex_options);

      vrna_fold_compound_free(fc);

      /* and print the results to stdout */
      for (hit_ptr = hits; hit_ptr->structure; hit_ptr++) {
          if (verbose) {
            printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n",
                    hit_ptr->structure,
                    hit_ptr->start_5,
                    hit_ptr->end_5,
                    hit_ptr->start_3,
                    hit_ptr->end_3,
                    hit_ptr->energy,
                    hit_ptr->energy - hit_ptr->dGint - hit_ptr->dGpk,
                    hit_ptr->dGint,
                    hit_ptr->dGpk);
          } else {
            printf("%s (%6.2f)\n", hit_ptr->structure, hit_ptr->energy);
          }
      }

      /*
       ########################################################
       # Generate Plot for the best structure
       ########################################################
       */

      /* make an EPS layout file for the MFE structure, i.e. first in the list */
      hit_ptr = hits;
      if (hit_ptr->start_5 > 0) { /* true hit with PK */
        annotation = vrna_strdup_printf("%u %u 13 1 0 0 omark\n"
                                        "%u %u 13 1 0 0 omark\n"
                                        "0 0 2 setrgbcolor\n"
                                        "2 setlinewidth\n"
                                        "%u cmark\n"
                                        "%u cmark\n"
                                        "1 setlinewidth",
                                        hit_ptr->start_5,
                                        hit_ptr->end_5,
                                        hit_ptr->start_3,
                                        hit_ptr->end_3,
                                        hit_ptr->start_5,
                                        hit_ptr->end_3);
        vrna_file_PS_rnaplot_a(s1, hit_ptr->structure, fname, annotation, "", &md);
        free(annotation);
        annotation = NULL;
      } else { /* no PK hit found, let's plot the MFE structure instead */
        vrna_file_PS_rnaplot_a(s1, hit_ptr->structure, fname, NULL, "", &md);
      }

      /*
       ########################################################
       # free memory
       ########################################################
       */
      for (hit_ptr = hits; hit_ptr->structure; hit_ptr++)
        free(hit_ptr->structure);

      free(hits);

      (void)fflush(stdout);
      i = access[0][0];
      while (--i > -1)
        free(access[i]);
      free(access);
      free(constraint);
    }

    free(s1);
    free(orig_s1);
    free(id_s1);
    free(plexstring);
    free(nonstandards);
    s1 = id_s1 = orig_s1 = NULL;

    plexstring = NULL;

    /* print user help for the next round if we get input from tty */
    if (istty)
      vrna_message_input_seq_simple();
  }
  return 0;
}


int
PlexHit_cmp(const void  *c1,
            const void  *c2)
{
  dupVar  *p1 = (dupVar *)c1;
  dupVar  *p2 = (dupVar *)c2;

  return p1->ddG >= p2->ddG;
}


int
PlexHit_cmp_energy(const void *c1,
                   const void *c2)
{
  dupVar  *p1 = (dupVar *)c1;
  dupVar  *p2 = (dupVar *)c2;

  if (p1->energy > p2->energy)
    return 1;
  else if (p1->energy < p2->energy)
    return -1;

  return 0;
}

int
PlexHit_cmp_active_energy(const void *c1,
                          const void *c2)
{
  dupVar  *p1 = (dupVar *)c1;
  dupVar  *p2 = (dupVar *)c2;

  if (p1->inactive > p2->inactive) {
    return 1;
  } else if (p1->inactive < p2->inactive) {
    return -1;
  } else if (!p2->inactive) {
    if (p1->energy < p2->energy)
      return -1;
    else if (p1->energy > p2->energy)
      return 1;
  }

  return 0;
}
