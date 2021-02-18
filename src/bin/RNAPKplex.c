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
      vrna_pkplex_t *hits, *hit_ptr;
      pf_scale = -1;

      pup       = (double **)vrna_alloc((length + 1) * sizeof(double *));
      pup[0]    = (double *)vrna_alloc(sizeof(double));   /*I only need entry 0*/
      pup[0][0] = unpaired;

      pUfp = spup = NULL;

      if (verbose)
        printf("Winsize = %d\nPairdist = %d\nUnpaired = %d\n", winsize, pairdist, unpaired);

      int tempdangles = dangles;
      dangles = 2;
      pl      = pfl_fold(s1, winsize, pairdist, cutoff, pup, &dpp, pUfp, spup);
      dangles = tempdangles;

      /*
       ########################################################
       # do Plex computations
       ########################################################
       */
      double  kT = (temperature + K0) * GASCONST / 1000.0;
      int     **access;
      /* prepare the accesibility array */
      access = (int **)vrna_alloc(sizeof(int *) * (unpaired + 2));
      for (i = 0; i < unpaired + 2; i++)
        access[i] = (int *)vrna_alloc(sizeof(int) * (length + 20));

      for (i = 0; i < length + 20; i++)
        for (j = 0; j < unpaired + 2; j++)
          access[j][i] = INF;

      for (i = 11; i < length + 11; i++) {
        for (j = 1; j < unpaired + 1; j++)
          if (pup[i - 10][j] > 0)
            access[j][i] = rint(100 * (-log(pup[i - 10][j])) * kT);
      }

      access[0][0] = unpaired + 2;

      plexstring = (char *)vrna_alloc(length + 1 + 20);
      strcpy(plexstring, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
      strcat(plexstring, s1);
      strcat(plexstring, "NNNNNNNNNN\0");

      if (verbose)
        printf("EnergyCutoff = %f\n", pk_penalty);

      hits  = PKLduplexfold_XS(plexstring,
                               (const int **)access,
                               (int)(-pk_penalty * 100) - 1,
                               MIN2(12, length - 3),
                               0);

      /*
       ########################################################
       # analyze Plex output
       ########################################################
       */
      size_t  NumberOfHits;
      char    *mfe_struct;
      double  mfe, mfe_pk;
      vrna_fold_compound_t  *fc;

      NumberOfHits  = 0;
      mfe_struct    = (char *)vrna_alloc(sizeof(char) * (length + 1));

      fc  = vrna_fold_compound(s1, &md, VRNA_OPTION_DEFAULT);
      mfe = mfe_pk = vrna_mfe(fc, mfe_struct);

      if (hits) {
        for (hit_ptr = hits; hit_ptr->structure; hit_ptr++) {
          /* output: */
          if (verbose) {
            printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n",
                    hit_ptr->structure,
                    hit_ptr->tb,
                    hit_ptr->te,
                    hit_ptr->qb,
                    hit_ptr->qe,
                    hit_ptr->ddG,
                    hit_ptr->energy,
                    hit_ptr->dG1,
                    hit_ptr->dG2);
          }
          NumberOfHits++;
        }

        /* first sort all the pre-results descending by their estimated energy */
        qsort(hits, NumberOfHits, sizeof(vrna_pkplex_t), PlexHit_cmp);

        /*  now we re-evaluate the energies and thereby prune the list of pre-results
        *  such that the re-avaluation process is not done too often.
        */

        constraint  = (char *)vrna_alloc(sizeof(char) * (length + 1));

        for (hit_ptr = hits; hit_ptr->structure; hit_ptr++) {
          /*
          * simple check whether we believe that this structure might achieve
          * a net free energy within our boundaries
          */
          double best_e = hit_ptr->energy +
                          mfe +
                          pk_penalty +
                          MAX2(hit_ptr->dG1, hit_ptr->dG2);

          if (best_e <= mfe_pk + subopts) {
            /* now for the exact evaluation of the structures energy incl. PKs */

            /* prepare the structure constraint for constrained folding */
            vrna_hc_init(fc);

            for (i = hit_ptr->tb; i <= hit_ptr->te; i++)
              vrna_hc_add_up(fc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
            for (i = hit_ptr->qb; i <= hit_ptr->qe; i++)
              vrna_hc_add_up(fc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);

            /* energy evaluation */
            constrainedEnergy = vrna_mfe(fc, constraint);

            /* compute net free energy */
            hit_ptr->energy += constrainedEnergy +
                               pk_penalty;

            /* check if this structure is worth keeping */
            if (hit_ptr->energy <= mfe_pk + subopts) {
              /* add pseudo-knot brackets to the structure */
              for (i = hit_ptr->tb - 1; i < hit_ptr->te; i++)
                if (hit_ptr->structure[i - hit_ptr->tb + 1] == '(')
                  constraint[i] = '[';

              for (i = hit_ptr->qb - 1; i < hit_ptr->qe; i++)
                if (hit_ptr->structure[i - hit_ptr->qb + 1 + 1 + 1 +
                                                hit_ptr->te - hit_ptr->tb] == ')')
                  constraint[i] = ']';

              if (hit_ptr->energy < mfe_pk)
                mfe_pk = hit_ptr->energy;

              free(hit_ptr->structure);
              hit_ptr->structure = constraint;
              hit_ptr->processed = 1;
              hit_ptr->inactive  = 0;

              constraint = (char *)vrna_alloc(sizeof(char) * (length + 1));
              continue;
            }
          }

          hit_ptr->inactive = 1;
        }

        free(constraint);
        constraint = NULL;
        vrna_fold_compound_free(fc);
      }

      /* add the MFE structure without PKs to list */
      hits = (vrna_pkplex_t *)vrna_realloc(hits, sizeof(vrna_pkplex_t) * (NumberOfHits + 2));
      hits[NumberOfHits].structure  = mfe_struct;
      hits[NumberOfHits].energy     = mfe;
      hits[NumberOfHits].tb         = 0;
      hits[NumberOfHits++].inactive = 0;
      hits[NumberOfHits].structure  = NULL;
      hits[NumberOfHits].energy     = (double)INF * 0.01;
      hits[NumberOfHits].inactive   = 1;

      /*
      * now go through the active hits again and filter those out that are above
      * the subopt threshold
      */
      for (hit_ptr = hits; hit_ptr->structure; hit_ptr++) {
        if ((!hit_ptr->inactive) &&
            (hit_ptr->energy > mfe_pk + subopts))
          hit_ptr->inactive = 1;
      }


      /* now sort the actual results again according to their energy */

      qsort(hits, NumberOfHits, sizeof(vrna_pkplex_t), PlexHit_cmp_active_energy);

      /* and print the results to stdout */
      for (hit_ptr = hits; !hit_ptr->inactive; hit_ptr++)
        printf("%s (%6.2f)\n", hit_ptr->structure, hit_ptr->energy);

      /*
       ########################################################
       # Generate Plot for the best structure
       ########################################################
       */

      /* make an EPS layout file for the MFE structure, i.e. first in the list */
      hit_ptr = hits;
      if (hit_ptr->tb > 0) { /* true hit with PK */
        annotation = vrna_strdup_printf("%d %d 13 1 0 0 omark\n"
                                        "%d %d 13 1 0 0 omark\n"
                                        "0 0 2 setrgbcolor\n"
                                        "2 setlinewidth\n"
                                        "%d cmark\n"
                                        "%d cmark\n"
                                        "1 setlinewidth",
                                        (int)hit_ptr->tb,
                                        hit_ptr->te,
                                        (int)hit_ptr->qb,
                                        hit_ptr->qe,
                                        hit_ptr->tb,
                                        hit_ptr->qe);
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

      free(pl);
      free(pup[0]);
      free(pup);
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
