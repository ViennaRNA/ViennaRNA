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
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_const.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/PKplex.h"
#include "RNAPKplex_cmdl.h"

int PlexHit_cmp(const void  *c1,
                const void  *c2);


int PlexHit_cmp_energy(const void *c1,
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
  int                     istty, i, j, noconv, length, pairdist, current, unpaired, winsize;
  float                   cutoff, constrainedEnergy;
  double                  **pup, subopts, pk_penalty;
  plist                   *pl, *dpp;
  vrna_md_t               md;
  vrna_param_t            *par;

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
  par         = NULL;
  annotation  = NULL;

  set_model_details(&md);

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
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

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
  while (!(vrna_file_fasta_read_record(&id_s1, &s1, &rest, NULL, options) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
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
    update_fold_params();
    if (length >= 5) {
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
      NumberOfHits  = 0;
      PlexHits      = (dupVar *)vrna_alloc(sizeof(dupVar) * PlexHitsArrayLength);
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
          if (pup[i - 10][j - 1 + 1] > 0)
            access[j][i] = rint(100 * (-log(pup[i - 10][j - 1 + 1])) * kT);
      }

      access[0][0] = unpaired + 2;

      plexstring = (char *)vrna_alloc(length + 1 + 20);
      strcpy(plexstring, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
      strcat(plexstring, s1);
      strcat(plexstring, "NNNNNNNNNN\0");

      if (verbose)
        printf("EnergyCutoff = %f\n", pk_penalty);

      PKLduplexfold_XS(plexstring, access, (int)(-pk_penalty * 100) - 1, MIN2(12, length - 3), 0);

      /*
       ########################################################
       # analyze Plex output
       ########################################################
       */

      /*adding empty hit*/
      PlexHits[NumberOfHits].tb         = 0;
      PlexHits[NumberOfHits].te         = 0;
      PlexHits[NumberOfHits].qb         = 0;
      PlexHits[NumberOfHits].qe         = 0;
      PlexHits[NumberOfHits].ddG        = 0;
      PlexHits[NumberOfHits].dG1        = 0;
      PlexHits[NumberOfHits].dG2        = 0;
      PlexHits[NumberOfHits].energy     = 0;
      PlexHits[NumberOfHits].structure  = NULL;
      NumberOfHits++;

      /* first sort all the pre-results descending by their estimated energy */
      qsort(PlexHits, NumberOfHits, sizeof(dupVar), PlexHit_cmp);

      /*  now we re-evaluate the energies and thereby prune the list of pre-results
       *  such that the re-avaluation process is not done too often.
       */
      double  mfe         = 0.;
      double  mfe_pk      = 0.;
      char    *mfe_struct = NULL;

      par         = vrna_params(&md);
      constraint  = (char *)vrna_alloc(length + 1);
      mfe_struct  = (char *)vrna_alloc(length + 1);

      mfe = mfe_pk = fold_par(s1, mfe_struct, par, 0, 0);
      /*      if(verbose)
       *        printf("%s (%6.2f) [mfe-pkfree]\n", mfe_struct, mfe);
       */
      for (current = 0; current < NumberOfHits; current++) {
        /* do evaluation for structures above the subopt threshold only */
        if (!PlexHits[current].inactive) {
          if (PlexHits[current].structure) {
            /* prepare the structure constraint for constrained folding */
            for (i = 0; i < length; i++)
              constraint[i] = '.';
            for (i = PlexHits[current].tb - 1; i < PlexHits[current].te; i++)
              constraint[i] = 'x';
            for (i = PlexHits[current].qb - 1; i < PlexHits[current].qe; i++)
              constraint[i] = 'x';
            constraint[length] = '\0';

            /* energy evaluation */
            constrainedEnergy = fold_par(s1, constraint, par, 1, 0);

            /* check if this structure is worth keeping */
            if (constrainedEnergy + PlexHits[current].ddG + pk_penalty <= mfe_pk + subopts) {
              /* add pseudo-knot brackets to the structure */
              for (i = PlexHits[current].tb - 1; i < PlexHits[current].te; i++)
                if (PlexHits[current].structure[i - PlexHits[current].tb + 1] == '(')
                  constraint[i] = '[';

              for (i = PlexHits[current].qb - 1; i < PlexHits[current].qe; i++)
                if (PlexHits[current].structure[i - PlexHits[current].qb + 1 + 1 + 1 + PlexHits[current].te - PlexHits[current].tb] == ')')
                  constraint[i] = ']';

              PlexHits[current].energy = constrainedEnergy + PlexHits[current].ddG + (float)pk_penalty + PlexHits[current].dG1 + PlexHits[current].dG2;
              if (PlexHits[current].energy < mfe_pk)
                mfe_pk = PlexHits[current].energy;

              free(PlexHits[current].structure);
              PlexHits[current].structure = strdup(constraint);
              PlexHits[current].processed = 1;
            } else {
              PlexHits[current].inactive = 1;
            }
          } else {
            PlexHits[current].energy = mfe;
            if (mfe > mfe_pk + subopts)
              PlexHits[current].inactive = 1;
          }

          /*
           * now go through the rest of the PlexHits array and mark all hits as inactive if they can
           * definetely not be within the results set according to current subopt settings
           */
          for (i = 0; i < NumberOfHits; i++) {
            if (!PlexHits[i].inactive) {
              if (PlexHits[i].structure) {
                if (!PlexHits[i].processed) {
                  double cost = PlexHits[i].dG1;
                  if (PlexHits[i].dG2 > cost)
                    cost = PlexHits[i].dG2;

                  cost += mfe + PlexHits[i].ddG + pk_penalty;
                  if (cost > mfe + subopts)
                    PlexHits[i].inactive = 1;
                } else {
                  if (PlexHits[i].energy > mfe_pk + subopts)
                    PlexHits[i].inactive = 1;
                }
              } else {
                if (mfe > mfe_pk + subopts)
                  PlexHits[i].inactive = 1;
              }
            }
          }
        }
      }
      constraint = NULL;
      free(par);

      /* now sort the actual results again according to their energy */

      qsort(PlexHits, NumberOfHits, sizeof(dupVar), PlexHit_cmp_energy);

      /* and print the results to stdout */
      for (i = 0; i < NumberOfHits; i++) {
        if (!PlexHits[i].inactive) {
          if (PlexHits[i].structure)
            printf("%s (%6.2f)\n", PlexHits[i].structure, PlexHits[i].energy);
          else
            printf("%s (%6.2f)\n", mfe_struct, mfe);
        }
      }

      /*
       ########################################################
       # Generate Plot for the best structure
       ########################################################
       */

      /* and print the results to stdout */
      for (i = 0; i < NumberOfHits; i++) {
        if (!PlexHits[i].inactive) {
          if (PlexHits[i].structure) {
            annotation = vrna_strdup_printf("%d %d 13 1 0 0 omark\n"
                                            "%d %d 13 1 0 0 omark\n"
                                            "0 0 2 setrgbcolor\n"
                                            "2 setlinewidth\n"
                                            "%d cmark\n"
                                            "%d cmark\n"
                                            "1 setlinewidth",
                                            (int)PlexHits[i].tb,
                                            PlexHits[i].te,
                                            (int)PlexHits[i].qb,
                                            PlexHits[i].qe,
                                            PlexHits[i].tb,
                                            PlexHits[i].qe);
            vrna_file_PS_rnaplot_a(s1, PlexHits[i].structure, fname, annotation, "", &md);
            free(annotation);
            annotation = NULL;
          } else {
            vrna_file_PS_rnaplot(s1, mfe_struct, fname, &md);
          }

          break;
        }
      }

#if 0
      if (verbose) {
        printf("\n");
        for (i = 0; i < NumberOfHits; i++)
          printf("%d\t%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", PlexHits[i].inactive, PlexHits[i].structure, PlexHits[i].tb, PlexHits[i].te, PlexHits[i].qb, PlexHits[i].qe, PlexHits[i].ddG, PlexHits[i].energy, PlexHits[i].dG1, PlexHits[i].dG2);
      }

      current = -1;

      while ((PlexHits[0].ddG + subopts >= PlexHits[current + 1].ddG) && (current + 1 < NumberOfHits)) {
        current++;

        /*
         ########################################################
         # Constrained RNAfold to reevaluate the actual energy
         ########################################################
         */
        constraint = (char *)vrna_alloc(length + 1);
        for (i = 0; i < length; i++) {
          if ((PlexHits[current].tb - 1 <= i) && (PlexHits[current].te - 1 >= i))
            constraint[i] = 'x';
          else if ((PlexHits[current].qb - 1 <= i) && (PlexHits[current].qe - 1 >= i))
            constraint[i] = 'x';
          else
            constraint[i] = '.';
        }
        constraint[length] = '\0';
        if (verbose)
          printf("Constrained structure:\n%s\n%s\n", orig_s1, constraint);

        fold_constrained  = 1;
        constrainedEnergy = fold(s1, constraint);
        if (verbose)
          printf("%s   %f\n", constraint, constrainedEnergy);

        /*
         ########################################################
         # Fusion Structure
         ########################################################
         */
        if (PlexHits[current].structure) {
          for (i = PlexHits[current].tb - 1; i <= PlexHits[current].te - 1; i++)
            if (PlexHits[current].structure[i - PlexHits[current].tb + 1] == '(')
              constraint[i] = '[';

          for (i = PlexHits[current].qb - 1; i <= PlexHits[current].qe - 1; i++)
            if (PlexHits[current].structure[i - PlexHits[current].qb + 1 + 1 + 1 + PlexHits[current].te - PlexHits[current].tb] == ')')
              constraint[i] = ']';
        }

        if (PlexHits[current].structure)
          printf("%s   (%3.2f)\n", constraint, constrainedEnergy + PlexHits[current].ddG - (float)pk_penalty);
        else
          printf("%s   (%3.2f)\n", constraint, constrainedEnergy + PlexHits[current].ddG);

        if (current == 0) {
          /*
           ########################################################
           # Generate Visualization
           ########################################################
           */

          annotation  = (char *)vrna_alloc(sizeof(char) * 300);
          temp        = (char *)vrna_alloc(sizeof(char) * 300);

          if (PlexHits[current].te) {
            int start = 0;
            int end;
            int stem = 1;
            for (i = 1; PlexHits[current].structure[i] != ')'; i++) {
              if ((stem) && (PlexHits[current].structure[i] != '(')) {
                end   = i - 1;
                stem  = 0;
                sprintf(temp, "%d %d 13 1 0 0 omark\n", (int)PlexHits[current].tb + start, PlexHits[current].tb + end);
                strcat(annotation, temp);
              }

              if ((!stem) && (PlexHits[current].structure[i] == '(')) {
                start = i;
                stem  = 1;
              }
            }
            stem  = 1;
            start = i;
            for (i; i <= strlen(PlexHits[current].structure); i++) {
              if ((stem) && (PlexHits[current].structure[i] != ')')) {
                end   = i - 1;
                stem  = 0;
                sprintf(temp, "%d %d 13 1 0 0 omark\n", PlexHits[current].qb + start - PlexHits[current].te + PlexHits[current].tb - 2, PlexHits[current].qb + end - PlexHits[current].te + PlexHits[current].tb - 2);
                strcat(annotation, temp);
              }

              if ((!stem) && (PlexHits[current].structure[i] == ')')) {
                start = i;
                stem  = 1;
              }
            }

            sprintf(temp, "0 0 2 setrgbcolor\n2 setlinewidth\n%d cmark\n%d cmark\n1 setlinewidth", PlexHits[current].tb, PlexHits[current].qe);
            strcat(annotation, temp);
            vrna_file_PS_rnaplot_a(s1, constraint, fname, annotation, "", &md);
            free(annotation);
            free(temp);
          } else {
            vrna_file_PS_rnaplot(s1, constraint, fname, &md);
          }
        }
      }
#endif

      /*
       ########################################################
       # free memory
       ########################################################
       */
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
    free(PlexHits);
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
