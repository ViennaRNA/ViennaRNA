/* Last changed Time-stamp: <2012-02-15 18:20:49 ivo> */
/*
                  Ineractive Access to folding Routines

                  c Ivo L Hofacker
                  Vienna RNA package
*/

/** \file
*** \brief RNAfold program source code
***
*** This code provides an interface for MFE and Partition function folding
*** of single linear or circular RNA molecules.
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
#include "read_epars.h"
#include "MEA.h"
#include "params.h"
#include "RNAfold_cmdl.h"



/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAfold.c,v 1.25 2009/02/24 14:22:21 ivo Exp $";

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[]){
  struct        RNAfold_args_info args_info;
  char          *buf, *rec_sequence, *rec_id, **rec_rest, *structure, *cstruc, *orig_sequence;
  char          fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH], *ParamFile;
  char          *ns_bases, *c;
  int           i, length, l, cl, sym, r, istty, pf, noPS, noconv, fasta;
  unsigned int  rec_type, read_opt;
  double        energy, min_en, kT, sfact;
  int           doMEA, circular, lucky;
  double        MEAgamma, bppmThreshold, betaScale;
  paramT          *mfe_parameters;
  pf_paramT       *pf_parameters;
  model_detailsT  md;

  rec_type      = read_opt = 0;
  rec_id        = buf = rec_sequence = structure = cstruc = orig_sequence = NULL;
  rec_rest      = NULL;
  ParamFile     = NULL;
  ns_bases      = NULL;
  pf_parameters = NULL;
  do_backtrack  = 1;
  pf            = 0;
  sfact         = 1.07;
  noPS          = 0;
  noconv        = 0;
  circular      = 0;
  gquad         = 0;
  fasta         = 0;
  cl            = l = length = 0;
  dangles       = 2;
  MEAgamma      = 1.;
  bppmThreshold = 1e-5;
  lucky         = 0;
  doMEA         = 0;
  betaScale     = 1.;

  set_model_details(&md);


  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* structure constraint */
  if(args_info.constraint_given)  fold_constrained=1;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     md.special_hp = tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      warn_user("required dangle model not implemented, falling back to default dangles=2");
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
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set pf scaling factor */
  if(args_info.pfScale_given)     sfact = args_info.pfScale_arg;
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)        circular=1;
  /* always look on the bright side of life */
  if(args_info.ImFeelingLucky_given)  lucky = pf = st_back = 1;
  /* set the bppm threshold for the dotplot */
  if(args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0.,args_info.bppmThreshold_arg));
  if(args_info.betaScale_given)   betaScale = args_info.betaScale_arg;
  /* do not produce postscript output */
  if(args_info.noPS_given)        noPS=1;
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != 1)
      do_backtrack = args_info.partfunc_arg;
  }
  /* MEA (maximum expected accuracy) settings */
  if(args_info.MEA_given){
    pf = doMEA = 1;
    if(args_info.MEA_arg != -1)
      MEAgamma = args_info.MEA_arg;
  }
  if(args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  /* free allocated memory of command line data structure */
  RNAfold_cmdline_parser_free (&args_info);


  /*
  #############################################
  # begin initializing
  #############################################
  */
  if(circular && gquad){
    nrerror("G-Quadruplex support is currently not available for circular RNA structures");
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (circular && noLonelyPairs)
    warn_user("depending on the origin of the circular sequence, some structures may be missed when using -noLP\nTry rotating your sequence a few times");

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

  /* print user help if we get input from tty */
  if(istty){
    if(fold_constrained){
      print_tty_constraint_full();
      print_tty_input_seq_str("Input sequence (upper or lower case) followed by structure constraint");
    }
    else print_tty_input_seq();
  }

  mfe_parameters = get_scaled_parameters(temperature, md);

  /* set options we wanna pass to read_record */
  if(istty)             read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  if(!fold_constrained) read_opt |= VRNA_INPUT_NO_REST;

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    if(rec_id){
      if(!istty) printf("%s\n", rec_id);
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    length  = (int)strlen(rec_sequence);
    structure = (char *)space(sizeof(char) *(length+1));

    /* parse the rest of the current dataset to obtain a structure constraint */
    if(fold_constrained){
      cstruc = NULL;
      unsigned int coptions = (rec_id) ? VRNA_CONSTRAINT_MULTILINE : 0;
      coptions |= VRNA_CONSTRAINT_ALL;
      getConstraint(&cstruc, (const char **)rec_rest, coptions);
      cl = (cstruc) ? (int)strlen(cstruc) : 0;

      if(cl == 0)           warn_user("structure constraint is missing");
      else if(cl < length)  warn_user("structure constraint is shorter than sequence");
      else if(cl > length)  nrerror("structure constraint is too long");
      if(cstruc) strncpy(structure, cstruc, sizeof(char)*(cl+1));
    }

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) str_DNA2RNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    str_uppercase(rec_sequence);

    if(istty) printf("length = %d\n", length);

    /*
    ########################################################
    # begin actual computations
    ########################################################
    */
    min_en = fold_par(rec_sequence, structure, mfe_parameters, fold_constrained, circular);

    if(!lucky){
      printf("%s\n%s", orig_sequence, structure);
      if (istty)
        printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
      else
        printf(" (%6.2f)\n", min_en);
      (void) fflush(stdout);

      if(fname[0] != '\0'){
        strcpy(ffname, fname);
        strcat(ffname, "_ss.ps");
      } else strcpy(ffname, "rna.ps");

      if(gquad){
        if (!noPS) (void) PS_rna_plot_a_gquad(orig_sequence, structure, ffname, NULL, NULL);
      } else {
        if (!noPS) (void) PS_rna_plot_a(orig_sequence, structure, ffname, NULL, NULL);
      }
    }
    if (length>2000) free_arrays();
    if (pf) {
      char *pf_struc = (char *) space((unsigned) length+1);
      if (md.dangles==1) {
          md.dangles=2;   /* recompute with dangles as in pf_fold() */
          min_en = (circular) ? energy_of_circ_struct_par(rec_sequence, structure, mfe_parameters, 0) : energy_of_struct_par(rec_sequence, structure, mfe_parameters, 0);
          md.dangles=1;
      }

      /* */

      kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      if (cstruc!=NULL) strncpy(pf_struc, cstruc, length+1);

      pf_parameters = get_boltzmann_factors(temperature, betaScale, md, pf_scale);
      energy = pf_fold_par(rec_sequence, pf_struc, pf_parameters, do_backtrack, fold_constrained, circular);

      if(lucky){
        init_rand();
        char *s = (circular) ? pbacktrack_circ(rec_sequence) : pbacktrack(rec_sequence);
        min_en = (circular) ? energy_of_circ_struct_par(rec_sequence, s, mfe_parameters, 0) : energy_of_struct_par(rec_sequence, s, mfe_parameters, 0);
        printf("%s\n%s", orig_sequence, s);
        if (istty)
          printf("\n free energy = %6.2f kcal/mol\n", min_en);
        else
          printf(" (%6.2f)\n", min_en);
        (void) fflush(stdout);
        if(fname[0] != '\0'){
          strcpy(ffname, fname);
          strcat(ffname, "_ss.ps");
        } else strcpy(ffname, "rna.ps");

        if (!noPS) (void) PS_rna_plot(orig_sequence, s, ffname);
        free(s);
      }
      else{
      
        if (do_backtrack) {
          printf("%s", pf_struc);
          if (!istty) printf(" [%6.2f]\n", energy);
          else printf("\n");
        }
        if ((istty)||(!do_backtrack))
          printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);


        if (do_backtrack) {
          plist *pl1,*pl2;
          char *cent;
          double dist, cent_en;
          FLT_OR_DBL *probs = export_bppm();

          if(gquad)
            assign_plist_gquad_from_pr(&pl1, length, bppmThreshold);
          else
            assign_plist_from_pr(&pl1, probs, length, bppmThreshold);

          assign_plist_from_db(&pl2, structure, 0.95*0.95);
          /* cent = centroid(length, &dist); <- NOT THREADSAFE */

          if(gquad){
            cent    = get_centroid_struct_gquad_pr(length, &dist);
            cent_en = energy_of_gquad_structure((const char *)rec_sequence, (const char *)cent, 0);
          } else {
            cent    = get_centroid_struct_pr(length, &dist, probs);
            cent_en = (circular) ? energy_of_circ_struct_par(rec_sequence, cent, mfe_parameters, 0) : energy_of_struct_par(rec_sequence, cent, mfe_parameters, 0);
          }

          printf("%s {%6.2f d=%.2f}\n", cent, cent_en, dist);
          free(cent);
          if (fname[0]!='\0') {
            strcpy(ffname, fname);
            strcat(ffname, "_dp.ps");
          } else strcpy(ffname, "dot.ps");
          (void) PS_dot_plot_list(orig_sequence, ffname, pl1, pl2, "");
          free(pl2);
          if (do_backtrack==2) {
            pl2 = stackProb(1e-5);
            if (fname[0]!='\0') {
              strcpy(ffname, fname);
              strcat(ffname, "_dp2.ps");
            } else strcpy(ffname, "dot2.ps");
            PS_dot_plot_list(orig_sequence, ffname, pl1, pl2,
                             "Probabilities for stacked pairs (i,j)(i+1,j-1)");
            free(pl2);
          }
          free(pl1);
          free(pf_struc);
          if(doMEA){
            float mea, mea_en;
            plist *pl;
            assign_plist_from_pr(&pl, probs, length, 1e-4/(1+MEAgamma));

            if(gquad){
              mea = MEA_seq(pl, rec_sequence, structure, MEAgamma, pf_parameters);
              mea_en = energy_of_gquad_structure((const char *)rec_sequence, (const char *)structure, 0);
              printf("%s {%6.2f MEA=%.2f}\n", structure, mea_en, mea);
            } else {
              mea = MEA(pl, structure, MEAgamma);
              mea_en = (circular) ? energy_of_circ_struct_par(rec_sequence, structure, mfe_parameters, 0) : energy_of_struct_par(rec_sequence, structure, mfe_parameters, 0);
              printf("%s {%6.2f MEA=%.2f}\n", structure, mea_en, mea);
            }

            free(pl);
          }
        }
        printf(" frequency of mfe structure in ensemble %g; ", exp((energy-min_en)/kT));
        if (do_backtrack)
          printf("ensemble diversity %-6.2f", mean_bp_distance(length));
        printf("\n");
      }
      free_pf_arrays();
      free(pf_parameters);
    }
    (void) fflush(stdout);

    /* clean up */
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
    rec_id = rec_sequence = structure = cstruc = NULL;
    rec_rest = NULL;

    /* print user help for the next round if we get input from tty */
    if(istty){
      if(fold_constrained){
        print_tty_constraint_full();
        print_tty_input_seq_str("Input sequence (upper or lower case) followed by structure constraint");
      }
      else print_tty_input_seq();
    }
  }
  
  free(mfe_parameters);
  return EXIT_SUCCESS;
}
