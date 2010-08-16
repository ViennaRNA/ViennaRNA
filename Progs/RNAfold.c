/* Last changed Time-stamp: <2007-12-05 13:55:42 ronny> */
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
#include "RNAfold_cmdl.h"

/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAfold.c,v 1.25 2009/02/24 14:22:21 ivo Exp $";

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[]){
  struct        RNAfold_args_info args_info;
  char          *string, *input_string, *structure=NULL, *cstruc=NULL;
  char          fname[80], ffname[80], gfname[80], *ParamFile=NULL;
  char          *ns_bases=NULL, *c;
  int           i, length, l, sym, r, istty, pf=0, noPS=0, noconv=0;
  unsigned int  input_type;
  double        energy, min_en, kT, sfact=1.07;
  int           doMEA=0, circular = 0;
  double        MEAgamma = 1.;

  do_backtrack  = 1;
  string        = NULL;

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
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)     dangles = args_info.dangles_arg;
  /* do not allow weak pairs */
  if(args_info.noLP_given)        noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) no_closingGU = 1;
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
  /* do not produce postscript output */
  if(args_info.noPS_given)        noPS=1;
  /* partition function settings */
  if(args_info.partfunc_given){
    pf = 1;
    if(args_info.partfunc_arg != -1)
      do_backtrack = args_info.partfunc_arg;
  }
  /* MEA (maximum expected accuracy) settings */
  if(args_info.MEA_given){
    pf = doMEA = 1;
    if(args_info.MEA_arg != -1)
      MEAgamma = args_info.MEA_arg;
  }

  /* free allocated memory of command line data structure */
  RNAfold_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
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

  if(fold_constrained && istty) print_tty_constraint_full();

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do{
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty) print_tty_input_seq();

    /* extract filename from fasta header if available */
    fname[0] = '\0';
    while((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER){
      (void) sscanf(input_string, "%42s", fname);
      printf(">%s\n", input_string); /* print fasta header if available */
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      length = (int)    strlen(input_string);
      string = strdup(input_string);
      free(input_string);
    }

    structure = (char *) space((unsigned) length+1);

    if(noconv)  str_RNA2RNA(string);
    else        str_DNA2RNA(string);

    if(istty) printf("length = %d\n", length);

    /* get structure constraint or break if necessary, entering an empty line results in a warning */
    if (fold_constrained) {
      input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
      if(input_type & VRNA_INPUT_QUIT){ break;}
      else if((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)){
        cstruc = strdup(input_string);
        free(input_string);
        strncpy(structure, cstruc, length);
      }
      else warn_user("constraints missing");
    }
    /*
    ########################################################
    # done with 'stdin' handling
    ########################################################
    */

    min_en = (circular) ? circfold(string, structure) : fold(string, structure);

    printf("%s\n%s", string, structure);
    if (istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);

    (void) fflush(stdout);

    if (fname[0]!='\0') {
      strcpy(ffname, fname);
      strcat(ffname, "_ss.ps");
      strcpy(gfname, fname);
      strcat(gfname, "_ss.g");
    } else {
      strcpy(ffname, "rna.ps");
      strcpy(gfname, "rna.g");
    }
    if (!noPS) (void) PS_rna_plot(string, structure, ffname);
    if (length>2000) free_arrays();
    if (pf) {
      char *pf_struc;
      pf_struc = (char *) space((unsigned) length+1);
        if (dangles==1) {
          dangles=2;   /* recompute with dangles as in pf_fold() */
          min_en = (circular) ? energy_of_circ_struct(string, structure) : energy_of_struct(string, structure);
          dangles=1;
      }

      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      if (cstruc!=NULL) strncpy(pf_struc, cstruc, length+1);

      energy = (circular) ? pf_circ_fold(string, pf_struc) : pf_fold(string, pf_struc);

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
        assign_plist_from_pr(&pl1, pr, length, 1e-5);
        assign_plist_from_db(&pl2, structure, 0.95*0.95);
        /* cent = centroid(length, &dist); <- NOT THREADSAFE */
        cent = get_centroid_struct_pl(length, &dist, pl1);
        cent_en = (circular) ? energy_of_circ_struct(string, cent) :energy_of_struct(string, cent);
        printf("%s {%6.2f d=%.2f}\n", cent, cent_en, dist);
        free(cent);
        if (fname[0]!='\0') {
          strcpy(ffname, fname);
          strcat(ffname, "_dp.ps");
        } else strcpy(ffname, "dot.ps");
        (void) PS_dot_plot_list(string, ffname, pl1, pl2, "");
        free(pl2);
        if (do_backtrack==2) {
          pl2 = stackProb(1e-5);
          if (fname[0]!='\0') {
            strcpy(ffname, fname);
            strcat(ffname, "_dp2.ps");
          } else strcpy(ffname, "dot2.ps");
          PS_dot_plot_list(string, ffname, pl1, pl2,
                           "Probabilities for stacked pairs (i,j)(i+1,j-1)");
          free(pl2);
        }
        free(pl1);
        free(pf_struc);
        if(doMEA){
          float mea, mea_en;
          plist *pl;
          assign_plist_from_pr(&pl, pr, length, 1e-4/(1+MEAgamma));
          mea = MEA(pl, structure, MEAgamma);
          mea_en = (circular) ? energy_of_circ_struct(string, structure) : energy_of_struct(string, structure);
          printf("%s {%6.2f MEA=%.2f}\n", structure, mea_en, mea);
          free(pl);
        }
      }
      printf(" frequency of mfe structure in ensemble %g; ", exp((energy-min_en)/kT));
      if (do_backtrack)
        printf("ensemble diversity %-6.2f", mean_bp_dist(length));
      printf("\n");

      free_pf_arrays();
    }
    if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    free(string);
    free(structure);
  } while (1);
  return 0;
}
