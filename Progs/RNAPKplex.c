/* Last changed Time-stamp: <2010-06-30 17:42:12 wolfgang> */
/*
             Compute pseudoknotted structure of an RNA

                           c Ivo L Hofacker
                          Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "fold_vars.h"
#include "utils.h"
#include "energy_const.h"
#include "LPfold.h"
#include "RNAPKplex_cmdl.h"
#include "PS_dot.h"
#include "fold.h"
#include "read_epars.h"
#include "PKplex.h"

int PlexHit_cmp (const void *c1, const void *c2);
short *make_pk_pair_table(const char *structure);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
  struct        PKplex_args_info args_info;
  char          *id_s1, *s1, *orig_s1, *ParamFile, *ns_bases, *c, *plexstring, *constraint;
  char          fname[FILENAME_MAX_LENGTH], *temp, *annotation, **rest;
  int           istty, l, i, j, noconv, length, pairdist, current, unpaired;
  int           winsize, openenergies, sym, energyCutoff;
  double        **pup = NULL; /*prob of being unpaired, lengthwise*/
  FILE          *pUfp = NULL, *spup = NULL;
  plist         *pl, *dpp = NULL;
  float         cutoff, constrainedEnergy;
  double        subopts=0;
  unsigned int  options=0;

  dangles      = 2;
  winsize      = 70;
  cutoff       = 0.01;
  pairdist     = 0;
  unpaired     = 0;
  noconv       = 0;
  openenergies = 1;
  energyCutoff = -810;
  ParamFile = ns_bases = NULL;
  s1 = id_s1 = orig_s1 = NULL;

  /*
  #############################################
  # check command line parameters
  #############################################
  */
  if(PKplex_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)              temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)           tetra_loop=0;
  /* do not allow weak pairs */
  if(args_info.noLP_given)              noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)              noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given)       no_closingGU = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)            noconv = 1;
  /* take another energy parameter set */
  if(args_info.paramFile_given)         ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)               ns_bases = strdup(args_info.nsp_arg);
  /* set the pair probability cutoff */
  if(args_info.cutoff_given)            cutoff = args_info.cutoff_arg;
  /* turn on verbose output (mainly for debugging) */
  if(args_info.verbose_given)           verbose = 1;
  /* set energy cutoff */
  if(args_info.energyCutoff_given)      energyCutoff = args_info.energyCutoff_arg;
  /* show suboptimal structures which are better than given value difference */
  if(args_info.subopts_given)           subopts = (double) args_info.subopts_arg/100;
  /* free allocated memory of command line data structure */
  PKplex_cmdline_parser_free(&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  if (ParamFile != NULL) {
    read_parameter_file(ParamFile);
  }

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
  if(istty) options |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  options |= VRNA_INPUT_NO_REST;
  if(istty) print_tty_input_seq();

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(!(read_record(&id_s1, &s1, &rest, options) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(id_s1){
      printf("%s\n", id_s1);
      (void) sscanf(id_s1, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else {
      strcpy(fname, "PKplex");
    }
    strcat(fname, ".ps");

    length=strlen(s1);
    winsize=pairdist=length;
    unpaired=MIN2(30, length-3);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) str_DNA2RNA(s1);
    /* store case-unmodified sequence */
    orig_s1 = strdup(s1);
    /* convert sequence to uppercase letters only */
    str_uppercase(s1);

    printf("%s\n", orig_s1);
    if (verbose) printf("length = %d\n", length);
    /*
    ########################################################
    # do PLfold computations
    ########################################################
    */
    update_fold_params();
    if (length >= 5){

      pf_scale  = -1;

      pup       =(double **)  space((length+1)*sizeof(double *));
      pup[0]    =(double *)   space(sizeof(double)); /*I only need entry 0*/
      pup[0][0] = unpaired;

      pUfp = spup = NULL;

      if (verbose) printf ("Winsize = %d\nPairdist = %d\nUnpaired = %d\n", winsize, pairdist, unpaired);
      int tempdangles = dangles;
      dangles = 2;
      pl = pfl_fold(s1, winsize, pairdist, cutoff, pup, &dpp, pUfp, spup);
      dangles = tempdangles;

    /*
    ########################################################
    # do Plex computations
    ########################################################
    */
      NumberOfHits=0;
      PlexHits = (dupVar *)space(sizeof(dupVar) * PlexHitsArrayLength);
      double kT= (temperature+K0)*GASCONST/1000.0;
      int **access;
      access = (int**) space(sizeof(int *) * (unpaired+2));
      for(i=0; i< unpaired+2; i++){
        access[i] =(int *) space(sizeof(int) * (length+20));
      }

      for(i=0;i<length+20;i++){
        for(j=0;j<unpaired+2;j++){
          access[j][i]=INF;
        }
      }

      for(i=11;i<length+11;i++){
        for(j=1;j<unpaired+1;j++){
          if (pup[i-10][j-1+1]>0) {
            access[j][i]=rint(100*(-log(pup[i-10][j-1+1]))*kT);
          }
        }
      }

      access[0][0]=unpaired+2;

      plexstring = (char *) space(length+1+20);
      strcpy(plexstring,"NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
      strcat(plexstring, s1);
      strcat(plexstring,"NNNNNNNNNN\0");

      if (verbose) printf("EnergyCutoff = %d\n", energyCutoff);
      PKLduplexfold_XS(plexstring, access, energyCutoff+1, MIN2(12, length-3), 0);

    /*
    ########################################################
    # analyze Plex output
    ########################################################
    */

      /*adding empty hit*/
      PlexHits[NumberOfHits].tb=0;
      PlexHits[NumberOfHits].te=0;
      PlexHits[NumberOfHits].qb=0;
      PlexHits[NumberOfHits].qe=0;
      PlexHits[NumberOfHits].ddG=0;
      PlexHits[NumberOfHits].dG1=0;
      PlexHits[NumberOfHits].dG2=0;
      PlexHits[NumberOfHits].energy=0;
      PlexHits[NumberOfHits].structure='\0';
      NumberOfHits++;

      qsort(PlexHits, NumberOfHits, sizeof(dupVar), PlexHit_cmp);

      if (verbose) {
        printf("\n");
        for(i=0;i<NumberOfHits;i++) {
          printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", PlexHits[i].structure, PlexHits[i].tb, PlexHits[i].te, PlexHits[i].qb, PlexHits[i].qe, PlexHits[i].ddG, PlexHits[i].energy, PlexHits[i].dG1, PlexHits[i].dG2);
        }
      }

      current=-1;
      while((PlexHits[0].ddG+subopts>=PlexHits[current+1].ddG) && (current+1<NumberOfHits)) {
        current++;

        /*
        ########################################################
        # Constrained RNAfold
        ########################################################
        */
        constraint = (char *) space(length+1);
        for(i=0; i<length; i++) {
          if((PlexHits[current].tb-1<=i) && (PlexHits[current].te-1>=i)) {
            constraint[i]='x';
          } else if((PlexHits[current].qb-1<=i) && (PlexHits[current].qe-1>=i)) {
            constraint[i]='x';
          } else {
            constraint[i]='.';
          }
        }
        constraint[length]='\0';
        if (verbose) printf("Constrained structure:\n%s\n%s\n", orig_s1, constraint);

        fold_constrained=1;
        constrainedEnergy=fold(s1, constraint);
        if (verbose) printf("%s   %f\n", constraint, constrainedEnergy);

        /*
        ########################################################
        # Fusion Structure
        ########################################################
        */
        if (PlexHits[current].structure) {
          for(i=PlexHits[current].tb-1; i<=PlexHits[current].te-1; i++) {
            if(PlexHits[current].structure[i-PlexHits[current].tb+1]=='(') {
              constraint[i]='[';
            }
          }
          for(i=PlexHits[current].qb-1; i<=PlexHits[current].qe-1; i++) {
            if(PlexHits[current].structure[i-PlexHits[current].qb+1+1+1+PlexHits[current].te-PlexHits[current].tb]==')') {
              constraint[i]=']';
            }
          }
        }

        if (PlexHits[current].structure) {
          printf("%s   (%3.2f)\n", constraint, constrainedEnergy+PlexHits[current].ddG-(float) energyCutoff/100);
        } else {
          printf("%s   (%3.2f)\n", constraint, constrainedEnergy+PlexHits[current].ddG);
        }

        if(current==0) {
        /*
        ########################################################
        # Generate Visualization
        ########################################################
        */

          annotation = (char *) space(sizeof(char)*300);
          temp = (char *) space(sizeof(char)*300);

          if (PlexHits[current].te) {
            int start=0;
            int end;
            int stem=1;
            for (i=1; PlexHits[current].structure[i]!=')'; i++) {
              if ((stem) && (PlexHits[current].structure[i]!='(')) {
                end=i-1;
                stem=0;
                sprintf(temp, "%d %d 13 1 0 0 omark\n", (int) PlexHits[current].tb+start, PlexHits[current].tb+end);
                strcat(annotation, temp);
              }
              if ((!stem) && (PlexHits[current].structure[i]=='(')) {
                start=i;
                stem=1;
              }
            }
            stem=1;
            start=i;
            for (i; i<=strlen(PlexHits[current].structure); i++) {
              if ((stem) && (PlexHits[current].structure[i]!=')')) {
                end=i-1;
                stem=0;
                sprintf(temp, "%d %d 13 1 0 0 omark\n", PlexHits[current].qb+start-PlexHits[current].te+PlexHits[current].tb-2, PlexHits[current].qb+end-PlexHits[current].te+PlexHits[current].tb-2);
                strcat(annotation, temp);
              }
              if ((!stem) && (PlexHits[current].structure[i]==')')) {
                start=i;
                stem=1;
              }
            }

            sprintf(temp, "0 0 2 setrgbcolor\n2 setlinewidth\n%d cmark\n%d cmark\n1 setlinewidth", PlexHits[current].tb, PlexHits[current].qe);
            strcat(annotation, temp);
            PS_rna_plot_a(s1, constraint, fname, annotation, "");
            free(annotation);
            free(temp);
          } else {
            PS_rna_plot(s1, constraint, fname);
          }
        }
      }

    /*
    ########################################################
    # free memory
    ########################################################
    */
      free(pl);
      free(pup[0]);
      free(pup);
      (void) fflush(stdout);
      i =  access[0][0];
      while(--i>-1){
        free(access[i]);
      }
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

/* print user help for the next round if we get input from tty */
    if(istty) print_tty_input_seq();
  }
  return 0;
}

int PlexHit_cmp (const void *c1, const void *c2) {
  dupVar *p1=(dupVar *)c1;
  dupVar *p2=(dupVar *)c2;
  return (p1->ddG >= p2->ddG);
}

short *make_pk_pair_table(const char *structure) {
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
  short i, j, hrund, heckig;
  short length;
  short *stackrund, *stackeckig;
  short *table;

  length = (short) strlen(structure);
  stackrund = (short *) space(sizeof(short)*(length+1));
  stackeckig = (short *) space(sizeof(short)*(length+1));
  table = (short *) space(sizeof(short)*(length+2));
  table[0] = length;

  for (hrund=0, heckig=0, i=1; i<=length; i++) {
    switch (structure[i-1]) {
      case '(':
        stackrund[hrund++]=i;
        break;
      case ')':
        j = stackrund[--hrund];
        if (hrund<0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced () brackets in make_pk_pair_table");
        }
        table[i]=j;
        table[j]=i;
        break;
      case '[':
        stackeckig[heckig++]=i;
        break;
      case ']':
        j = stackeckig[--heckig];
        if (heckig<0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced [] brackets in make_pk_pair_table");
        }
        table[i]=j;
        table[j]=i;
        break;
      default:   /* unpaired base, usually '.' */
        table[i]= 0;
      break;
    }
  }
  if (hrund!=0) {
    fprintf(stderr, "%s\n", structure);
    nrerror("unbalanced () brackets in make_pk_pair_table");
  }
  if (heckig!=0) {
    fprintf(stderr, "%s\n", structure);
    nrerror("unbalanced [] brackets in make_pk_pair_table");
  }

  free(stackrund);
  free(stackeckig);
  return(table);
}
