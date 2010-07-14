/* Last changed Time-stamp: <2008-12-03 16:38:01 ivo> */
/*
                Ineractive access to suboptimal folding

                           c Ivo L Hofacker
                          Vienna RNA package
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "part_func.h"
#include "fold.h"
#include "cofold.h"
#include "fold_vars.h"
#include "utils.h"
#include "read_epars.h"
#include "subopt.h"
#include "RNAsubopt_cmdl.h"

/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAsubopt.c,v 1.20 2008/12/03 16:55:44 ivo Exp $";

PRIVATE char *tokenize(char *line);
PRIVATE void putoutzuker(SOLUTION* zukersolution);

int main(int argc, char *argv[]){
  struct        RNAsubopt_args_info args_info;
  unsigned int  input_type;
  char          fname[80], *cstruc, *sequence, *c, *input_string;
  char          *structure = NULL, *ParamFile = NULL, *ns_bases = NULL;
  int           i, length, l, sym, istty;
  double        deltaf, deltap;
  int           delta, n_back, noconv, circular, dos, zuker;

  do_backtrack  = 1;
  dangles       = 2;
  delta         = 100;
  deltap = n_back = noconv = circular = dos = zuker = 0;
  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAsubopt_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
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
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* energy range */
  if(args_info.deltaEnergy_given) delta = (int) (0.1+args_info.deltaEnergy_arg*100);
  /* energy range after post evaluation */
  if(args_info.deltaEnergyPost_given) deltap = args_info.deltaEnergyPost_arg;
  /* sorted output */
  if(args_info.sorted_given)      subopt_sorted = 1;
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)        circular=1;
  /* stochastic backtracking */
  if(args_info.stochBT_given){
    n_back = args_info.stochBT_arg;
    init_rand();
  }
  /* density of states */
  if(args_info.dos_given){
    dos = 1;
    print_energy = -999999;
  }
  /* logarithmic multiloop energies */
  if(args_info.logML_given) logML = 1;
  /* zuker subopts */
  if(args_info.zuker_given) zuker = 1;

  if(zuker){
    if(circular){
      warn_user("Sorry, zuker subopts not yet implemented for circfold");
      RNAsubopt_cmdline_parser_print_help();
      exit(1);
    }
    else if(n_back>0){
      warn_user("Can't do zuker subopts and stochastic subopts at the same time");
      RNAsubopt_cmdline_parser_print_help();
      exit(1);
    }
  }

  /* free allocated memory of command line data structure */
  RNAsubopt_cmdline_parser_free(&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */

  if (ParamFile != NULL) read_parameter_file(ParamFile);

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

  if(fold_constrained && istty) print_tty_constraint(VRNA_CONSTRAINT_DOT | VRNA_CONSTRAINT_X);


  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do {
    cut_point = -1;
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty){ 
      if (!zuker)
        printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      print_tty_input_seq();
    }
    /* extract filename from fasta header if available */
    fname[0] = '\0';
    while((input_type = get_input_line(&input_string, 0)) == VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string);
      (void) sscanf(input_string, "%42s", fname);
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      sequence  = tokenize(input_string); /* frees input_string */
      length    = (int) strlen(sequence);
    }
    structure = (char *) space((unsigned) length+1);

    if(noconv)  str_RNA2RNA(sequence);
    else        str_DNA2RNA(sequence);

    if(istty){
      if (cut_point == -1)
        printf("length = %d\n", length);
      else
        printf("length1 = %d\nlength2 = %d\n", cut_point-1, length-cut_point+1);
    }

    /* get structure constraint or break if necessary, entering an empty line results in a warning */
    if (fold_constrained) {
      input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
      if(input_type & VRNA_INPUT_QUIT){ break;}
      else if((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)){
        cstruc = tokenize(input_string);
        strncpy(structure, cstruc, length);
        for (i=0; i<length; i++)
          if (structure[i]=='|')
            nrerror("constraints of type '|' not allowed");
        free(cstruc);
      }
      else warn_user("constraints missing");
    }
    /*
    ########################################################
    # done with 'stdin' handling, now init everything properly
    ########################################################
    */

    if((logML != 0 || dangles==1 || dangles==3) && dos == 0)
      if(deltap<=0) deltap = delta/100. + 0.001;
    if (deltap>0)
      print_energy = deltap;

    /* first lines of output (suitable  for sort +1n) */
    if (fname[0] != '\0')
      printf("> %s [%d]\n", fname, delta);

    /* stochastic backtracking */
    if(n_back>0){
      double mfe, kT;
      char *ss;
      st_back=1;
      ss = (char *) space(strlen(sequence)+1);
      strncpy(ss, structure, length);
      mfe = fold(sequence, ss);
      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      pf_scale = exp(-(1.03*mfe)/kT/length);
      strncpy(ss, structure, length);
      /* ignore return value, we are not interested in the free energy */
      (circular) ? (void) pf_circ_fold(sequence, ss) : (void) pf_fold(sequence, ss);
      free(ss);
      for (i=0; i<n_back; i++) {
        char *s;
        s =(circular) ? pbacktrack_circ(sequence) : pbacktrack(sequence);
        printf("%s\n", s);
        free(s);
      }
      free_pf_arrays();
    }
    /* normal subopt */
    else if(!zuker){
      (circular) ? subopt_circ(sequence, structure, delta, stdout) : subopt(sequence, structure, delta, stdout);
      if (dos) {
        int i;
        for (i=0; i<= MAXDOS && i<=delta/10; i++) {
          printf("%4d %6d\n", i, density_of_states[i]);
        }
      }
    }
    /* Zuker suboptimals */
    else{
      SOLUTION *zr;
      int i;
      if (cut_point!=-1) {
        nrerror("Sorry, zuker subopts not yet implemented for cofold\n");
      }
      zr = zukersubopt(sequence);
      putoutzuker(zr);
      (void)fflush(stdout);
      for (i=0; zr[i].structure; i++) {
        free(zr[i].structure);
      }
      free(zr);
    }
    (void)fflush(stdout);
    free(sequence);
    free(structure);
  } while (1);
  return 0;
}

PRIVATE char *tokenize(char *line)
{
  char *pos, *copy;
  int cut = -1;

  copy = (char *) space(strlen(line)+1);
  (void) sscanf(line, "%s", copy);
  pos = strchr(copy, '&');
  if (pos) {
    cut = (int) (pos-copy)+1;
    if (cut >= strlen(copy)) cut = -1;
    if (strchr(pos+1, '&')) nrerror("more than one cut-point in input");
    for (;*pos;pos++) *pos = *(pos+1); /* splice out the & */
  }
  if (cut > -1) {
    if (cut_point==-1) cut_point = cut;
    else if (cut_point != cut) {
      fprintf(stderr,"cut_point = %d cut = %d\n", cut_point, cut);
      nrerror("Sequence and Structure have different cut points.");
    }
  }

  free(line);
  return copy;
}
PRIVATE void putoutzuker(SOLUTION* zukersolution) {
  int i;
  printf("%s [%.2f]\n",zukersolution[0].structure,zukersolution[0].energy/100.);
  for(i=1; zukersolution[i].structure; i++) {
    printf("%s [%.2f]\n", zukersolution[i].structure, zukersolution[i].energy/100.);
  }
  return;
}
