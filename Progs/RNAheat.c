/*
                      Heat Capacity of RNA molecule

                    c Ivo Hofacker and Peter Stadler
                          Vienna RNA package


            calculates specific heat using C = - T d^2/dT^2 G(T)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "read_epars.h"
#include "RNAheat_cmdl.h"

/*@unused@*/
static char rcsid[] = "$Id: RNAheat.c,v 1.12 2000/09/28 11:23:14 ivo Rel $";

#define MAXWIDTH  201
#define GASCONST  1.98717  /* in [cal/K] */
#define K0        273.15

PRIVATE float F[MAXWIDTH];
PRIVATE float ddiff(float f[], float h, int m);
PRIVATE void heat_capacity(char *string, float T_min, float T_max, float h, int m);

int main(int argc, char *argv[]){
  struct  RNAheat_args_info args_info;
  char                      *string, *input_string, *ns_bases, *c, *ParamFile;
  int                       i, length, l, sym;
  float                     T_min, T_max, h;
  int                       mpoints, istty, noconv = 0;
  unsigned int              input_type;

  string = ParamFile = ns_bases = NULL;
  T_min    = 0.;
  T_max    = 100.;
  h        = 1;
  mpoints  = 2;
  dangles  = 2;   /* dangles can be 0 (no dangles) or 2, default is 2 */

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAheat_cmdline_parser(argc, argv, &args_info) != 0) exit(1);

  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given){
    dangles = args_info.dangles_arg;
    if(dangles % 2){
      warn_user("using default dangles = 2");
      dangles = 2;
    }
  }
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
  /* Tmin */
  if(args_info.Tmin_given)        T_min = args_info.Tmin_arg;
  /* Tmax */
  if(args_info.Tmax_given)        T_max = args_info.Tmax_arg;
  /* step size */
  if(args_info.stepsize_given)    h = args_info.stepsize_arg;
  /* ipoints */
  if(args_info.ipoints_given){
    mpoints = args_info.ipoints_arg;
    if (mpoints < 1)    mpoints = 1;
    if (mpoints > 100)  mpoints = 100;
  }

  /* free allocated memory of command line data structure */
  RNAheat_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  if (ParamFile!=NULL) read_parameter_file(ParamFile);

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

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do {
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty) print_tty_input_seq();

    /* skip fasta headers and comments */
    while((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string);
      /* skip fasta headers and comments */
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      string = strdup(input_string);
      length = (int) strlen(string);
      free(input_string);
    }

    if(noconv)  str_RNA2RNA(string);
    else        str_DNA2RNA(string);

    if (istty) printf("length = %d\n", length);

    heat_capacity(string, T_min, T_max, h, mpoints);
    free(string);
    (void) fflush(stdout);
  } while (1);
  return 0;
}

PRIVATE void heat_capacity(char *string, float T_min, float T_max,
                          float h, int m)
{
   int length, i;
   char *structure;
   float hc, kT, min_en;
   
   length = (int) strlen(string);
   
   do_backtrack = 0;   

   temperature = T_min -m*h;
   /* initialize_fold(length); <- obsolete */
   structure = (char *) space((unsigned) length+1);
   min_en = fold(string, structure);
   free(structure); free_arrays();
   kT = (temperature+K0)*GASCONST/1000;    /* in kcal */
   pf_scale = exp(-(1.07*min_en)/kT/length );
   /* init_pf_fold(length); <- obsolete */
   
   for (i=0; i<2*m+1; i++) {
      F[i] = pf_fold(string, NULL);   /* T_min -2h */
      temperature += h;
      kT = (temperature+K0)*GASCONST/1000;
      pf_scale=exp(-(F[i]/length +h*0.00727)/kT); /* try to extrapolate F */
      update_pf_params(length); 
   }
   while (temperature <= (T_max+m*h+h)) {
      
      hc = - ddiff(F,h,m)* (temperature +K0 - m*h -h); 
      printf("%g   %g\n", (temperature-m*h-h), hc);  
      
      for (i=0; i<2*m; i++)
         F[i] = F[i+1];
      F[2*m] = pf_fold(string, NULL); 
      temperature += h;
      kT = (temperature+K0)*GASCONST/1000;
      pf_scale=exp(-(F[i]/length +h*0.00727)/kT);
      update_pf_params(length); 
   }
   free_pf_arrays();
}

/* ------------------------------------------------------------------------- */

PRIVATE float ddiff(float f[], float h, int m)
{
   float fp;
   int i;
   float A, B;
   A = (float)(m*(m+1)*(2*m+1)/3 );                            /* 2*sum(x^2) */
   B = (float)(m*(m+1)*(2*m+1) ) * (float)(3*m*m+3*m-1) /15.;  /* 2*sum(x^4) */
   
   fp=0.;
   for (i=0; i<2*m+1; i++)
      fp += f[i]*( A - (float) ( (2*m+1)*(i-m)*(i-m)) );
   
   fp /= ( ( A*A - B*( (float)(2*m+1) ) )*h*h/2. );
   return (float)fp;
   
}

