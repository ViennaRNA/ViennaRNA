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
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"

#define PRIVATE      static
#define PUBLIC
#define MAXWIDTH     201

PRIVATE float F[MAXWIDTH];
PRIVATE float ddiff(float f[], float h, int m);
PRIVATE void  usage(void);

#define GASCONST 1.98717  /* in [cal/K] */
#define K0 273.15

PRIVATE void heat_capacity(char *string, float T_min, float T_max,
			  float h, int m)
{
   extern int pf_dangl;
   int length, i;
   char *structure=NULL;
   float hc, kT, min_en;
   
   length = strlen(string);
   
   do_backtrack = 0;   

   temperature = T_min -m*h;
   initialize_fold(length);
   structure = (char *) space(length+1);
   min_en = fold(string, structure);
   pf_dangl=1;
   min_en = energy_of_struct(string,structure);
   free(structure); structure=NULL; free_arrays();
   kT = (temperature+K0)*GASCONST/1000;    /* in kcal */
   pf_scale = exp((1.07*min_en)/kT/length );
   init_pf_fold(length);
   
   for (i=0; i<2*m+1; i++) {
      F[i] = pf_fold(string, structure);   /* T_min -2h */
      temperature += h;
      kT = (temperature+K0)*GASCONST/1000;
      pf_scale=exp((F[i]/length)/kT);
      update_pf_params(length); 
   }
   while (temperature <= (T_max+m*h+h)) {
      
      hc = - ddiff(F,h,m)* (temperature +K0 - m*h -h); 
      printf("%g   %g\n", (temperature-m*h-h), hc);  
      
      for (i=0; i<2*m; i++)
	 F[i] = F[i+1];
      F[2*m] = pf_fold(string, structure); 
      temperature += h;
      kT = (temperature+K0)*GASCONST/1000;
      pf_scale=exp((F[i]/length)/kT);
      update_pf_params(length); 
   }
   free_pf_arrays();
}

/* -------------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------------- */

char  scale1[] = "....,....1....,....2....,....3....,....4";
char  scale2[] = "....,....5....,....6....,....7....,....8";

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   char *string, *line;
   char  ns_bases[33]="", *c;
   char  ParamFile[256]="";
   int  i, length, l, sym;
   float T_min, T_max, h;
   int m_points;
   int istty;
   
   T_min=0.; T_max=100.; h=1; m_points=2;
   string=NULL;
   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') 
	 switch ( argv[i][1] ) {
	  case 'T':
	    if (strncmp(argv[i], "-Tmin", 5)==0) {
	       if (sscanf(argv[++i], "%f", &T_min)==)
		 usage();
	    }
	    if (strncmp(argv[i], "-Tmax",5)==0) {
	       if (sscanf(argv[++i], "%f", &T_max)==2)
		 usage();
	    }
	    break;
	  case 'h':
	    if (i>argc-2) usage();
	    if (sscanf(argv[++i],"%f",&h)==2)
	      usage();
	    break;
	  case 'n':
	    if ( strcmp(argv[i], "-noGU" )==0) noGU=1;
	    if ( strcmp(argv[i], "-noCloseGU" ) ==0) no_closingGU=1;
	    if ( strcmp(argv[i], "-nsp") ==0) {
	       if (sscanf(argv[++i], "%32s", ns_bases)==0)
		 usage();
	    }
	    break;
	  case '4':
	    tetra_loop=0;
	    break;
	  case 'e':
	    if (sscanf(argv[++i],"%d", &energy_set)==)
	      usage();
	    break;
	  case 'm':
	    if (sscanf(argv[++i],"%d", &m_points)==0)
	      usage();
	    if (m_points<1) m_points=1;
	    if (m_points>100) m_points=100;
	    break;
	  case 'd': dangles=0;
	    break;
	  case 'P':
	    if (sscanf(argv[++i], "%255s", ParamFile)==0)
	      usage();
	    break;
	  default: usage();
	 }
   }

   if (ParamFile[0])
     read_parameter_file(ParamFile);
   
   if (ns_bases[0]) {
     nonstandards = space(33);
     c=ns_bases;
     i=sym=0;
     if (*c=='-') {
       sym=1; c++;
     }
     while (*c) {
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
   
   istty = isatty(fileno(stdout));
   
   do {
      if (istty) {
	 printf("\nInput string (upper or lower case); @ to quit\n");
	 printf("%s%s\n", scale1, scale2);
      }
      
      if ((line = get_line(stdin))==NULL) break;
      
      /* skip comment lines and get filenames */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	 printf("%s\n", line);
	 free(line);
	 if ((line = get_line(stdin))==NULL) line = "@";
      } 
      
      if (strcmp(line, "@") == 0) break;
      
      string = (char *) space(strlen(line)+1);
      sscanf(line,"%s",string);
      free(line);
      length = strlen(string);
       
      for (l = 0; l < length; l++) string[l] = toupper(string[l]);
      if (istty)
	 printf("length = %d\n", length);
      
      heat_capacity(string, T_min, T_max, h, m_points);
      fflush(stdout);
   } while (1);
   return 0;
}

PRIVATE void usage(void)
{
   nrerror("usage: RNAheat [-Tmin t1] [-Tmax t2] [-h stepsize] "
	   "[-m ipoints] [-4]\n"
	   "               [-noGU] [noCloseGU] [-d] [-nsp pairs] "
	   "[-e e_set] \n");
}
