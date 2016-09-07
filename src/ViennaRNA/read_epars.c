/*
                  read energy parameters from a file

                      Stephan Kopp, Ivo Hofacker
                          Vienna RNA Package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_const.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/read_epars.h"

#define PUBLIC
#define PRIVATE   static
#define PARSET 20

#define DEF -50
#define NST 0

PRIVATE FILE *fp;

PRIVATE void  display_array(int *p, int size, int line, FILE *fp);
PRIVATE char  *get_array1(int *arr, int size);
PRIVATE void  ignore_comment(char *line);
PRIVATE void  check_symmetry(void);
PRIVATE void  update_nst(int array[NBPAIRS+1][NBPAIRS+1][5][5][5][5]);

/**
*** read a 1dimensional array from file
*** \param array  a pointer to the first element in the array
*** \param dim    the size of the array
*** \param shift  the first position the new values will be written in
**/
PRIVATE void  rd_1dim(int *array, int dim, int shift);
PRIVATE void  rd_1dim_slice(int *array, int dim, int shift, int post);
PRIVATE void  rd_2dim(int *array,
                      int dim1, int dim2,
                      int shift1, int shift2);
PRIVATE void  rd_2dim_slice(int *array,
                      int dim1, int dim2,
                      int shift1, int shift2,
                      int post1, int post2);
PRIVATE void  rd_3dim(int *array,
                      int dim1, int dim2, int dim3,
                      int shift1, int shift2, int shift3);
PRIVATE void  rd_3dim_slice(int *array,
                      int dim1, int dim2, int dim3,
                      int shift1, int shift2, int shift3,
                      int post1, int post2, int post3);
PRIVATE void  rd_4dim(int *array,
                      int dim1, int dim2, int dim3, int dim4,
                      int shift1, int shift2, int shift3, int shift4);
PRIVATE void  rd_4dim_slice(int *array,
                      int dim1, int dim2, int dim3, int dim4,
                      int shift1, int shift2, int shift3, int shift4,
                      int post1, int post2, int post3, int post4);
PRIVATE void  rd_5dim(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5,
                      int shift1, int shift2, int shift3, int shift4, int shift5);
PRIVATE void  rd_5dim_slice(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5,
                      int shift1, int shift2, int shift3, int shift4, int shift5,
                      int post1, int post2, int post3, int post4, int post5);
PRIVATE void  rd_6dim(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5, int dim6,
                      int shift1, int shift2, int shift3, int shift4, int shift5, int shift6);
PRIVATE void  rd_6dim_slice(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5, int dim6,
                      int shift1, int shift2, int shift3, int shift4, int shift5, int shift6,
                      int post1, int post2, int post3, int post4, int post5, int post6);
PRIVATE void  rd_Tetraloop37(void);
PRIVATE void  rd_Triloop37(void);
PRIVATE void  rd_Hexaloop37(void);

/*------------------------------------------------------------*/
PUBLIC void read_parameter_file(const char fname[]){
  char        *line, ident[256];
  enum parset type;
  int         r;

  if (!(fp=fopen(fname,"r"))) {
    vrna_message_warning("\nread_parameter_file:\n"
                                "\t\tcan't open file %s\n"
                                "\t\tusing default parameters instead.",
                                fname);
    return;
  }

  if (!(line = vrna_read_line(fp))) {
    vrna_message_warning(" File %s is improper.\n", fname);
    fclose(fp);
    return;
  }

  if (strncmp(line,"## RNAfold parameter file v2.0",30)!=0) {
    vrna_message_warning( "Missing header line in file.\n"
                          "May be this file has not v2.0 format.\n"
                          "Use INTERRUPT-key to stop.");
  }
  free(line);

  while((line=vrna_read_line(fp))) {

    r = sscanf(line, "# %255s", ident);
    if (r==1) {
      type = gettype(ident);
      switch (type){
        case QUIT:    break;
        case S:       rd_2dim(&(stack37[0][0]), NBPAIRS+1, NBPAIRS+1, 1, 1);
                      break;
        case S_H:     rd_2dim(&(stackdH[0][0]), NBPAIRS+1, NBPAIRS+1, 1, 1);
                      break;
        case HP:      rd_1dim(&(hairpin37[0]), 31, 0);
                      break;
        case HP_H:    rd_1dim(&(hairpindH[0]), 31, 0);
                      break;
        case B:       rd_1dim(&(bulge37[0]), 31, 0);
                      break;
        case B_H:     rd_1dim(&(bulgedH[0]), 31, 0);
                      break;
        case IL:      rd_1dim(&(internal_loop37[0]), 31, 0);
                      break;
        case IL_H:    rd_1dim(&(internal_loopdH[0]), 31, 0);
                      break;
        case MME:     rd_3dim(&(mismatchExt37[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MME_H:   rd_3dim(&(mismatchExtdH[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMH:     rd_3dim(&(mismatchH37[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMH_H:   rd_3dim(&(mismatchHdH[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMI:     rd_3dim(&(mismatchI37[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMI_H:   rd_3dim(&(mismatchIdH[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMI1N:   rd_3dim(&(mismatch1nI37[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMI1N_H: rd_3dim(&(mismatch1nIdH[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMI23:   rd_3dim(&(mismatch23I37[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMI23_H: rd_3dim(&(mismatch23IdH[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMM:     rd_3dim(&(mismatchM37[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case MMM_H:   rd_3dim(&(mismatchMdH[0][0][0]),
                          NBPAIRS+1, 5, 5,
                          1, 0, 0);
                      break;
        case INT11:   rd_4dim(&(int11_37[0][0][0][0]),
                          NBPAIRS+1, NBPAIRS+1, 5, 5,
                          1, 1, 0, 0);
                      break;
        case INT11_H: rd_4dim(&(int11_dH[0][0][0][0]),
                          NBPAIRS+1, NBPAIRS+1, 5, 5,
                          1, 1, 0, 0);
                      break;
        case INT21:   rd_5dim(&(int21_37[0][0][0][0][0]),
                          NBPAIRS+1, NBPAIRS+1, 5, 5, 5,
                          1, 1, 0, 0, 0);
                      break;
        case INT21_H: rd_5dim(&(int21_dH[0][0][0][0][0]),
                          NBPAIRS+1, NBPAIRS+1, 5, 5, 5,
                          1, 1, 0, 0, 0);
                      break;
        case INT22:   rd_6dim_slice(&(int22_37[0][0][0][0][0][0]),
                          NBPAIRS+1, NBPAIRS+1, 5, 5, 5, 5,
                          1, 1, 1, 1, 1, 1,
                          1, 1, 0, 0, 0, 0);
                      update_nst(int22_37);
                      break;
        case INT22_H: rd_6dim_slice(&(int22_dH[0][0][0][0][0][0]),
                          NBPAIRS+1, NBPAIRS+1, 5, 5, 5, 5,
                          1, 1, 1, 1, 1, 1,
                          1, 1, 0, 0, 0, 0);
                      update_nst(int22_dH);
                      break;
        case D5:      rd_2dim(&(dangle5_37[0][0]), NBPAIRS+1, 5, 1, 0);
                      break;
        case D5_H:    rd_2dim(&(dangle5_dH[0][0]), NBPAIRS+1, 5, 1, 0);
                      break;
        case D3:      rd_2dim(&(dangle3_37[0][0]), NBPAIRS+1, 5, 1, 0);
                      break;
        case D3_H:    rd_2dim(&(dangle3_dH[0][0]), NBPAIRS+1, 5, 1, 0);
                      break;
        case ML:      {
                        int values[6];
                        rd_1dim(&values[0], 6, 0);
                        ML_BASE37     = values[0];
                        ML_BASEdH     = values[1];
                        ML_closing37  = values[2];
                        ML_closingdH  = values[3];
                        ML_intern37   = values[4];
                        ML_interndH   = values[5];
                      }
                      break;
        case NIN:     {
                        int values[3];
                        rd_1dim(&values[0], 3, 0);
                        ninio37 = values[0];
                        niniodH = values[1];
                        MAX_NINIO  = values[2];
                      }
                      break;
        case MISC:    {
                        int values[4];
                        rd_1dim(&values[0], 4, 0);
                        DuplexInit37 = values[0];
                        DuplexInitdH = values[1];
                        TerminalAU37 = values[2];
                        TerminalAUdH = values[3];
                      }
                      break;
        case TL:      rd_Tetraloop37();
                      break;
        case TRI:     rd_Triloop37();
                      break;
        case HEX:     rd_Hexaloop37();
                      break;
        default:      /* do nothing but complain */
                      vrna_message_warning("read_epars: Unknown field identifier in `%s'", line);
      }
    } /* else ignore line */
    free(line);
  }
  fclose(fp);

  check_symmetry();
  return;
}

/*------------------------------------------------------------*/

PRIVATE void display_array(int *p, int size, int nl, FILE *fp)
{
  int i;

  for (i=1; i<=size; i++, p++) {
    switch(*p)
      {
      case  INF: fprintf(fp,"   INF"); break;
      case -INF: fprintf(fp,"  -INf"); break;
      case  DEF: fprintf(fp,"   DEF"); break;
      default:   fprintf(fp,"%6d",  *p); break;
      }
    if ((i%nl)==0) fprintf(fp,"\n");
  }
  if (size%nl) fprintf(fp,"\n");
  return;
}

/*------------------------------------------------------------*/

PRIVATE char *get_array1(int *arr, int size)
{
  int    i, p, pos, pp, r, last;
  char  *line, buf[16];


  i = last = 0;
  while( i<size ) {
    line = vrna_read_line(fp);
    if (!line) vrna_message_error("unexpected end of file in get_array1");
    ignore_comment(line);
    pos=0;
    while ((i<size)&&(sscanf(line+pos,"%15s%n", buf, &pp)==1)) {
      pos += pp;
      if (buf[0]=='*') {i++; continue;}
      else if (buf[0]=='x') { /* should only be used for loop parameters */
        if (i==0) vrna_message_error("can't extrapolate first value");
        p = arr[last] + (int) (0.5+ lxc37*log(((double) i)/(double)(last)));
      }
      else if (strcmp(buf,"DEF") == 0) p = DEF;
      else if (strcmp(buf,"INF") == 0) p = INF;
      else if (strcmp(buf,"NST") == 0) p = NST;
      else {
        r=sscanf(buf,"%d", &p);
        if (r!=1) {
          return line+pos;
          vrna_message_error("can't interpret `%s' in get_array1", buf);
          exit(1);
        }
        last = i;
      }
      arr[i++]=p;
    }
    free(line);
  }

  return NULL;
}

PRIVATE void rd_1dim(int *array, int dim, int shift){
  rd_1dim_slice(array, dim, shift, 0);
}

PRIVATE void rd_1dim_slice(int *array, int dim, int shift, int post){
  char *cp;
  cp   = get_array1(array+shift, dim-shift-post);

  if (cp) {
    vrna_message_error("\nrd_1dim: %s", cp);
    exit(1);
  }
  return;
}

PRIVATE void  rd_2dim(int *array, int dim1, int dim2, int shift1, int shift2){
  rd_2dim_slice(array, dim1, dim2, shift1, shift2, 0, 0);
}

PRIVATE void  rd_2dim_slice(int *array,
                      int dim1, int dim2,
                      int shift1, int shift2,
                      int post1, int post2){
  int i;
  int delta_pre   = shift1 + shift2;
  int delta_post  = post1 + post2;

  if(delta_pre + delta_post == 0){
    rd_1dim(array, dim1 * dim2, 0);
    return;
  }
  for (i=shift1; i<dim1 - post1; i++)
    rd_1dim_slice(array + (i*dim2), dim2, shift2, post2);
  return;
}

PRIVATE void  rd_3dim(int *array, int dim1, int dim2, int dim3, int shift1, int shift2, int shift3){
  rd_3dim_slice(array,
                dim1, dim2, dim3,
                shift1, shift2, shift3,
                0, 0, 0);
}

PRIVATE void  rd_3dim_slice(int *array,
                            int dim1, int dim2, int dim3,
                            int shift1, int shift2, int shift3,
                            int post1, int post2, int post3){
  int    i;
  int delta_pre   = shift1 + shift2 + shift3;
  int delta_post  = post1 + post2 + post3;

  if(delta_pre + delta_post == 0){
    rd_1dim(array, dim1 * dim2 * dim3, 0);
    return;
  }
  for (i=shift1; i<dim1 - post1; i++){
    rd_2dim_slice(array + (i * dim2 * dim3),
            dim2, dim3,
            shift2, shift3,
            post2, post3);
  }
  return;
}

PRIVATE void  rd_4dim(int *array,
                      int dim1, int dim2, int dim3, int dim4,
                      int shift1, int shift2, int shift3, int shift4){
  rd_4dim_slice(array,
                dim1, dim2, dim3, dim4,
                shift1, shift2, shift3, shift4,
                0, 0, 0, 0);
}

PRIVATE void  rd_4dim_slice(int *array,
                      int dim1, int dim2, int dim3, int dim4,
                      int shift1, int shift2, int shift3, int shift4,
                      int post1, int post2, int post3, int post4){
  int i;
  int delta_pre   = shift1 + shift2 + shift3 + shift4;
  int delta_post  = post1 + post2 + post3 + post4;

  if(delta_pre + delta_post == 0){
    rd_1dim(array, dim1 * dim2 * dim3 * dim4, 0);
    return;
  }
  for(i=shift1; i<dim1 - post1; i++){
    rd_3dim_slice(array + (i * dim2 * dim3 * dim4),
            dim2, dim3, dim4,
            shift2, shift3, shift4,
            post2, post3, post4);
  }
  return;
}

PRIVATE void  rd_5dim(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5,
                      int shift1, int shift2, int shift3, int shift4, int shift5){
  rd_5dim_slice(array,
                dim1, dim2, dim3, dim4, dim5,
                shift1, shift2, shift3, shift4, shift5,
                0, 0, 0, 0, 0);
}

PRIVATE void  rd_5dim_slice(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5,
                      int shift1, int shift2, int shift3, int shift4, int shift5,
                      int post1, int post2, int post3, int post4, int post5){
  int i;
  int delta_pre   = shift1 + shift2 + shift3 + shift4 + shift5;
  int delta_post  = post1 + post2 + post3 + post4 + post5;

  if(delta_pre + delta_post == 0){
    rd_1dim(array, dim1 * dim2 * dim3 * dim4 * dim5, 0);
    return;
  }
  for(i=shift1; i<dim1 - post1; i++)
    rd_4dim_slice(array + (i * dim2 * dim3 * dim4 * dim5),
            dim2, dim3, dim4, dim5,
            shift2, shift3, shift4, shift5,
            post2, post3, post4, post5);
  return;
}

/**
*** \param dim1   The size of the first dimension
*** \param shift1 The pre shift for the first dimension
**/
PRIVATE void  rd_6dim(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5, int dim6,
                      int shift1, int shift2, int shift3, int shift4, int shift5, int shift6){
  rd_6dim_slice(array,
                dim1, dim2, dim3, dim4, dim5, dim6,
                shift1, shift2, shift3, shift4, shift5, shift6,
                0, 0, 0, 0, 0, 0);
}

PRIVATE void  rd_6dim_slice(int *array,
                      int dim1, int dim2, int dim3, int dim4, int dim5, int dim6,
                      int shift1, int shift2, int shift3, int shift4, int shift5, int shift6,
                      int post1, int post2, int post3, int post4, int post5, int post6){
  int i;
  int delta_pre   = shift1 + shift2 + shift3 + shift4 + shift5 + shift6;
  int delta_post  = post1 + post2 + post3 + post4 + post5 + post6;

  if(delta_pre + delta_post == 0){
    rd_1dim(array, dim1 * dim2 * dim3 * dim4 * dim5 * dim6, 0);
    return;
  }
  for(i=shift1; i<dim1 - post1; i++)
    rd_5dim_slice(array + (i * dim2 * dim3 * dim4 * dim5 * dim6),
            dim2, dim3, dim4, dim5, dim6,
            shift2, shift3, shift4, shift5, shift6,
            post2, post3, post4, post5, post6);
  return;
}


/*------------------------------------------------------------*/
PRIVATE void  rd_Tetraloop37(void)
{
  int    i, r;
  char   *buf;

  i=0;
  /* erase old tetraloop entries */
  memset(&Tetraloops, 0, 281);
  memset(&Tetraloop37, 0, sizeof(int)*40);
  memset(&TetraloopdH, 0, sizeof(int)*40);
  do {
    buf = vrna_read_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%6s %d %d", &Tetraloops[7*i], &Tetraloop37[i], &TetraloopdH[i]);
    strcat(Tetraloops, " ");
    free(buf);
    i++;
  } while((r==3)&&(i<40));
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_Hexaloop37(void)
{
  int    i, r;
  char   *buf;

  i=0;
  /* erase old hexaloop entries */
  memset(&Hexaloops, 0, 361);
  memset(&Hexaloop37, 0, sizeof(int)*40);
  memset(&HexaloopdH, 0, sizeof(int)*40);
  do {
    buf = vrna_read_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%8s %d %d", &Hexaloops[9*i], &Hexaloop37[i], &HexaloopdH[i]);
    strcat(Hexaloops, " ");
    free(buf);
    i++;
  } while((r==3)&&(i<40));
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_Triloop37(void)
{
  int    i, r;
  char   *buf;

  i=0;
  /* erase old hexaloop entries */
  memset(&Triloops,   0, 241);
  memset(&Triloop37,  0, sizeof(int)*40);
  memset(&TriloopdH,  0, sizeof(int)*40);
  do {
    buf = vrna_read_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%5s %d %d", &Triloops[6*i], &Triloop37[i], &TriloopdH[i]);
    strcat(Triloops, " ");
    free(buf);
    i++;
  } while((r==3)&&(i<40));
  return;
}

/*------------------------------------------------------------*/


PRIVATE void ignore_comment(char * line)
{
  /* excise C style comments */
  /* only one comment per line, no multiline comments */
  char *cp1, *cp2;

  if ((cp1=strstr(line, "/*"))) {
    cp2 = strstr(cp1, "*/");
    if (cp2==NULL)
      vrna_message_error("unclosed comment in parameter file");
    /* can't use strcpy for overlapping strings */
    for (cp2+=2; *cp2!='\0'; cp2++, cp1++)
      *cp1 = *cp2;
    *cp1 = '\0';
  }

  return;
}
/*------------------------------------------------------------*/

PUBLIC char *settype(enum parset s){
  switch(s){
    case        S:  return "stack";
    case      S_H:  return "stack_enthalpies";
    case       HP:  return "hairpin";
    case     HP_H:  return "hairpin_enthalpies";
    case        B:  return "bulge";
    case      B_H:  return "bulge_enthalpies";
    case       IL:  return "interior";
    case     IL_H:  return "interior_enthalpies";
    case      MME:  return "mismatch_exterior";
    case    MME_H:  return "mismatch_exterior_enthalpies";
    case      MMH:  return "mismatch_hairpin";
    case    MMH_H:  return "mismatch_hairpin_enthalpies";
    case      MMI:  return "mismatch_interior";
    case    MMI_H:  return "mismatch_interior_enthalpies";
    case    MMI1N:  return "mismatch_interior_1n";
    case  MMI1N_H:  return "mismatch_interior_1n_enthalpies";
    case    MMI23:  return "mismatch_interior_23";
    case  MMI23_H:  return "mismatch_interior_23_enthalpies";
    case      MMM:  return "mismatch_multi";
    case    MMM_H:  return "mismatch_multi_enthalpies";
    case       D5:  return "dangle5";
    case     D5_H:  return "dangle5_enthalpies";
    case       D3:  return "dangle3";
    case     D3_H:  return "dangle3_enthalpies";
    case    INT11:  return "int11";
    case  INT11_H:  return "int11_enthalpies";
    case    INT21:  return "int21";
    case  INT21_H:  return "int21_enthalpies";
    case    INT22:  return "int22";
    case  INT22_H:  return "int22_enthalpies";
    case       ML:  return "ML_params";
    case      NIN:  return "NINIO";
    case      TRI:  return "Triloops";
    case       TL:  return "Tetraloops";
    case      HEX:  return "Hexaloops";
    case     QUIT:  return "END";
    case     MISC:  return "Misc";
    default: vrna_message_error("\nThe answer is: 42\n");
  }
  return "";
}
/*------------------------------------------------------------*/

PUBLIC enum parset gettype(const char *ident){
  if      (strcmp(ident,"stack") == 0)                            return S;
  else if (strcmp(ident,"stack_enthalpies") == 0)                 return S_H;
  else if (strcmp(ident,"hairpin") == 0)                          return HP;
  else if (strcmp(ident,"hairpin_enthalpies") == 0)               return HP_H;
  else if (strcmp(ident,"bulge") == 0)                            return B;
  else if (strcmp(ident,"bulge_enthalpies") == 0)                 return B_H;
  else if (strcmp(ident,"interior") == 0)                         return IL;
  else if (strcmp(ident,"interior_enthalpies") == 0)              return IL_H;
  else if (strcmp(ident,"mismatch_exterior") == 0)                return MME;
  else if (strcmp(ident,"mismatch_exterior_enthalpies") == 0)     return MME_H;
  else if (strcmp(ident,"mismatch_hairpin") == 0)                 return MMH;
  else if (strcmp(ident,"mismatch_hairpin_enthalpies") == 0)      return MMH_H;
  else if (strcmp(ident,"mismatch_interior") == 0)                return MMI;
  else if (strcmp(ident,"mismatch_interior_enthalpies") == 0)     return MMI_H;
  else if (strcmp(ident,"mismatch_interior_1n") == 0)             return MMI1N;
  else if (strcmp(ident,"mismatch_interior_1n_enthalpies") == 0)  return MMI1N_H;
  else if (strcmp(ident,"mismatch_interior_23") == 0)             return MMI23;
  else if (strcmp(ident,"mismatch_interior_23_enthalpies") == 0)  return MMI23_H;
  else if (strcmp(ident,"mismatch_multi") == 0)                   return MMM;
  else if (strcmp(ident,"mismatch_multi_enthalpies") == 0)        return MMM_H;
  else if (strcmp(ident,"int11") == 0)                            return INT11;
  else if (strcmp(ident,"int11_enthalpies") == 0)                 return INT11_H;
  else if (strcmp(ident,"int21") == 0)                            return INT21;
  else if (strcmp(ident,"int21_enthalpies") == 0)                 return INT21_H;
  else if (strcmp(ident,"int22") == 0)                            return INT22;
  else if (strcmp(ident,"int22_enthalpies") == 0)                 return INT22_H;
  else if (strcmp(ident,"dangle5")== 0)                           return D5;
  else if (strcmp(ident,"dangle5_enthalpies")== 0)                return D5_H;
  else if (strcmp(ident,"dangle3")== 0)                           return D3;
  else if (strcmp(ident,"dangle3_enthalpies")== 0)                return D3_H;
  else if (strcmp(ident,"ML_params")== 0)                         return ML;
  else if (strcmp(ident,"NINIO") == 0)                            return NIN;
  else if (strcmp(ident,"Triloops") == 0)                         return TRI;
  else if (strcmp(ident,"Tetraloops") == 0)                       return TL;
  else if (strcmp(ident,"Hexaloops") == 0)                        return HEX;
  else if (strcmp(ident,"Misc") == 0)                             return MISC;
  else if (strcmp(ident,"END") == 0)                              return QUIT;
  else return UNKNOWN;
}

/*---------------------------------------------------------------*/

PUBLIC void write_parameter_file(const char fname[]){
  FILE *outfp;
  int c;
  char *pnames[] = {"NP", "CG", "GC", "GU", "UG", "AU", "UA", " @"};
  char bnames[] = "@ACGU";
  outfp = fopen(fname, "w");
  if (!outfp) {
    vrna_message_error("can't open file %s", fname);
    exit(1);
  }
  fprintf(outfp,"## RNAfold parameter file v2.0\n");

  fprintf(outfp,"\n# %s\n", settype(S));
  fprintf(outfp,"/*  CG    GC    GU    UG    AU    UA    @  */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(stack37[c]+1,NBPAIRS,NBPAIRS, outfp);

  fprintf(outfp,"\n# %s\n", settype(S_H));
  fprintf(outfp,"/*  CG    GC    GU    UG    AU    UA    @  */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(stackdH[c]+1,NBPAIRS,NBPAIRS, outfp);

  fprintf(outfp,"\n# %s\n", settype(MMH));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchH37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMH_H));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchHdH[k][i],5,5, outfp);

  }

  fprintf(outfp,"\n# %s\n", settype(MMI));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchI37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMI_H));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchIdH[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMI1N));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatch1nI37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMI1N_H));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatch1nIdH[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMI23));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatch23I37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMI23_H));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatch23IdH[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMM));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchM37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MMM_H));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchMdH[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MME));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchExt37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(MME_H));
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++)
      display_array(mismatchExtdH[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# %s\n", settype(D5));
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(dangle5_37[c], 5, 5, outfp);

  fprintf(outfp,"\n# %s\n", settype(D5_H));
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(dangle5_dH[c], 5, 5, outfp);

  fprintf(outfp,"\n# %s\n", settype(D3));
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(dangle3_37[c], 5, 5, outfp);

  fprintf(outfp,"\n# %s\n", settype(D3_H));
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(dangle3_dH[c], 5, 5, outfp);


  /* dont print "no pair" entries for interior loop arrays */
  fprintf(outfp,"\n# %s\n", settype(INT11));
  { int i,k,l;
  for (k=1; k<NBPAIRS+1; k++)
    for (l=1; l<NBPAIRS+1; l++) {
      fprintf(outfp, "/* %2s..%2s */\n", pnames[k], pnames[l]);
      for (i=0; i<5; i++)
        display_array(int11_37[k][l][i], 5, 5, outfp);
    }
  }

  fprintf(outfp,"\n# %s\n", settype(INT11_H));
  { int i,k,l;
  for (k=1; k<NBPAIRS+1; k++)
    for (l=1; l<NBPAIRS+1; l++) {
      fprintf(outfp, "/* %2s..%2s */\n", pnames[k], pnames[l]);
      for (i=0; i<5; i++)
        display_array(int11_dH[k][l][i],5,5, outfp);
    }
  }

  fprintf(outfp,"\n# %s\n", settype(INT21));
  { int p1, p2, i, j;
  for (p1=1; p1<NBPAIRS+1; p1++)
    for (p2=1; p2<NBPAIRS+1; p2++)
      for (i=0; i<5; i++) {
        fprintf(outfp, "/* %2s.%c..%2s */\n",
                pnames[p1], bnames[i], pnames[p2]);
        for (j=0; j<5; j++)
          display_array(int21_37[p1][p2][i][j],5,5, outfp);
      }
  }

  fprintf(outfp,"\n# %s\n", settype(INT21_H));
  { int p1, p2, i, j;
  for (p1=1; p1<NBPAIRS+1; p1++)
    for (p2=1; p2<NBPAIRS+1; p2++)
      for (i=0; i<5; i++) {
        fprintf(outfp, "/* %2s.%c..%2s */\n",
                pnames[p1], bnames[i], pnames[p2]);
        for (j=0; j<5; j++)
          display_array(int21_dH[p1][p2][i][j],5,5, outfp);
      }
  }

  fprintf(outfp,"\n# %s\n", settype(INT22));
  { int p1, p2, i, j, k;
  for (p1=1; p1<NBPAIRS; p1++)
    for (p2=1; p2<NBPAIRS; p2++)
      for (i=1; i<5; i++)
        for (j=1; j<5; j++) {
          fprintf(outfp, "/* %2s.%c%c..%2s */\n",
                  pnames[p1], bnames[i], bnames[j], pnames[p2]);
          for (k=1; k<5; k++)
            display_array(int22_37[p1][p2][i][j][k]+1,4,5, outfp);
        }
  }

  fprintf(outfp,"\n# %s\n", settype(INT22_H));
  { int p1, p2, i, j, k;
  for (p1=1; p1<NBPAIRS; p1++)
    for (p2=1; p2<NBPAIRS; p2++)
      for (i=1; i<5; i++)
        for (j=1; j<5; j++) {
          fprintf(outfp, "/* %2s.%c%c..%2s */\n",
                  pnames[p1], bnames[i], bnames[j], pnames[p2]);
          for (k=1; k<5; k++)
            display_array(int22_dH[p1][p2][i][j][k]+1,4,5, outfp);
        }
  }

  fprintf(outfp,"\n# %s\n", settype(HP));
  display_array(hairpin37, 31, 10, outfp);

  fprintf(outfp,"\n# %s\n", settype(HP_H));
  display_array(hairpindH, 31, 10, outfp);

  fprintf(outfp,"\n# %s\n", settype(B));
  display_array(bulge37, 31, 10, outfp);

  fprintf(outfp,"\n# %s\n", settype(B_H));
  display_array(bulgedH, 31, 10, outfp);

  fprintf(outfp,"\n# %s\n", settype(IL));
  display_array(internal_loop37, 31, 10, outfp);

  fprintf(outfp,"\n# %s\n", settype(IL_H));
  display_array(internal_loopdH, 31, 10, outfp);

  fprintf(outfp,"\n# %s\n", settype(ML));
  fprintf(outfp,"/* F = cu*n_unpaired + cc + ci*loop_degree (+TermAU) */\n");
  fprintf(outfp,"/*\t    cu\t cu_dH\t    cc\t cc_dH\t    ci\t ci_dH  */\n");
  fprintf(outfp,"\t%6d\t%6d\t%6d\t%6d\t%6d\t%6d\n", ML_BASE37, ML_BASEdH, ML_closing37, ML_closingdH, ML_intern37, ML_interndH);

  fprintf(outfp,"\n# %s\n", settype(NIN));
  fprintf(outfp,"/* Ninio = MIN(max, m*|n1-n2| */\n"
              "/*\t    m\t  m_dH     max  */\n"
              "\t%6d\t%6d\t%6d\n", ninio37, niniodH, MAX_NINIO);

  fprintf(outfp,"\n# %s\n", settype(MISC));
  fprintf(outfp,"/* all parameters are pairs of 'energy enthalpy' */\n");
  fprintf(outfp,"/*    DuplexInit     TerminalAU      LXC */\n");
  fprintf(outfp,"   %6d %6d %6d  %6d %3.6f %6d\n", DuplexInit37, DuplexInitdH, TerminalAU37, TerminalAUdH, lxc37, 0);

  fprintf(outfp,"\n# %s\n", settype(HEX));
  for (c=0; c< strlen(Hexaloops)/9; c++)
    fprintf(outfp,"\t%.8s %6d %6d\n", Hexaloops+c*9, Hexaloop37[c], HexaloopdH[c]);

  fprintf(outfp,"\n# %s\n", settype(TL));
  for (c=0; c< strlen(Tetraloops)/7; c++)
    fprintf(outfp,"\t%.6s %6d %6d\n", Tetraloops+c*7, Tetraloop37[c], TetraloopdH[c]);

  fprintf(outfp,"\n# %s\n", settype(TRI));
  for (c=0; c< strlen(Triloops)/6; c++)
    fprintf(outfp,"\t%.5s %6d %6d\n", Triloops+c*6, Triloop37[c], TriloopdH[c]);

  fprintf(outfp,"\n# %s\n", settype(QUIT));
  fclose(outfp);
}

PRIVATE void check_symmetry(void) {
  int i,j,k,l;

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      if (stack37[i][j] != stack37[j][i])
        vrna_message_warning("stacking energies not symmetric");

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      if (stackdH[i][j] != stackdH[j][i])
        vrna_message_warning("stacking enthalpies not symmetric");


  /* interior 1x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++)
          if (int11_37[i][j][k][l] != int11_37[j][i][l][k])
            vrna_message_warning("int11 energies not symmetric (%d,%d,%d,%d) (%d vs. %d)",
                                        i, j, k, l, int11_37[i][j][k][l], int11_37[j][i][l][k]);

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++)
          if (int11_dH[i][j][k][l] != int11_dH[j][i][l][k])
            vrna_message_warning("int11 enthalpies not symmetric");

  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)
              if (int22_37[i][j][k][l][m][n] != int22_37[j][i][m][n][k][l])
                vrna_message_warning("int22 energies not symmetric");
        }

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)
              if (int22_dH[i][j][k][l][m][n] != int22_dH[j][i][m][n][k][l])
                vrna_message_warning("int22 enthalpies not symmetric: %d %d %d %d %d %d",
                                            i,j,k,l,m,n);
        }
}

/* update nonstandard nucleotide/basepair involved contributions for int22 */
PRIVATE void update_nst(int array[NBPAIRS+1][NBPAIRS+1][5][5][5][5]){
  int    i, j, k, l, m, n;
  int max, max2, max3, max4, max5, max6;

  /* get maxima for one nonstandard nucleotide */
  for (i=1; i<NBPAIRS; i++){
    for (j=1; j<NBPAIRS; j++){
      for (k=1; k<5; k++){
        for (l=1; l<5; l++){
          for (m=1; m<5; m++){
            max = max2 = max3 = max4 = -INF; /* max of {CGAU} */
            for(n=1;n<5;n++){
              max   = MAX2(max,   array[i][j][k][l][m][n]);
              max2  = MAX2(max2,  array[i][j][k][l][n][m]);
              max3  = MAX2(max3,  array[i][j][k][n][l][m]);
              max4  = MAX2(max4,  array[i][j][n][k][l][m]);
            }
            array[i][j][k][l][m][0] = max;
            array[i][j][k][l][0][m] = max2;
            array[i][j][k][0][l][m] = max3;
            array[i][j][0][k][l][m] = max4;
          }
        }
      }
    }
  }
  /* get maxima for two nonstandard nucleotides */
  for (i=1; i<NBPAIRS; i++){
    for (j=1; j<NBPAIRS; j++){
      for (k=1; k<5; k++){
        for (l=1; l<5; l++){
          max = max2 = max3 = max4 = max5 = max6 = -INF; /* max of {CGAU} */
          for (m=1; m<5; m++){
            max   = MAX2(max,   array[i][j][k][l][m][0]);
            max2  = MAX2(max2,  array[i][j][k][m][0][l]);
            max3  = MAX2(max3,  array[i][j][m][0][k][l]);
            max4  = MAX2(max4,  array[i][j][0][k][l][m]);
            max5  = MAX2(max5,  array[i][j][0][k][m][l]);
            max6  = MAX2(max6,  array[i][j][k][0][l][m]);
          }
          array[i][j][k][l][0][0] = max;
          array[i][j][k][0][0][l] = max2;
          array[i][j][0][0][k][l] = max3;
          array[i][j][k][0][l][0] = max6;
          array[i][j][0][k][0][l] = max5;
          array[i][j][0][k][l][0] = max4;
        }
      }
    }
  }
  /* get maxima for three nonstandard nucleotides */
  for (i=1; i<NBPAIRS; i++){
    for (j=1; j<NBPAIRS; j++){
      for (k=1; k<5; k++){
        max = max2 = max3 = max4 = -INF; /* max of {CGAU} */
        for (l=1; l<5; l++){
          /* should be arbitrary where index l resides in last 3 possible locations */
          max   = MAX2(max,   array[i][j][k][l][0][0]);
          max2  = MAX2(max2,  array[i][j][0][k][l][0]);
          max3  = MAX2(max3,  array[i][j][0][0][k][l]);
          max4  = MAX2(max4,  array[i][j][0][0][l][k]);
        }
        array[i][j][k][0][0][0] = max;
        array[i][j][0][k][0][0] = max2;
        array[i][j][0][0][k][0] = max3;
        array[i][j][0][0][0][k] = max4;
      }
    }
  }
  /* get maxima for 4 nonstandard nucleotides */
  for (i=1; i<NBPAIRS; i++){
    for (j=1; j<NBPAIRS; j++){
      max = -INF; /* max of {CGAU} */
      for (k=1; k<5; k++){
        max   = MAX2(max,   array[i][j][k][0][0][0]);
      }
      array[i][j][0][0][0][0] = max;
    }
  }

  /* now compute contributions for nonstandard base pairs ... */
  /* first, 1 nonstandard bp */
  for (i=1; i<NBPAIRS; i++){
    for (k=0; k<5; k++){
      for (l=0; l<5; l++){
        for (m=0; m<5; m++){
          for(n=0;n<5;n++){
            max = max2 = -INF;
            for(j=1;j<NBPAIRS;j++){
              max   = MAX2(max, array[i][j][k][l][m][n]);
              max2  = MAX2(max2, array[j][i][k][l][m][n]);
            }
            array[i][NBPAIRS][k][l][m][n] = max;
            array[NBPAIRS][i][k][l][m][n] = max2;
          }
        }
      }
    }
  }

  /* now 2 nst base pairs */
  for (k=0; k<5; k++){
    for (l=0; l<5; l++){
      for (m=0; m<5; m++){
        for(n=0;n<5;n++){
          max = -INF;
          for(j=1;j<NBPAIRS;j++){
            max   = MAX2(max, array[NBPAIRS][j][k][l][m][n]);
          }
          array[NBPAIRS][NBPAIRS][k][l][m][n] = max;
        }
      }
    }
  }

}
