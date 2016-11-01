/*
###################################
# convert energy parameter files  #
# from ViennaRNAPackage 1.8.4 to  #
# 2.0 format                      #
#                                 #
# Ronny Lorenz                    #
###################################
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/pair_mat.h"

#include "1.8.4_epars.h"
#include "1.8.4_intloops.h"

#include "ViennaRNA/convert_epars.h"

enum parset_184 {UNKNOWN_184= -1, QUIT_184, S_184, SH_184, HP_184, B_184, IL_184, MMI_184, MMH_184, MMM_184, MM_H_184,
             DE5_184, DE3_184, DE5_H_184, DE3_H_184, ML_184, TL_184, TRI_184, TE_184, NIN_184, MISC_184,
             INT11_184, INT11_H_184, INT21_184, INT21_H_184, INT22_184, INT22_H_184};


PRIVATE unsigned int  read_old_parameter_file(FILE *ifile, int skip_header);
PRIVATE void          write_new_parameter_file(FILE *ofile, unsigned int options);
PRIVATE void          rd_stacks(int stack[NBPAIRS+1][NBPAIRS+1], FILE *fp);
PRIVATE void          rd_loop(int looparray[31], FILE *fp);
PRIVATE void          rd_mismatch(int mismatch[NBPAIRS+1][5][5], FILE *fp);
PRIVATE void          rd_int11(int int11[NBPAIRS+1][NBPAIRS+1][5][5], FILE *fp);
PRIVATE void          rd_int21(int int21[NBPAIRS+1][NBPAIRS+1][5][5][5], FILE *fp);
PRIVATE void          rd_int22(int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5], FILE *fp);
PRIVATE void          rd_dangle(int dangles[NBPAIRS+1][5], FILE *fp);
PRIVATE void          rd_MLparams(FILE *fp);
PRIVATE void          rd_misc(FILE *fp);
PRIVATE void          rd_ninio(FILE *fp);
PRIVATE void          rd_Tetra_loop(FILE *fp);
PRIVATE void          rd_Tri_loop(FILE *fp);
PRIVATE void          check_symmetry(void);
PRIVATE enum          parset_184 gettype_184(char ident[]);
PRIVATE char          *get_array1(int *arr, int size, FILE *fp);
PRIVATE void          ignore_comment(char *line);
PRIVATE void          display_array(int *p, int size, int line, FILE *fp);


PUBLIC void convert_parameter_file(const char *iname, const char *oname, unsigned int options){
  FILE          *ifile, *ofile;
  unsigned int  old_options = 0;
  int           skip_input_header = 0;

  if(options & VRNA_CONVERT_OUTPUT_DUMP){
    if(oname == NULL) oname = iname;
    skip_input_header = 1;
  }
  else{
    if(iname == NULL){
      ifile = stdin;
      skip_input_header = 1;
    }
    else if(!(ifile=fopen(iname,"r"))){
      vrna_message_warning("convert_epars: can't open file %s", iname);
      return;
    }
    /* read old (1.8.4 format) parameter file */
    old_options = read_old_parameter_file(ifile, skip_input_header);
    if(ifile != stdin) fclose(ifile);
    check_symmetry();
  }

  if(options & VRNA_CONVERT_OUTPUT_VANILLA)
    options = old_options;

  if(oname == NULL) ofile = stdout;
  else if(!(ofile=fopen(oname,"a+"))){
    vrna_message_warning("convert_epars: can't open file %s for writing", oname);
    return;
  }
  write_new_parameter_file(ofile, options);
  if(ofile != stdout) fclose(ofile);
}


/*------------------------------------------------------------*/
PRIVATE unsigned int read_old_parameter_file(FILE *ifile, int skip_header){
  char                  *line, ident[32];
  enum      parset_184  type;
  int                   r, last;
  unsigned  int         read_successfully = 0;

  if (!(line = vrna_read_line(ifile))) {
    vrna_message_warning("convert_epars: can't read input parameter file");
    return 0;
  }
  if(!skip_header){
    if (strncmp(line,"## RNAfold parameter file",25)!=0){
      vrna_message_warning("convert_epars: Missing header line in input parameter file.\n"
                "May be this file has incorrect format?");
      free(line);
      return 0;
    }
    free(line);
    line = vrna_read_line(ifile);
  }
  last = 0;
  do{
    r = sscanf(line, "# %31s", ident);
    if (r==1){
      type = gettype_184(ident);
      switch (type){
        case QUIT_184:    if(ifile == stdin){
                            vrna_message_info(stderr, "press ENTER to continue...");
                            fflush(stderr);
                          }
                          last = 1;
                          break;
        case SH_184:      rd_stacks(enthalpies_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_STACK;
                          break;
        case S_184:       rd_stacks(stack37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_STACK;
                          break;
        case HP_184:      rd_loop(hairpin37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_HP;
                          break;
        case B_184:       rd_loop(bulge37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_BULGE;
                          break;
        case IL_184:      rd_loop(internal_loop37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT;
                          break;
        case MMH_184:     rd_mismatch(mismatchH37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_MM_HP;
                          break;
        case MMI_184:     rd_mismatch(mismatchI37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_MM_INT
                                              |VRNA_CONVERT_OUTPUT_MM_INT_1N  /* since 1:n-interior loop mismatches are treated seperately in 2.0 */
                                              |VRNA_CONVERT_OUTPUT_MM_INT_23; /* since 2:3-interior loop mismatches are treated seperately in 2.0 */
                          break;
        case MMM_184:     rd_mismatch(mismatchM37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_MM_MULTI;
                          break;
        case MM_H_184:    rd_mismatch(mism_H_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_MM_HP      /* since hairpin mismatches are treated seperately in 2.0 */
                                              |VRNA_CONVERT_OUTPUT_MM_INT     /* since interior loop  mismatches are treated seperately in 2.0 */
                                              |VRNA_CONVERT_OUTPUT_MM_INT_1N  /* since 1:n-interior loop mismatches are treated seperately in 2.0 */
                                              |VRNA_CONVERT_OUTPUT_MM_INT_23  /* since 2:3-interior loop mismatches are treated seperately in 2.0 */
                                              |VRNA_CONVERT_OUTPUT_MM_MULTI;  /* since multi loop mismatches are treated seperately in 2.0 */
                          break;
        case INT11_184:   rd_int11(int11_37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT_11;
                          break;
        case INT11_H_184: rd_int11(int11_H_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT_11;
                          break;
        case INT21_184:   rd_int21(int21_37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT_21;
                          break;
        case INT21_H_184: rd_int21(int21_H_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT_21;
                          break;
        case INT22_184:   rd_int22(int22_37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT_22;
                          break;
        case INT22_H_184: rd_int22(int22_H_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_INT_22;
                          break;
        case DE5_184:     rd_dangle(dangle5_37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_DANGLE5
                                              |VRNA_CONVERT_OUTPUT_MM_MULTI /* since multi loop mismatches were treated as dangle contribution */
                                              |VRNA_CONVERT_OUTPUT_MM_EXT;  /* since exterior loop mismatches were treated as dangle contribution */
                          break;
        case DE5_H_184:   rd_dangle(dangle5_H_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_DANGLE5
                                              |VRNA_CONVERT_OUTPUT_MM_MULTI /* since multi loop mismatches were treated as dangle contribution */
                                              |VRNA_CONVERT_OUTPUT_MM_EXT;  /* since exterior loop mismatches were treated as dangle contribution */
                          break;
        case DE3_184:     rd_dangle(dangle3_37_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_DANGLE3
                                              |VRNA_CONVERT_OUTPUT_MM_MULTI /* since multi loop mismatches were treated as dangle contribution */
                                              |VRNA_CONVERT_OUTPUT_MM_EXT;  /* since exterior loop mismatches were treated as dangle contribution */
                          break;
        case DE3_H_184:   rd_dangle(dangle3_H_184, ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_DANGLE3
                                              |VRNA_CONVERT_OUTPUT_MM_MULTI /* since multi loop mismatches were treated as dangle contribution */
                                              |VRNA_CONVERT_OUTPUT_MM_EXT;  /* since exterior loop mismatches were treated as dangle contribution */
                          break;
        case ML_184:      rd_MLparams(ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_ML
                                              |VRNA_CONVERT_OUTPUT_MISC;    /* since TerminalAU went to "misc" section */
                          break;
        case NIN_184:     rd_ninio(ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_NINIO;
                          break;
        case TL_184:      rd_Tetra_loop(ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_SPECIAL_HP;
                          break;
        case TRI_184:     rd_Tri_loop(ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_SPECIAL_HP;
                          break;
        case MISC_184:    rd_misc(ifile);
                          read_successfully |= VRNA_CONVERT_OUTPUT_MISC;
                          break;
        default:          /* do nothing but complain */
                          vrna_message_warning("convert_parameter_file: Unknown field identifier in `%s'", line);
      }
    } /* else ignore line */
    free(line);
  } while((line=vrna_read_line(ifile)) && !last);
  return read_successfully;
}

PRIVATE void display_array(int *p, int size, int nl, FILE *fp){
  int i;
  for (i=1; i<=size; i++, p++) {
    switch(*p){
      case  INF: fprintf(fp,"   INF");    break;
      case -INF: fprintf(fp,"  -INf");    break;
      case  DEF: fprintf(fp,"   DEF");    break;
      default:   fprintf(fp,"%6d",  *p);  break;
    }
    if ((i%nl)==0) fprintf(fp,"\n");
  }
  if (size%nl) fprintf(fp,"\n");
  return;
}

PRIVATE char *get_array1(int *arr, int size, FILE *fp){
  int    i, p, pos, pp, r, last;
  char  *line, buf[16];
  i = last = 0;
  while( i<size ) {
    line = vrna_read_line(fp);
    if (!line) vrna_message_error("convert_epars: unexpected end of file in get_array1");
    ignore_comment(line);
    pos=0;
    while ((i<size)&&(sscanf(line+pos,"%15s%n", buf, &pp)==1)) {
      pos += pp;
      if (buf[0]=='*') {i++; continue;}
      else if (buf[0]=='x') { /* should only be used for loop parameters */
        if (i==0) vrna_message_error("convert_epars: can't extrapolate first value");
        p = arr[last] + (int) (0.5+ lxc37_184*log(((double) i)/(double)(last)));
      }
      else if (strcmp(buf,"DEF") == 0) p = DEF;
      else if (strcmp(buf,"INF") == 0) p = INF;
      else if (strcmp(buf,"NST") == 0) p = NST;
      else {
        r=sscanf(buf,"%d", &p);
        if (r!=1) {
          return line+pos;
          vrna_message_error("convert_epars: can't interpret `%s' in get_array1", buf);
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
/*------------------------------------------------------------*/

PRIVATE void  rd_stacks(int stacks[NBPAIRS+1][NBPAIRS+1], FILE *fp)
{
  int    i;
  char  *cp;
  for (i=1; i<=NBPAIRS; i++) {
    cp = get_array1(stacks[i]+1,NBPAIRS, fp);
    if (cp) {
      vrna_message_error("convert_epars: \nrd_stacks: %s", cp);
      exit(1);
    }
  }
  return;
}
/*------------------------------------------------------------*/

PRIVATE void rd_loop(int loop[31], FILE *fp)
{
  char *cp;

  cp   = get_array1(loop, 31, fp);

  if (cp) {
    vrna_message_error("convert_epars: \nrd_loop: %s", cp);
    exit(1);
  }
  return;
}
/*------------------------------------------------------------*/

PRIVATE void rd_mismatch(int mismatch[NBPAIRS+1][5][5], FILE *fp)
{
  char  *cp;
  int    i;

  for (i=1; i<NBPAIRS+1; i++) {
    cp = get_array1(mismatch[i][0],5*5, fp);
    if (cp) {
      vrna_message_error("convert_epars: rd_mismatch: in field mismatch[%d]\n\t%s", i, cp);
      exit(1);
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_int11(int int11[NBPAIRS+1][NBPAIRS+1][5][5], FILE *fp)
{
  char  *cp;
  int    i, j;

  for (i=1; i<NBPAIRS+1; i++) {
    for (j=1; j<NBPAIRS+1; j++) {
      cp = get_array1(int11[i][j][0],5*5, fp);
      if (cp) {
        vrna_message_error("convert_epars: rd_int11: in field int11[%d][%d]\n\t%s", i, j, cp);
        exit(1);
      }
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_int21(int int21[NBPAIRS+1][NBPAIRS+1][5][5][5], FILE *fp)
{
  char  *cp;
  int    i, j, k;

  for (i=1; i<NBPAIRS+1; i++) {
    for (j=1; j<NBPAIRS+1; j++) {
      for (k=0; k<5; k++) {
        cp = get_array1(int21[i][j][k][0],5*5, fp);
        if (cp) {
          vrna_message_error("convert_epars: rd_int21: in field int21[%d][%d][%d]\n\t%s",
                                    i, j, k, cp);
          exit(1);
        }
      }
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_int22(int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5], FILE *fp)
{
  char  *cp;
  int    i, j, k, l, m;

  for (i=1; i<NBPAIRS+1; i++)
    for (j=1; j<NBPAIRS+1; j++)
      for (k=1; k<5; k++)
        for (l=1; l<5; l++)
          for (m=1; m<5; m++) {
            cp = get_array1(int22[i][j][k][l][m]+1,4, fp);
            if (cp) {
              vrna_message_error("convert_epars: rd_int22: in field "
                                        "int22[%d][%d][%d][%d][%d]\n\t%s",
                                        i, j, k, l, m, cp);
              exit(1);
            }
          }

  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_dangle(int dangle[NBPAIRS+1][5], FILE *fp)
{
  int   i;
  char *cp;

  for (i=0; i< NBPAIRS+1; i++) {
    cp = get_array1(dangle[i],5, fp);
    if (cp) {
      vrna_message_error("convert_epars: \nrd_dangle: %s", cp);
      exit(1);
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_MLparams(FILE *fp)
{
  char *cp;
  int values[4];

  cp   = get_array1(values,4, fp);
  if (cp) {
    vrna_message_error("convert_epars: rd_MLparams: %s", cp);
    exit(1);
  }

  ML_BASE37_184     = values[0];
  ML_closing37_184  = values[1];
  ML_intern37_184   = values[2];
  TerminalAU_184    = values[3];

  return;
}

/*------------------------------------------------------------*/

PRIVATE void  rd_misc(FILE *fp)
{
  char *cp;
  int values[1]; /* so far just one */

  cp   = get_array1(values,1, fp);
  if (cp) {
    vrna_message_error("convert_epars: rd_misc: %s", cp);
    exit(1);
  }

  DuplexInit_184 = values[0];

  return;
}

/*------------------------------------------------------------*/

PRIVATE void  rd_ninio(FILE *fp)
{
  char  *cp;
  int temp[2];

  cp = get_array1(temp, 2, fp);

  if (cp) {
    vrna_message_error("convert_epars: rd_F_ninio: %s", cp);
    exit(1);
  }
  F_ninio37_184[2] = temp[0];
  MAX_NINIO_184  = temp[1];
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_Tetra_loop(FILE *fp)
{
  int    i, r;
  char   *buf;

  i=0;
  memset(&Tetraloops_184, 0, 1400);
  memset(&TETRA_ENERGY37_184, 0, sizeof(int)*200);
  do {
    buf = vrna_read_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%6s %d", &Tetraloops_184[7*i], &TETRA_ENERGY37_184[i]);
    strcat(Tetraloops_184, " ");
    free(buf);
    i++;
  } while((r==2)&&(i<200));
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_Tri_loop(FILE *fp)
{
  int    i, r;
  char   *buf;

  i=0;
  memset(&Triloops_184, 0, 241);
  memset(&Triloop_E37_184, 0, sizeof(int)*40);
  do {
    buf = vrna_read_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%5s %d", &Triloops_184[6*i], &Triloop_E37_184[i]);
    Triloops_184[6*i+5]=' ';
    free(buf);
    i++;
  } while((r==2)&&(i<40));
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
      vrna_message_error("convert_epars: unclosed comment in parameter file");
    /* can't use strcpy for overlapping strings */
    for (cp2+=2; *cp2!='\0'; cp2++, cp1++)
      *cp1 = *cp2;
    *cp1 = '\0';
  }

  return;
}

PRIVATE enum parset_184 gettype_184(char ident[]){
  if (strcmp(ident,"stack_enthalpies") == 0)          return SH_184;
  else if (strcmp(ident,"stack_energies") == 0)       return S_184;
  else if (strcmp(ident,"hairpin" ) == 0)             return HP_184;
  else if (strcmp(ident,"bulge") == 0)                return B_184;
  else if (strcmp(ident,"internal_loop") == 0)        return IL_184;
  else if (strcmp(ident,"mismatch_hairpin") == 0)     return MMH_184;
  else if (strcmp(ident,"mismatch_interior") == 0)    return MMI_184;
  else if (strcmp(ident,"mismatch_multi") == 0)       return MMM_184;
  else if (strcmp(ident,"mismatch_enthalpies") == 0)  return MM_H_184;
  else if (strcmp(ident,"int11_energies") == 0)       return INT11_184;
  else if (strcmp(ident,"int11_enthalpies") == 0)     return INT11_H_184;
  else if (strcmp(ident,"int21_energies") == 0)       return INT21_184;
  else if (strcmp(ident,"int21_enthalpies") == 0)     return INT21_H_184;
  else if (strcmp(ident,"int22_energies") == 0)       return INT22_184;
  else if (strcmp(ident,"int22_enthalpies") == 0)     return INT22_H_184;
  else if (strcmp(ident,"dangle5")== 0)               return DE5_184;
  else if (strcmp(ident,"dangle3")== 0)               return DE3_184;
  else if (strcmp(ident,"dangle5_enthalpies")== 0)    return DE5_H_184;
  else if (strcmp(ident,"dangle3_enthalpies")== 0)    return DE3_H_184;
  else if (strcmp(ident,"ML_params")== 0)             return ML_184;
  else if (strcmp(ident,"NINIO") == 0)                return NIN_184;
  else if (strcmp(ident,"Tetraloops") == 0)           return TL_184;
  else if (strcmp(ident,"Triloops") == 0)             return TRI_184;
  else if (strcmp(ident,"END") == 0)                  return QUIT_184;
  else return UNKNOWN_184;
}

PRIVATE void write_new_parameter_file(FILE *ofile, unsigned int option_bits){
  int           c;
  char          *pnames[] = {"NP", "CG", "GC", "GU", "UG", "AU", "UA", " @"};
  char          bnames[]  = "@ACGU";
  unsigned  int options   = 0;

  options = (option_bits & VRNA_CONVERT_OUTPUT_ALL) ?
              VRNA_CONVERT_OUTPUT_HP
            | VRNA_CONVERT_OUTPUT_STACK
            | VRNA_CONVERT_OUTPUT_MM_HP
            | VRNA_CONVERT_OUTPUT_MM_INT
            | VRNA_CONVERT_OUTPUT_MM_INT_1N
            | VRNA_CONVERT_OUTPUT_MM_INT_23
            | VRNA_CONVERT_OUTPUT_MM_MULTI
            | VRNA_CONVERT_OUTPUT_MM_EXT
            | VRNA_CONVERT_OUTPUT_DANGLE5
            | VRNA_CONVERT_OUTPUT_DANGLE3
            | VRNA_CONVERT_OUTPUT_INT_11
            | VRNA_CONVERT_OUTPUT_INT_21
            | VRNA_CONVERT_OUTPUT_INT_22
            | VRNA_CONVERT_OUTPUT_BULGE
            | VRNA_CONVERT_OUTPUT_INT
            | VRNA_CONVERT_OUTPUT_ML
            | VRNA_CONVERT_OUTPUT_MISC
            | VRNA_CONVERT_OUTPUT_SPECIAL_HP
            | VRNA_CONVERT_OUTPUT_NINIO
            :
              option_bits;

  make_pair_matrix(); /* needed for special loop energy contributions */

  fprintf(ofile,"## RNAfold parameter file v2.0\n");

  if(options & VRNA_CONVERT_OUTPUT_STACK){
    fprintf(ofile,"\n# %s\n", settype(S));
    fprintf(ofile,"/*  CG    GC    GU    UG    AU    UA    @  */\n");
    for (c=1; c<NBPAIRS+1; c++)
      display_array(stack37_184[c]+1,NBPAIRS,NBPAIRS, ofile);
    fprintf(ofile,"\n# %s\n", settype(S_H));
    fprintf(ofile,"/*  CG    GC    GU    UG    AU    UA    @  */\n");
    for (c=1; c<NBPAIRS+1; c++)
      display_array(enthalpies_184[c]+1,NBPAIRS,NBPAIRS, ofile);
  }

  if(options & VRNA_CONVERT_OUTPUT_MM_HP){
    fprintf(ofile,"\n# %s\n", settype(MMH));
    { int i,k;
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++)
          display_array(mismatchH37_184[k][i],5,5, ofile);
    }
    fprintf(ofile,"\n# %s\n", settype(MMH_H));
    { int i,k;
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++)
          display_array(mism_H_184[k][i],5,5, ofile);
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_MM_INT){
    fprintf(ofile,"\n# %s\n", settype(MMI));
    { int i,k;
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++)
          display_array(mismatchI37_184[k][i],5,5, ofile);
    }
    fprintf(ofile,"\n# %s\n", settype(MMI_H));
    { int i,k;
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++)
          display_array(mism_H_184[k][i],5,5, ofile);
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_MM_INT_1N){
    fprintf(ofile,"\n# %s\n", settype(MMI1N));
    { int i,k;
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++)
          display_array(mismatchI37_184[k][i],5,5, ofile);
    }
    fprintf(ofile,"\n# %s\n", settype(MMI1N_H));
    { int i,k;
    for (k=1; k<NBPAIRS+1; k++)
      for (i=0; i<5; i++)
        display_array(mism_H_184[k][i],5,5, ofile);
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_MM_INT_23){
    fprintf(ofile,"\n# %s\n", settype(MMI23));
    { int i,k;
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++)
          display_array(mismatchI37_184[k][i],5,5, ofile);
    }
    fprintf(ofile,"\n# %s\n", settype(MMI23_H));
    { int i,k;
    for (k=1; k<NBPAIRS+1; k++)
      for (i=0; i<5; i++)
        display_array(mism_H_184[k][i],5,5, ofile);
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_MM_MULTI){
    fprintf(ofile,"\n# %s\n", settype(MMM));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    { int i,j,k;
      int bla[5];
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++){
          for(j=0;j<5; j++)
            bla[j] = ((dangle5_37_184[k][i] == INF) ? 0 : dangle5_37_184[k][i]) + ((dangle3_37_184[k][j] == INF) ? 0 : dangle3_37_184[k][j]);
          display_array(bla,5,5, ofile);
        }
    }
    fprintf(ofile,"\n# %s\n", settype(MMM_H));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    { int i,j,k,bla[5];
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++){
          for(j=0;j<5; j++)
            bla[j] = ((dangle5_H_184[k][i] == INF) ? 0 : dangle5_H_184[k][i]) + ((dangle3_H_184[k][j] == INF) ? 0 : dangle3_H_184[k][j]);
          display_array(bla,5,5, ofile);
        }
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_MM_EXT){
    fprintf(ofile,"\n# %s\n", settype(MME));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    { int i,j,k;
      int bla[5];
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++){
          for(j=0;j<5; j++)
            bla[j] = ((dangle5_37_184[k][i] == INF) ? 0 : dangle5_37_184[k][i]) + ((dangle3_37_184[k][j] == INF) ? 0 : dangle3_37_184[k][j]);
          display_array(bla,5,5, ofile);
        }
    }
    fprintf(ofile,"\n# %s\n", settype(MME_H));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    { int i,j,k,bla[5];
      for (k=1; k<NBPAIRS+1; k++)
        for (i=0; i<5; i++){
          for(j=0;j<5; j++)
            bla[j] = ((dangle5_37_184[k][i] == INF) ? 0 : dangle5_H_184[k][i]) + ((dangle3_H_184[k][j] == INF) ? 0 : dangle3_H_184[k][j]);
          display_array(bla,5,5, ofile);
        }
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_DANGLE5){
    fprintf(ofile,"\n# %s\n", settype(D5));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    for (c=1; c<NBPAIRS+1; c++)
      display_array(dangle5_37_184[c], 5, 5, ofile);
    fprintf(ofile,"\n# %s\n", settype(D5_H));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    for (c=1; c<NBPAIRS+1; c++)
      display_array(dangle5_H_184[c], 5, 5, ofile);
  }

  if(options & VRNA_CONVERT_OUTPUT_DANGLE3){
    fprintf(ofile,"\n# %s\n", settype(D3));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    for (c=1; c<NBPAIRS+1; c++)
      display_array(dangle3_37_184[c], 5, 5, ofile);
    fprintf(ofile,"\n# %s\n", settype(D3_H));
    fprintf(ofile,"/*  @     A     C     G     U   */\n");
    for (c=1; c<NBPAIRS+1; c++)
      display_array(dangle3_H_184[c], 5, 5, ofile);
  }

  if(options & VRNA_CONVERT_OUTPUT_INT_11){
    /* don't print "no pair" entries for interior loop arrays */
    fprintf(ofile,"\n# %s\n", settype(INT11));
    { int i,k,l;
      for (k=1; k<NBPAIRS+1; k++)
        for (l=1; l<NBPAIRS+1; l++){
          fprintf(ofile, "/* %2s..%2s */\n", pnames[k], pnames[l]);
          for (i=0; i<5; i++)
            display_array(int11_37_184[k][l][i], 5, 5, ofile);
        }
    }
    fprintf(ofile,"\n# %s\n", settype(INT11_H));
    { int i,k,l;
      for (k=1; k<NBPAIRS+1; k++)
        for (l=1; l<NBPAIRS+1; l++){
          fprintf(ofile, "/* %2s..%2s */\n", pnames[k], pnames[l]);
          for (i=0; i<5; i++)
            display_array(int11_H_184[k][l][i],5,5, ofile);
        }
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_INT_21){
    fprintf(ofile,"\n# %s\n", settype(INT21));
    { int p1, p2, i, j;
      for (p1=1; p1<NBPAIRS+1; p1++)
        for (p2=1; p2<NBPAIRS+1; p2++)
          for (i=0; i<5; i++){
            fprintf(ofile, "/* %2s.%c..%2s */\n", pnames[p1], bnames[i], pnames[p2]);
            for (j=0; j<5; j++)
              display_array(int21_37_184[p1][p2][i][j],5,5, ofile);
          }
    }
    fprintf(ofile,"\n# %s\n", settype(INT21_H));
    { int p1, p2, i, j;
      for (p1=1; p1<NBPAIRS+1; p1++)
        for (p2=1; p2<NBPAIRS+1; p2++)
          for (i=0; i<5; i++){
            fprintf(ofile, "/* %2s.%c..%2s */\n", pnames[p1], bnames[i], pnames[p2]);
            for (j=0; j<5; j++)
              display_array(int21_H_184[p1][p2][i][j],5,5, ofile);
          }
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_INT_22){
    fprintf(ofile,"\n# %s\n", settype(INT22));
    { int p1, p2, i, j, k;
      for (p1=1; p1<NBPAIRS; p1++)
        for (p2=1; p2<NBPAIRS; p2++)
          for (i=1; i<5; i++)
            for (j=1; j<5; j++){
              fprintf(ofile, "/* %2s.%c%c..%2s */\n", pnames[p1], bnames[i], bnames[j], pnames[p2]);
              for (k=1; k<5; k++)
                display_array(int22_37_184[p1][p2][i][j][k]+1,4,5, ofile);
            }
    }
    fprintf(ofile,"\n# %s\n", settype(INT22_H));
    { int p1, p2, i, j, k;
      for (p1=1; p1<NBPAIRS; p1++)
        for (p2=1; p2<NBPAIRS; p2++)
          for (i=1; i<5; i++)
            for (j=1; j<5; j++){
              fprintf(ofile, "/* %2s.%c%c..%2s */\n", pnames[p1], bnames[i], bnames[j], pnames[p2]);
              for (k=1; k<5; k++)
                display_array(int22_H_184[p1][p2][i][j][k]+1,4,5, ofile);
            }
    }
  }

  if(options & VRNA_CONVERT_OUTPUT_HP){
    fprintf(ofile,"\n# %s\n", settype(HP));
    display_array(hairpin37_184, 31, 10, ofile);
    /* we had no hairpin enthalpies before, so
    *  we just pretend to have had some with value 0
    */
    fprintf(ofile,"\n# %s\n", settype(HP_H));
    {
      fprintf(ofile, "   INF   INF   INF");
      for(c=4;c<=31; c++){
        fprintf(ofile, "%6d", 0);
        if(c%10 == 0) fprintf(ofile, "\n");
      }
    }
    fprintf(ofile,"\n");
  }

  if(options & VRNA_CONVERT_OUTPUT_BULGE){
    fprintf(ofile,"\n# %s\n", settype(B));
    display_array(bulge37_184, 31, 10, ofile);

    /* we had no bulge enthalpies before, so
    *  we just pretend to have had some with value 0
    */
    fprintf(ofile,"\n# %s\n", settype(B_H));
    {
      fprintf(ofile, "   INF");
      for(c=2;c<=31; c++){
        fprintf(ofile, "%6d", 0);
        if(c%10 == 0) fprintf(ofile, "\n");
      }
    }
    fprintf(ofile,"\n");
  }

  if(options & VRNA_CONVERT_OUTPUT_INT){
    fprintf(ofile,"\n# %s\n", settype(IL));
    display_array(internal_loop37_184, 31, 10, ofile);

    /* we had no internal_loop enthalpies before, so
    *  we just pretend to have had some with value 0
    */
    fprintf(ofile,"\n# %s\n", settype(IL_H));
    {
      fprintf(ofile, "   INF   INF   INF   INF");
      for(c=5;c<=31; c++){
        fprintf(ofile, "%6d", 0);
        if(c%10 == 0) fprintf(ofile, "\n");
      }
    }
    fprintf(ofile,"\n");
    fprintf(ofile,"\n# %s\n"
                  "/* Ninio = MIN(max, m*|n1-n2| */\n"
                  "/*\t    m\t  m_dH     max  */\n"
                  "\t%6d\t%6d\t%6d\n", settype(NIN), F_ninio37_184[2], 0, MAX_NINIO_184);
  }

  if(options & VRNA_CONVERT_OUTPUT_ML){
    fprintf(ofile,"\n# %s\n", settype(ML));
    fprintf(ofile,"/* F = cu*n_unpaired + cc + ci*loop_degree (+TermAU) */\n");
    fprintf(ofile,"/*\t    cu\t cu_dH\t    cc\t cc_dH\t    ci\t ci_dH  */\n");
    fprintf(ofile,"\t%6d\t%6d\t%6d\t%6d\t%6d\t%6d\n", ML_BASE37_184, 0, ML_closing37_184, 0, ML_intern37_184, 0);
  }

  if(options & VRNA_CONVERT_OUTPUT_MISC){
    fprintf(ofile,"\n# %s\n", settype(MISC));
    fprintf(ofile,"/* all parameters are pairs of 'energy enthalpy' */\n");
    fprintf(ofile,"/*    DuplexInit     TerminalAU   LXC  */\n");
    fprintf(ofile,"   %6d %6d %6d %6d   %3.6f %6d\n", DuplexInit_184, 0, TerminalAU_184, 0, lxc37_184, 0);
  }

  if(options & VRNA_CONVERT_OUTPUT_SPECIAL_HP){
    fprintf(ofile,"\n# %s\n", settype(TRI));
    {
      int base_en = hairpin37_184[3];
      int base_dH = TETRA_ENTH37_184;
      for (c=0; c< (int)strlen(Triloops_184)/6; c++){
        int en = base_en;
        char bla[5];
        strncpy(bla, Triloops_184+c*6, 5);
        int type = pair[(short)encode_char(toupper(bla[0]))][(short)encode_char(toupper(bla[4]))];
        if(type > 2) en += TerminalAU_184;
        fprintf(ofile,"\t%.5s %6d %6d\n", Triloops_184+c*6, Triloop_E37_184[c] + en, base_dH);
      }
    }

    /* since the old hairpin loop function treated the tabulated tetraloop energy as bonus
    *  and the new one takes this tabulated energy as a total energy, we have to compute some
    *  things now...
    */
    fprintf(ofile,"\n# %s\n", settype(TL));
    {
      int base_en = hairpin37_184[4];
      int base_dH = TETRA_ENTH37_184;
      for (c=0; c< (int)strlen(Tetraloops_184)/7; c++){
        char bla[6];
        int en = base_en;
        int dH = base_dH;
        strncpy(bla, Tetraloops_184+c*7, 6);
        short si  = (short)encode_char(toupper(bla[1]));
        short sj  = (short)encode_char(toupper(bla[4]));
        int type  = pair[(short)encode_char(toupper(bla[0]))][(short)encode_char(toupper(bla[5]))];
        en   += mismatchH37_184[type][si][sj];
        dH   += mism_H_184[type][si][sj];
        fprintf(ofile,"\t%.6s %6d %6d\n", Tetraloops_184+c*7, en + TETRA_ENERGY37_184[c], dH);
      }
    }
    fprintf(ofile,"\n# %s\n", settype(HEX));
    {
      fprintf(ofile, "\n");
    }
  }

  fprintf(ofile, "\n# %s\n", settype(QUIT));
}

PRIVATE void check_symmetry(void) {
  int i,j,k,l;

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      if (stack37_184[i][j] != stack37_184[j][i])
        vrna_message_warning("stacking energies not symmetric");

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      if (enthalpies_184[i][j] != enthalpies_184[j][i])
        vrna_message_warning("stacking enthalpies not symmetric");


  /* interior 1x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++)
          if (int11_37_184[i][j][k][l] != int11_37_184[j][i][l][k])
            vrna_message_warning("int11 energies not symmetric");

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++)
          if (int11_H_184[i][j][k][l] != int11_H_184[j][i][l][k])
            vrna_message_warning("int11 enthalpies not symmetric");

  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)
              if (int22_37_184[i][j][k][l][m][n] != int22_37_184[j][i][m][n][k][l])
                vrna_message_warning("int22 energies not symmetric");
        }

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)
              if (int22_H_184[i][j][k][l][m][n] != int22_H_184[j][i][m][n][k][l])
                vrna_message_warning("int22 enthalpies not symmetric: %d %d %d %d %d %d", i,j,k,l,m,n);
        }
}
