/*
		  read energy parameters from a file

		      Stephan Kopp, Ivo Hofacker
			  Vienna RNA Package
*/
     

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "energy_const.h"
#include "energy_par.h"

static char rcsid[] = "$Id: read_epars.c,v 1.10 2004/12/10 16:32:35 ivo Exp $";

#define PUBLIC
#define PRIVATE   static
#define PARSET 20
enum parset {UNKNOWN= -1, QUIT, S, SH, HP, B, IL, MMI, MMH, MMM, MM_H,
	     DE5, DE3, DE5_H, DE3_H, ML, TL, TRI, TE, NIN, MISC,
	     INT11, INT11_H, INT21, INT21_H, INT22, INT22_H, DUMP, HELP}; 


/*------------------  identifiers ----------------------------*/
#define DEF -50
#define NST 0

#define DEF_TEMP   37.0    /* default temperature */

/*----------------------- prototypes -------------------------*/
PUBLIC  void  read_parameter_file(const char fname[]);
PUBLIC  void  write_parameter_file(const char fname[]);
  
PRIVATE void  rd_stacks(int stack[NBPAIRS+1][NBPAIRS+1]);
PRIVATE void  rd_loop(int looparray[31]);
PRIVATE void  rd_mismatch(int mismatch[NBPAIRS+1][5][5]);
PRIVATE void  rd_int11(int int11[NBPAIRS+1][NBPAIRS+1][5][5]);
PRIVATE void  rd_int21(int int21[NBPAIRS+1][NBPAIRS+1][5][5][5]);
PRIVATE void  rd_int22(int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5]);
PRIVATE void  rd_dangle(int dangles[NBPAIRS+1][5]);
PRIVATE void  rd_MLparams(void);
PRIVATE void  rd_misc(void);
PRIVATE void  rd_ninio(void);
PRIVATE void  rd_Tetra_loop(void);
PRIVATE void  rd_Tri_loop(void);
PRIVATE void check_symmetry(void);

PRIVATE enum parset gettype(char ident[]);
PRIVATE char *get_array1(int *arr, int size);

PRIVATE void  ignore_comment(char *line);

PRIVATE  void  display_array(int *p, int size, int line, FILE *fp);

/*------------------------------------------------------------*/
PRIVATE FILE *fp;
PRIVATE float rtemp=DEF_TEMP;

/*------------------------------------------------------------*/
PUBLIC void read_parameter_file(const char fname[])
{
  char    *line, ident[32];
  enum parset type;
  int      r, changed=0;

  if (!(fp=fopen(fname,"r"))) {
    fprintf(stderr,
	    "\nread_parameter_file:\n"
	    "\t\tcan't open file %s\n"
	    "\t\tusing default parameters instead.\n", fname);
    return;
  }

  if (!(line = get_line(fp))) {
    fprintf(stderr," File %s is inproper.\n", fname);
    fclose(fp);
    return;
  }

  if (strncmp(line,"## RNAfold parameter file",25)!=0) {
    fprintf(stderr,
	    "Missing header line in file.\n"
	    "May be this file has incorrect format.\n"
	    "Use INTERRUPT-key to stop.\n");
  }
  free(line);
  
  while((line=get_line(fp))) {
    
    r = sscanf(line, "# %31s", ident);
    if (r==1) {
      type = gettype(ident);
      switch (type)
	{
	case QUIT: break;
	case SH:     rd_stacks(enthalpies);    changed |= SH; break;
	case S:      rd_stacks(stack37);       changed |= S;  break;
	case HP:     rd_loop(hairpin37);       changed |= HP; break;
	case B:      rd_loop(bulge37);         changed |= B;  break;
	case IL:     rd_loop(internal_loop37); changed |= IL; break;
	case MMH:    rd_mismatch(mismatchH37); changed |= MMH; break;
	case MMI:    rd_mismatch(mismatchI37); changed |= MMI; break;
	case MMM:    rd_mismatch(mismatchM37); changed |= MMM; break;
	case MM_H:   rd_mismatch(mism_H);      changed |= MM_H; break;
	case INT11:  rd_int11(int11_37);       changed |= INT11; break;
	case INT11_H:rd_int11(int11_H);        changed |= INT11_H; break;
	case INT21:  rd_int21(int21_37);       changed |= INT21; break;
	case INT21_H:rd_int21(int21_H);        changed |= INT21_H; break;
	case INT22:  rd_int22(int22_37);       changed |= INT22; break;
	case INT22_H:rd_int22(int22_H);        changed |= INT22_H; break;
	case DE5:    rd_dangle(dangle5_37);    changed |= DE5;  break;
	case DE5_H:  rd_dangle(dangle5_H);     changed |= DE5_H;  break;
	case DE3:    rd_dangle(dangle3_37);    changed |= DE3; break;
	case DE3_H:  rd_dangle(dangle3_H);     changed |= DE3_H; break;
	case ML:     rd_MLparams();	       changed |= ML;  break;
	case NIN:    rd_ninio();	       changed |= NIN; break;
	case TL:     rd_Tetra_loop();          changed |= TL; break;
	case TRI:    rd_Tri_loop();            changed |= TRI; break;
	case MISC:   rd_misc();                changed |= MISC; break;
	  
	default: /* maybe it's a temperature */
	  r=sscanf(ident, "%f", &rtemp);
	  if (r!=1) fprintf(stderr," Unknown field identifier in `%s'\n", line);
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
    line = get_line(fp);
    if (!line) nrerror("unexpected end of file in get_array1");
    ignore_comment(line);
    pos=0;
    while ((i<size)&&(sscanf(line+pos,"%15s%n", buf, &pp)==1)) {
      pos += pp;
      if (buf[0]=='*') {i++; continue;}
      else if (buf[0]=='x') { /* should only be used for loop parameters */
	if (i==0) nrerror("can't extrapolate first value");
	p = arr[last] + (int) (0.5+ lxc37*log(((double) i)/(double)(last)));
      }
      else if (strcmp(buf,"DEF") == 0) p = DEF;
      else if (strcmp(buf,"INF") == 0) p = INF;
      else if (strcmp(buf,"NST") == 0) p = NST;
      else {
	r=sscanf(buf,"%d", &p);
	if (r!=1) {
	  return line+pos;
	  fprintf(stderr, "can't interpret `%s' in get_array1\n", buf);
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

PRIVATE void  rd_stacks(int stacks[NBPAIRS+1][NBPAIRS+1])
{
  int    i;
  char  *cp; 
  for (i=1; i<=NBPAIRS; i++) {
    cp = get_array1(stacks[i]+1,NBPAIRS);
    if (cp) {
      fprintf(stderr,"\nrd_stacks: %s\n", cp);
      exit(1);
    }
  }
  return;
}
/*------------------------------------------------------------*/

PRIVATE void rd_loop(int loop[31])
{
  char *cp;
  
  cp   = get_array1(loop, 31);
  
  if (cp) {
    fprintf(stderr,"\nrd_loop: %s\n", cp);
    exit(1);
  }
  return;
}
/*------------------------------------------------------------*/

PRIVATE void rd_mismatch(int mismatch[NBPAIRS+1][5][5])
{
  char  *cp;
  int    i;

  for (i=1; i<NBPAIRS+1; i++) {
    
    cp = get_array1(mismatch[i][0],5*5);
    if (cp) {
      fprintf(stderr, "rd_mismatch: in field mismatch[%d]\n\t%s\n", i, cp);
      exit(1);
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_int11(int int11[NBPAIRS+1][NBPAIRS+1][5][5])
{
  char  *cp;
  int    i, j;

  for (i=1; i<NBPAIRS+1; i++) {
    for (j=1; j<NBPAIRS+1; j++) {
      
      cp = get_array1(int11[i][j][0],5*5);
      if (cp) {
	fprintf(stderr, "rd_int11: in field int11[%d][%d]\n\t%s\n", i, j, cp);
	exit(1);
      }
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_int21(int int21[NBPAIRS+1][NBPAIRS+1][5][5][5])
{
  char  *cp;
  int    i, j, k;
  
  for (i=1; i<NBPAIRS+1; i++) {
    for (j=1; j<NBPAIRS+1; j++) {
      for (k=0; k<5; k++) {
	cp = get_array1(int21[i][j][k][0],5*5);
	if (cp) {
	  fprintf(stderr, "rd_int21: in field int21[%d][%d][%d]\n\t%s\n",
		  i, j, k, cp);
	  exit(1);
	}
      }
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_int22(int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5])
{
  char  *cp;
  int    i, j, k, l, m;
  
  for (i=1; i<NBPAIRS+1; i++) 
    for (j=1; j<NBPAIRS+1; j++)
      for (k=1; k<5; k++) 
	for (l=1; l<5; l++)
	  for (m=1; m<5; m++) {
	    cp = get_array1(int22[i][j][k][l][m]+1,4);
	    if (cp) {
	      fprintf(stderr, "rd_int22: in field "
		      "int22[%d][%d][%d][%d][%d]\n\t%s\n",
		      i, j, k, l, m, cp);
	      exit(1);
	    }
	  }

  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_dangle(int dangle[NBPAIRS+1][5])
{
  int   i;
  char *cp;

  for (i=0; i< NBPAIRS+1; i++) {
    cp = get_array1(dangle[i],5);
    if (cp) {
      fprintf(stderr,"\nrd_dangle: %s\n", cp);
      exit(1);
    }
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_MLparams(void)
{
  char *cp;
  int values[4];

  cp   = get_array1(values,4);
  if (cp) {
    fprintf(stderr,"rd_MLparams: %s\n", cp);
    exit(1);
  }

  ML_BASE37 = values[0];
  ML_closing37 = values[1];
  ML_intern37  = values[2];
  TerminalAU   = values[3];
  
  return;
}

/*------------------------------------------------------------*/

PRIVATE void  rd_misc(void)
{
  char *cp;
  int values[1]; /* so far just one */

  cp   = get_array1(values,1);
  if (cp) {
    fprintf(stderr,"rd_misc: %s\n", cp);
    exit(1);
  }

  DuplexInit = values[0];
  
  return;
}

/*------------------------------------------------------------*/

PRIVATE void  rd_ninio(void)
{
  char  *cp;
  int temp[2];

  cp = get_array1(temp, 2);

  if (cp) {
    fprintf(stderr,"rd_F_ninio: %s\n", cp);
    exit(1);
  }
  F_ninio37[2] = temp[0];
  MAX_NINIO  = temp[1];
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_Tetra_loop(void)
{
  int    i, r;
  char   *buf;

  i=0;
  do {
    buf = get_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%6s %d", &Tetraloops[7*i], &TETRA_ENERGY37[i]);
    strcat(Tetraloops, " ");
    free(buf);
    i++;
  } while((r==2)&&(i<200));
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_Tri_loop(void)
{
  int    i, r;
  char   *buf;

  i=0;
  do {
    buf = get_line(fp);
    if (buf==NULL) break;
    r = sscanf(buf,"%5s %d", &Triloops[6*i], &Triloop_E37[i]);
    Triloops[6*i+5]=' ';
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
      nrerror("unclosed comment in parameter file");
    /* can't use strcpy for overlapping strings */
    for (cp2+=2; *cp2!='\0'; cp2++, cp1++)  
      *cp1 = *cp2;
    *cp1 = '\0';
  }
  
  return;
}
/*------------------------------------------------------------*/  

/*@unused@*/
PRIVATE char *settype(enum parset s)
{
  switch(s)
    {
    case   SH: return "stack_enthalpies";
    case    S: return "stack_energies";
    case   HP: return "hairpin";
    case    B: return "bulge";
    case   IL: return "internal_loop";
    case  MMH: return "mismatch_hairpin";
    case  MMI: return "mismatch_interior";
    case  MMM: return "mismatch_multi";
    case MM_H: return "mismatch_enthalpies";
    case  DE5: return "dangle5";
    case  DE3: return "dangle3";
    case  DE5_H: return "dangle5_enthalpies";
    case  DE3_H: return "dangle3_enthalpies";
    case INT11:  return " int11_energies";
    case INT11_H:  return " int11_enthalpies";  
    case INT21:  return " int21_energies";
    case INT21_H:  return " int21_enthalpies";  
    case INT22:  return " int22_energies";
    case INT22_H:  return " int22_enthalpies";  
    case   ML: return "ML_params";
    case  NIN: return "NINIO";
    case   TL: return "Tetraloops";
    case  TRI: return "Triloops";  
    case QUIT:
    case DUMP:
    case HELP: return "";
    default: fprintf(stderr,"^8723300-3338111\n"); exit(-1);
    }
  return "";
}
/*------------------------------------------------------------*/ 

PRIVATE enum parset gettype(char ident[])
{
  if (strcmp(ident,"stack_enthalpies") == 0)         return SH;   
  else if (strcmp(ident,"stack_energies") == 0)      return S;    
  else if (strcmp(ident,"hairpin" ) == 0) 	     return HP;   
  else if (strcmp(ident,"bulge") == 0) 		     return B;    
  else if (strcmp(ident,"internal_loop") == 0) 	     return IL;   
  else if (strcmp(ident,"mismatch_hairpin") == 0)    return MMH;  
  else if (strcmp(ident,"mismatch_interior") == 0)   return MMI;
  else if (strcmp(ident,"mismatch_multi") == 0)      return MMM;  
  else if (strcmp(ident,"mismatch_enthalpies") == 0) return MM_H;
  else if (strcmp(ident,"int11_energies") == 0)      return INT11;  
  else if (strcmp(ident,"int11_enthalpies") == 0)    return INT11_H;
  else if (strcmp(ident,"int21_energies") == 0)      return INT21;  
  else if (strcmp(ident,"int21_enthalpies") == 0)    return INT21_H;
  else if (strcmp(ident,"int22_energies") == 0)      return INT22;  
  else if (strcmp(ident,"int22_enthalpies") == 0)    return INT22_H;
  else if (strcmp(ident,"dangle5")== 0) 	     return DE5;  
  else if (strcmp(ident,"dangle3")== 0)		     return DE3;  
  else if (strcmp(ident,"dangle5_enthalpies")== 0)   return DE5_H;  
  else if (strcmp(ident,"dangle3_enthalpies")== 0)   return DE3_H;  
  else if (strcmp(ident,"ML_params")== 0)	     return ML;  
  else if (strcmp(ident,"NINIO") == 0)	             return NIN;   
  else if (strcmp(ident,"Tetraloops") == 0)	     return TL;
  else if (strcmp(ident,"Triloops") == 0)	     return TRI;     
  else if (strcmp(ident, "END") == 0) 		     return QUIT;
  else return UNKNOWN;
}

/*---------------------------------------------------------------*/

PUBLIC void write_parameter_file(const char fname[]) {
  FILE *outfp;
  int c;
  char *pnames[] = {"NP", "CG", "GC", "GU", "UG", "AU", "UA", " @"};
  char bnames[] = "@ACGU";
  outfp = fopen(fname, "w");
  if (!outfp) {
    fprintf(stderr, "can't open file %s\n", fname);
    exit(1);
  }
  fprintf(outfp,"## RNAfold parameter file\n");
  
  fprintf(outfp,"\n# stack_energies\n");
  fprintf(outfp,"/*  CG    GC    GU    UG    AU    UA    @  */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(stack37[c]+1,NBPAIRS,NBPAIRS, outfp);
  
  fprintf(outfp,"\n# stack_enthalpies\n");
  fprintf(outfp,"/*  CG    GC    GU    UG    AU    UA    @  */\n");
  for (c=1; c<NBPAIRS+1; c++)
    display_array(enthalpies[c]+1,NBPAIRS,NBPAIRS, outfp);
  
  fprintf(outfp,"\n# mismatch_hairpin\n");
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++) 
      display_array(mismatchH37[k][i],5,5, outfp);
  }
  
  fprintf(outfp,"\n# mismatch_interior\n");
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++) 
      display_array(mismatchI37[k][i],5,5, outfp);
  }

  fprintf(outfp,"\n# mismatch_multi\n");
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++) 
      display_array(mismatchM37[k][i],5,5, outfp);
  }
  
  fprintf(outfp,"\n# mismatch_enthalpies\n");
  { int i,k;
  for (k=1; k<NBPAIRS+1; k++)
    for (i=0; i<5; i++) 
      display_array(mism_H[k][i],5,5, outfp);
  }
  fprintf(outfp,"\n# dangle5\n");
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=0; c<NBPAIRS+1; c++)
    display_array(dangle5_37[c], 5, 5, outfp);
  
  fprintf(outfp,"\n# dangle3\n");
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=0; c<NBPAIRS+1; c++)
    display_array(dangle3_37[c], 5, 5, outfp);
  
  fprintf(outfp,"\n# dangle5_enthalpies\n");
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=0; c<NBPAIRS+1; c++)
    display_array(dangle5_H[c], 5, 5, outfp);
  
  fprintf(outfp,"\n# dangle3_enthalpies\n");
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for (c=0; c<NBPAIRS+1; c++)
    display_array(dangle3_H[c], 5, 5, outfp);


  /* don;t print "no pair" entries for interior loop arrays */
  fprintf(outfp,"\n# int11_energies\n");
  { int i,k,l;
  for (k=1; k<NBPAIRS+1; k++)
    for (l=1; l<NBPAIRS+1; l++) {
      fprintf(outfp, "/* %2s..%2s */\n", pnames[k], pnames[l]);
      for (i=0; i<5; i++)
	display_array(int11_37[k][l][i], 5, 5, outfp);
    }
  }

  fprintf(outfp,"\n# int11_enthalpies\n");
  { int i,k,l;
  for (k=1; k<NBPAIRS+1; k++)
    for (l=1; l<NBPAIRS+1; l++) {
      fprintf(outfp, "/* %2s..%2s */\n", pnames[k], pnames[l]);
      for (i=0; i<5; i++) 
	display_array(int11_H[k][l][i],5,5, outfp);
    }
  }

  fprintf(outfp,"\n# int21_energies\n");
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

  fprintf(outfp,"\n# int21_enthalpies\n");
  { int p1, p2, i, j;
  for (p1=1; p1<NBPAIRS+1; p1++) 
    for (p2=1; p2<NBPAIRS+1; p2++)
      for (i=0; i<5; i++) {
	fprintf(outfp, "/* %2s.%c..%2s */\n",
		pnames[p1], bnames[i], pnames[p2]);
	for (j=0; j<5; j++)
	  display_array(int21_H[p1][p2][i][j],5,5, outfp);
      }
  }

  fprintf(outfp,"\n# int22_energies\n");
  { int p1, p2, i, j, k;
  for (p1=1; p1<NBPAIRS+1; p1++) 
    for (p2=1; p2<NBPAIRS+1; p2++)
      for (i=1; i<5; i++)
	for (j=1; j<5; j++) {
	  fprintf(outfp, "/* %2s.%c%c..%2s */\n",
		  pnames[p1], bnames[i], bnames[j], pnames[p2]);
	  for (k=1; k<5; k++) 
	    display_array(int22_37[p1][p2][i][j][k]+1,4,5, outfp);
	}
  }
  
  fprintf(outfp,"\n# int22_enthalpies\n");
  { int p1, p2, i, j, k;
  for (p1=1; p1<NBPAIRS+1; p1++) 
    for (p2=1; p2<NBPAIRS+1; p2++)
      for (i=1; i<5; i++)
	for (j=1; j<5; j++) {
	  fprintf(outfp, "/* %2s.%c%c..%2s */\n",
		  pnames[p1], bnames[i], bnames[j], pnames[p2]);
	  for (k=1; k<5; k++) 
	    display_array(int22_H[p1][p2][i][j][k]+1,4,5, outfp);
	}
  }
  
  fprintf(outfp,"\n# hairpin\n");
  display_array(hairpin37, 31, 10, outfp);
  
  fprintf(outfp,"\n# bulge\n");
  display_array(bulge37, 31, 10, outfp);
  
  fprintf(outfp,"\n# internal_loop\n");
  display_array(internal_loop37, 31, 10, outfp);
  
  fprintf(outfp,"\n# ML_params\n");
  fprintf(outfp,"/* F = cu*n_unpaired + cc + ci*loop_degree (+TermAU) */\n");
  fprintf(outfp,"/*\t    cu\t    cc\t    ci\t TerminalAU */\n");
  fprintf(outfp,"\t%6d\t%6d\t%6d\t%6d\n",
	  ML_BASE37, ML_closing37, ML_intern37, TerminalAU);
  
  fprintf(outfp,"\n# NINIO\n"
	  "/* Ninio = MIN(max, m*|n1-n2| */\n"
	  "/*       m   max              */\n"
	  "\t%3d %4d\n", F_ninio37[2], MAX_NINIO);

  fprintf(outfp,"\n# Tetraloops\n");
  for (c=0; c< strlen(Tetraloops)/7; c++)
    fprintf(outfp,"\t%.6s\t%4d\n", Tetraloops+c*7, TETRA_ENERGY37[c]);

  fprintf(outfp,"\n# Triloops\n");
  for (c=0; c< strlen(Triloops)/6; c++)
    fprintf(outfp,"\t%.5s\t%4d\n", Triloops+c*6, Triloop_E37[c]);

  fprintf(outfp, "\n#END\n");
  fclose(outfp);
}

PRIVATE void check_symmetry(void) {
  int i,j,k,l;

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      if (stack37[i][j] != stack37[j][i])
	fprintf(stderr, "WARNING: stacking energies not symmetric\n");

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      if (enthalpies[i][j] != enthalpies[j][i])
	fprintf(stderr, "WARNING: stacking enthalpies not symmetric\n");

  
  /* interior 1x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) 
          if (int11_37[i][j][k][l] != int11_37[j][i][l][k])
	    fprintf(stderr, "WARNING: int11 energies not symmetric\n");

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) 
          if (int11_H[i][j][k][l] != int11_H[j][i][l][k])
	    fprintf(stderr, "WARNING: int11 enthalpies not symmetric\n");

  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)      
              if (int22_37[i][j][k][l][m][n] != int22_37[j][i][m][n][k][l])
		fprintf(stderr, "WARNING: int22 energies not symmetric\n");
        }

  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)      
              if (int22_H[i][j][k][l][m][n] != int22_H[j][i][m][n][k][l])
		fprintf(stderr, "WARNING: int22 enthalpies not symmetric: %d %d %d %d %d %d\n", i,j,k,l,m,n);
        }
}
