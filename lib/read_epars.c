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

#ifdef dmalloc
#include "/usr/local/include/dmalloc.h"
#define space(X) calloc(1,(X))
#endif

static char rcsid[] = "$Id: read_epars.c,v 1.4 1997/11/06 16:46:07 ivo Rel $";

#define PUBLIC
#define PRIVATE   static
#define PARSET 20
enum parset {UNKNOWN= -1, QUIT, S, SH, HP, B, IL, MMI, MMH, MM_H,
	     DE5, DE3, ML, TL, TE, FN, MN, DUMP, HELP};

PRIVATE char names[PARSET];

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
PRIVATE void  rd_dangle(int dangles[NBPAIRS+1][5]);
PRIVATE void  rd_MLparams(void);
PRIVATE void  rd_MAX_ninio(void);
PRIVATE void  rd_F_ninio(void);
PRIVATE void  rd_Tetra_loop(void);

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
  
  while(line=get_line(fp)) {
    
    r = sscanf(line, "# %31s", ident);
    if (r==1) {
      type = gettype(ident);
      switch (type)
	{
	case QUIT: break;
	case SH:   rd_stacks(enthalpies);    changed |= SH; break;
	case S:    rd_stacks(stack37);       changed |= S;  break;    
	case HP:   rd_loop(hairpin37);       changed |= HP; break;
	case B:    rd_loop(bulge37);         changed |= B;  break;
	case IL:   rd_loop(internal_loop37); changed |= IL; break;
	case MMH:  rd_mismatch(mismatchH37); changed |= MMH; break;
	case MMI:  rd_mismatch(mismatchI37); changed |= MMI; break;
	case MM_H: rd_mismatch(mism_H);      changed |= MM_H; break;
	case DE5:  rd_dangle(dangle5_37);    changed |= DE5;  break;  
	case DE3:  rd_dangle(dangle3_37);    changed |= DE3; break;  
	case ML:   rd_MLparams();	     changed |= ML;  break;  
	case MN:   rd_MAX_ninio();	     changed |= MN; break; 
	case FN:   rd_F_ninio();	     changed |= FN; break;
	case TL:   rd_Tetra_loop();          changed |= TL; break;
	default: /* maybe it's a temperature */
	  r=sscanf(ident, "%f", &rtemp);
	  if (r!=1) fprintf(stderr," Unknown field identifier in `%s'\n", line);
	}
    } /* else ignore line */
    free(line);  
  }
  
  fclose(fp);
  
  return;
}

/*------------------------------------------------------------*/

PRIVATE void display_array(int *p, int size, int nl, FILE *fp)
{
  int i;
  
  for(i=1;i<=size;i++, p++) {
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
	p = arr[last] + rint(lxc37*log(((double) i)/(double)(last)));
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

  for(i=0; i<NBPAIRS+1;i++) {
    cp = get_array1(stacks[i],NBPAIRS+1);
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

  for(i=0;i<NBPAIRS+1;i++) {
    
    cp = get_array1(mismatch[i][0],5*5);
    if (cp) {
      fprintf(stderr, "rd_mismatch: in field mismatch[%d]\n\t%s\n", i, cp);
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

  for(i=0; i< NBPAIRS+1; i++) {
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
  int values[3];

  cp   = get_array1(values,3);
  if (cp) {
    fprintf(stderr,"rd_MLparams: %s\n", cp);
    exit(1);
  }

  ML_BASE37 = values[0];
  ML_closing37 = values[1];
  ML_intern37  = values[2];
  
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_MAX_ninio(void)
{
  char  *cp;
  
  cp = get_array1(&MAX_NINIO, 1);

  if (cp) {
    fprintf(stderr,"rd_MAX_ninio: %s\n", cp);
    exit(1);
  }
  return;
}

/*------------------------------------------------------------*/
PRIVATE void  rd_F_ninio(void)
{
  char  *cp;

  cp = get_array1(F_ninio37, 5);

  if (cp) {
    fprintf(stderr,"rd_F_ninio: %s\n", cp);
    exit(1);
  }
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
    r = sscanf(buf,"%4s %d", Tetraloops+5*i, &TETRA_ENERGY37[i]);
    Tetraloops[5*i+4]=' ';
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
  
  if (cp1=strstr(line, "/*")) {
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
    case MM_H: return "mismatch_enthalpies";
    case  DE5: return "dangle5";
    case  DE3: return "dangle3";
    case   ML: return "ML_params";
    case   MN: return "MAX_NINIO";
    case   FN: return "F_ninio";
    case   TE: return "TETRA_ENERGY";
    case   TL: return "Tetraloops";
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
  else if (strcmp(ident,"mismatch_enthalpies") == 0) return MM_H; 
  else if (strcmp(ident,"dangle5")== 0) 	     return DE5;  
  else if (strcmp(ident,"dangle3")== 0)		     return DE3;  
  else if (strcmp(ident,"ML_params")== 0)	     return ML;  
  else if (strcmp(ident,"MAX_NINIO") == 0)	     return MN;   
  else if (strcmp(ident,"F_ninio") == 0) 	     return FN;   
  else if (strcmp(ident,"Tetraloops") == 0)	     return TL;   
  else if (strcmp(ident, "END") == 0) 		     return QUIT;
  else return UNKNOWN;
}

/*---------------------------------------------------------------*/

PUBLIC void write_parameter_file(const char fname[]) {
  FILE *outfp;
  int c;

  outfp = fopen(fname, "w");
  if (!outfp) {
    fprintf(stderr, "can't open file %s\n", fname);
    exit(1);
  }
  fprintf(outfp,"## RNAfold parameter file\n");
  
  fprintf(outfp,"\n# stack_energies\n");
  fprintf(outfp,"/*       CG    GC    GU    UG    AU    UA    @  */\n");
  for(c=0;c<NBPAIRS+1;c++)
    display_array(stack37[c],NBPAIRS+1,NBPAIRS+1, outfp);
  
  fprintf(outfp,"\n# stack_enthalpies\n");
  fprintf(outfp,"/*       CG    GC    GU    UG    AU    UA    @  */\n");
  for(c=0;c<NBPAIRS+1;c++)
    display_array(enthalpies[c],NBPAIRS+1,NBPAIRS+1, outfp);

  fprintf(outfp,"\n# hairpin\n");
  display_array(hairpin37, 31, 10, outfp);
  
  fprintf(outfp,"\n# bulge\n");
  display_array(bulge37, 31, 10, outfp);
  
  fprintf(outfp,"\n# internal_loop\n");
  display_array(internal_loop37, 31, 10, outfp);
  
  fprintf(outfp,"\n# mismatch_hairpin\n");
  { int i,k;
  for(k=0;k<NBPAIRS+1;k++)
    for(i=0;i<5;i++) 
      display_array(mismatchH37[k][i],5,5, outfp);
  }
  
  fprintf(outfp,"\n# mismatch_interior\n");
  { int i,k;
  for(k=0;k<NBPAIRS+1;k++)
    for(i=0;i<5;i++) 
      display_array(mismatchI37[k][i],5,5, outfp);
  }
  
  fprintf(outfp,"\n# mismatch_enthalpies\n");
  { int i,k;
  for(k=0;k<NBPAIRS+1;k++)
    for(i=0;i<5;i++) 
      display_array(mism_H[k][i],5,5, outfp);
  }
  
  fprintf(outfp,"\n# dangle5\n");
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for(c=0;c<NBPAIRS+1; c++)
    display_array(dangle5_37[c], 5, 5, outfp);
  
  fprintf(outfp,"\n# dangle3\n");
  fprintf(outfp,"/*  @     A     C     G     U   */\n");
  for(c=0;c<NBPAIRS+1; c++)
    display_array(dangle3_37[c], 5, 5, outfp);
  
  fprintf(outfp,"\n# ML_params\n");
  fprintf(outfp,"/* F = cu*n_unpaired + ci*loop_degree + cc  */\n");
  fprintf(outfp,"/*\t    cu\t    cc\t    ci */\n");
  fprintf(outfp,"\t%6d\t%6d\t%6d\n", ML_BASE37, ML_closing37, ML_intern37);
  
  fprintf(outfp,"\n# MAX_NINIO\n\t%d\n", MAX_NINIO);
  
  fprintf(outfp,"\n# Tetraloops\n");
  for(c=0; c< strlen(Tetraloops)/5; c++)
    fprintf(outfp,"\t%.4s\t%d\n", Tetraloops+c*5, TETRA_ENERGY37[c]);

  fclose(outfp);
}
