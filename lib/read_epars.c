/*

   read_epars.c:
     read user defined energy parametrs from file
     with format that is produced by RNAedit

   (c) S.Kopp, TBI, Univ. of Vienna, Austria, Mar 1997
   email: robo@tbi.univie.ac.at


*/
     


#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#ifdef SGI
#include <malloc.h>
#endif

#include "utils.h"

#include "energy_const.h"
#include "energy_par.h"

/*------------------  identifiers ----------------------------*/
#define H        "enthalpies"
#define S        "entopies"
#define HP       "hairpin37"
#define B        "bulge37"
#define IL       "internal_loop37"
#define MM       "mismatch"
#define DE5      "dangle5"
#define DE3      "dangle3"
#define MLB      "ML_BASE37"
#define MLC      "ML_closing37"
#define MLI      "ML_intern37"
#define MN       "MAX_NINIO"
#define FN       "F_ninio37"
#define TE       "TETRA_ENERGY37"
#define TL       "Tetraloops"
#define END      "END"
/*------------------------------------------------------------*/

#define PUBLIC
#define PRIVATE   static

#define DEF_TEMP   37.0    /* default temperature */

/*----------------------- prototypes -------------------------*/
PUBLIC  void  read_energy_parameter(char fname[]);

PRIVATE void  rd_enthalpies(void);
PRIVATE void  rd_entropies(void);
PRIVATE void  rd_hairpin(void);
PRIVATE void  rd_bulge(void);
PRIVATE void  rd_internal_loop(void);
PRIVATE void  rd_mismatch(void);
PRIVATE void  rd_dangle5(void);
PRIVATE void  rd_dangle3(void);
PRIVATE void  rd_ML_base(void);
PRIVATE void  rd_ML_closing(void);
PRIVATE void  rd_ML_intern(void);
PRIVATE void  rd_MAX_ninio(void);
PRIVATE void  rd_F_ninio(void);
PRIVATE void  rd_Tetra_energy(void);
PRIVATE void  rd_Tetra_loop(void);
PRIVATE float rd_temp(FILE *f1);

PRIVATE char *get_array1(int *arr, int size, float temp);

PUBLIC  void  ignore_blanks(FILE *f1);
PUBLIC  void  ignore_comment(FILE *f1);
PUBLIC  void  read_to_eol(FILE *f1);

PRIVATE void  rescale(int *arr, char *ctrl, int size, float temp);
PRIVATE void  extrapolate(int *arr, char *to_be_x, int size);

PUBLIC  void  display_array(int *p, int size, int line);

/*------------------------------------------------------------*/
PRIVATE FILE *fp;

/*------------------------------------------------------------*/
PUBLIC void read_energy_parameter(char fname[])
{
  char    *line;
  int      c;
  
  if(!(fp=fopen(fname,"r"))){
    fprintf(stderr,
	    "\nread_energy_parameter:\n"
	    "\t\tcan't open file %s\n"
	    "\t\tusing default parameters instead.\n", fname);
    return;
  }

  line = get_line(fp);

  if(strcmp(line,"## RNAfold parameter file")!=0){
    fprintf(stderr,
	    "Missing header line in file.\n"
	    "May be this file has incorrect format.\n"
	    "Use INTERRUPT-key to stop.\n");
  }
  free(line);
  
  while((c=fgetc(fp))!=EOF){
    
    if(c=='#'){
      ignore_blanks(fp);
      line = get_line(fp);
  
      if(strcmp(line,END)==0)        break;
      else if(strcmp(line,H)  == 0) rd_enthalpies();    
      else if(strcmp(line,S)  == 0) rd_entropies();    
      else if(strcmp(line,HP) == 0) rd_hairpin();
      else if(strcmp(line,B)  == 0) rd_bulge(); 
      else if(strcmp(line,IL) == 0) rd_internal_loop();
      else if(strcmp(line,MM) == 0) rd_mismatch();
      else if(strcmp(line,DE5)== 0) rd_dangle5();
      else if(strcmp(line,DE3)== 0) rd_dangle3();    
      else if(strcmp(line,MLB)== 0) rd_ML_base();
      else if(strcmp(line,MLC)== 0) rd_ML_closing();
      else if(strcmp(line,MLI)== 0) rd_ML_intern();
      else if(strcmp(line,MN) == 0) rd_MAX_ninio();
      else if(strcmp(line,FN) == 0) rd_F_ninio();
      else if(strcmp(line,TE) == 0) rd_Tetra_energy();
      else if(strcmp(line,TL) == 0) rd_Tetra_loop();
      else
           fprintf(stderr," Unknown field identifier `%s'\n", line);

      free(line);  
    }
  }
  
  fclose(fp);
#ifdef CONTROL
#define WAIT    fprintf(stderr,"press enter key."); getchar()
  
  fprintf(stderr,"\n***********************************\n");

  fprintf(stderr,"Enthalpies:\n");
  for(c=0;c<NBPAIRS+1;c++)
    display_array(enthalpies[c],NBPAIRS+1,NBPAIRS+1);
  WAIT;
  
  fprintf(stderr,"\nEntropies:\n");
  for(c=0;c<NBPAIRS+1;c++)
    display_array(entropies[c],NBPAIRS+1,NBPAIRS+1);
  WAIT;
  
  fprintf(stderr,"\nHairpin:\n");
  display_array(hairpin37, 31, 10);
  WAIT;

  fprintf(stderr,"\nBulge:\n");
  display_array(bulge37, 31, 10);
  WAIT;
  
  fprintf(stderr,"\nInternal Loop:\n");
  display_array(internal_loop37, 31, 10);
  WAIT;
  
  fprintf(stderr,"\nTerminal Mismatch:\n");
  { int i,k,j;
  for(k=0;k<5;k++)
    for(j=0;j<5;j++){
      for(i=0;i<5;i++)
	display_array(mismatch[k][j][i],5,10);
      fprintf(stderr,"\n");
    }
  }
  WAIT;

  fprintf(stderr,"\nDangle5:\n");
  for(c=0;c<NBPAIRS+1; c++)
    display_array(dangle5[c], 5, 5);
  WAIT;

  fprintf(stderr,"\nDangle5:\n");
  for(c=0;c<NBPAIRS+1; c++)
    display_array(dangle3[c], 5, 5);
  WAIT;

  fprintf(stderr,"\nML_BASE37 = %d\n", ML_BASE37);
  fprintf(stderr,"ML_closing37 = ");
  display_array(ML_closing37,NBPAIRS+1,NBPAIRS+1);
  fprintf(stderr,"ML_intern37 = ");
  display_array(ML_intern37,NBPAIRS+1,NBPAIRS+1);
  WAIT;
  
  fprintf(stderr,"\nMAX_NINIO = %d\n", MAX_NINIO);
  fprintf(stderr,"F_ninio37 = ");
  display_array(F_ninio37,5, 5);
  WAIT;
  
  fprintf(stderr,"\nTetra_Energy37 = %d\n", TETRA_ENERGY37);
  fprintf(stderr,"Tetra Loops:\n");
  for(c=0; c< N_TETRALOOPS; c++)
    fprintf(stderr,"%s\n", Tetraloops[c]);
  WAIT;
  
		
#endif
  
  return;
}

/*------------------------------------------------------------*/
PUBLIC void display_array(int *p, int size, int nl)
{
  int i;
  
  for(i=1;i<=size;i++, p++){
    switch(*p)
      {
      case  INF: fprintf(stderr,"   INF"); break;
      case -INF: fprintf(stderr,"  -INf"); break;
      case  DEF: fprintf(stderr,"   DEF"); break;
      default:   fprintf(stderr,"%6d",  *p); break;
      }
    if(!(i%nl)) fprintf(stderr,"\n");
  }
  return;
}
/*------------------------------------------------------------*/
PRIVATE char *get_array1(int *arr, int size, float temp)
{
  int    c, d, i,
        *p, sign;
  char  *back,
        *ctrl,
        *extrapol,
         buf[16],
         X;

  back     = (char *)space(512);
  extrapol = (char *)space(size);
  ctrl     = (char *)space(size);

  memset(ctrl, 1, size*sizeof(char));

  p = arr;
  i = 0;
  while( i<size ){
    sign = 1;

    while((c=fgetc(fp)) == ' ');  /* skip blanks without ungetc!!! */
    if(c!='\n'){
      if(c =='-'){                 /* if neg val, eat one more char */      
	c=fgetc(fp);   
	sign = -1;
      }
      if(ungetc(c,fp)!=c){        /* pushback that char */
	  sprintf(back,"can't ungetc(%c)", c);
	  free(extrapol);
	  free(ctrl);
	  return back;
	}
      
           
      if(c=='*'){ ctrl[i] = 0; fgetc(fp);}  /* use default value */
      else if(c=='x'){
	X           = 1;
	ctrl[i]     = 0;  /* don't try to rescale this value */ 
	extrapol[i] = 1;
	fgetc(fp);
      }
      else if(isdigit(c))
	fscanf(fp,"%d", p);
      else {
	fscanf(fp,"%s", buf);
	if     (strcmp(buf,"DEF") == 0) *p = DEF;
	else if(strcmp(buf,"INF") == 0) *p = INF;
	else if(strcmp(buf,"NST") == 0) *p = NST;
	else if(strcmp(buf,"NSM") == 0) *p = NSM;
	else if(strcmp(buf,"/*")  == 0) /* */
	  read_to_eol(fp); 
	else {
	  sprintf(back,"unknown field entry '%s'", buf);
	  free(extrapol);
	  free(ctrl);
	  return back;
	}
      }
      *p++ *= sign;
      i++;
    }
  }

  if(temp<DEF_TEMP || temp > DEF_TEMP)
    rescale(arr,ctrl,size,temp);
  if(X) /* for bulge and interior loop only */
    extrapolate(arr,extrapol,size);
  
  back[0] = '\0';
  free(extrapol);
  free(ctrl);
  return back;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_enthalpies(void)
{
  int    i;
  char  *cp;
  float  temp;

  temp = rd_temp(fp);
  ignore_comment(fp);
  for(i=0; i<NBPAIRS+1;i++){
    cp = get_array1(enthalpies[i],NBPAIRS+1, temp);
    if(strlen(cp)){
       fprintf(stderr,"\nrd_enthalpies: %s\n", cp);
       exit(1);
    }
  free(cp);
  }
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_entropies(void)
{
  int   i;
  char *cp;
  float temp;

  temp = rd_temp(fp);
  ignore_comment(fp);
  for(i=0; i<NBPAIRS+1;i++){
    cp = get_array1(entropies[i],NBPAIRS+1, temp);
    if(strlen(cp)){
       fprintf(stderr,"\nrd_entropies: %s\n", cp);
       exit(1);
    }
    free(cp);
  }
  return;
}
/*------------------------------------------------------------*/
PRIVATE void rd_hairpin(void)
{
  char *cp;
  float temp;

  temp = rd_temp(fp);
  cp   = get_array1(hairpin37, 31, temp);
 
  if(strlen(cp)){
    fprintf(stderr,"\nrd_hairpin: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_bulge(void)
{
  char *cp;
  float temp;
  
  temp = rd_temp(fp);
  cp   = get_array1(bulge37, 31, temp);

  if(strlen(cp)){
    fprintf(stderr,"\nrd_bulge: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void rd_internal_loop(void)
{
  char *cp;
  float temp;

  temp = rd_temp(fp);
  
  cp = get_array1(internal_loop37, 31, temp);
  if(strlen(cp)){
    fprintf(stderr,"\nrd_internal_loop: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}

/*------------------------------------------------------------*/
PRIVATE void rd_mismatch(void)
{
  char  *cp;
  int    i, j, k,
         c;
  float  temp;

  temp = rd_temp(fp);
  ignore_comment(fp);          /* skip first line */

  for(i=0;i<5;i++)
    for(j=0;j<5;j++){
      while((c=fgetc(fp)) != '\n');
      if(c!='\n'){
	if(ungetc(c,fp)!=c){
	   fprintf(stderr,
		   "unpleasent error while unget(%c)\n", c);
	   exit(-1);
	 }
      }
      ignore_comment(fp);
      for(k=0;k<5;k++){
	cp = get_array1(mismatch[i][j][k],5,temp);
	if(strlen(cp)){
	  fprintf(stderr,
		  "rd_mismatch: in field mismatch[%d][%d][%d]\n"
		  "\t%s\n", i, j, k, cp);
	  exit(1);
	}
	free(cp);
      }
    }
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_dangle5(void)
{
  int   i;
  char *cp;
  float temp;

  temp = rd_temp(fp);
  ignore_comment(fp);
  for(i=0; i< NBPAIRS+1; i++){
    cp = get_array1(dangle5[i],5,temp);
    if(strlen(cp)){
      fprintf(stderr,"\nrd_dangle5: %s\n", cp);
      exit(1);
    }
    read_to_eol(fp);
    free(cp);
  }
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_dangle3(void)
{
  int   i;
  char *cp;
  float temp;

  temp = rd_temp(fp);
  ignore_comment(fp);
  for(i=0; i< NBPAIRS+1; i++){
    cp = get_array1(dangle3[i],5,temp);
    if(strlen(cp)){
      fprintf(stderr,"\nrd_dangle3: %s\n", cp);
      exit(1);
    }
    read_to_eol(fp);
    free(cp);
  }
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_ML_base(void)
{
  char *cp;
  float temp;

  temp = rd_temp(fp);
  cp   = get_array1(&ML_BASE37,1, temp);

  if(strlen(cp)){
    fprintf(stderr,"rd_ML_base: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_ML_closing(void)
{
  char  *cp;
  float temp;

  temp = rd_temp(fp);
  cp   = get_array1(ML_closing37,NBPAIRS+1,temp);

  if(strlen(cp)){
    fprintf(stderr,"rd_ML_closing: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_ML_intern(void)
{
  char  *cp;
  float  temp;

  temp = rd_temp(fp);
  cp   = get_array1(ML_intern37,NBPAIRS+1,temp);

  if(strlen(cp)){
    fprintf(stderr,"rd_ML_intern: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_MAX_ninio(void)
{
  char  *cp;
  
  cp = get_array1(&MAX_NINIO, 1, DEF_TEMP);

  if(strlen(cp)){
    fprintf(stderr,"rd_MAX_ninio: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_F_ninio(void)
{
  char  *cp;
  float  temp;

  temp = rd_temp(fp);
  cp = get_array1(F_ninio37, 5, temp);

  if(strlen(cp)){
    fprintf(stderr,"rd_F_ninio: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_Tetra_energy(void)
{
  char  *cp;
  float  temp;

  temp = rd_temp(fp);
  cp   = get_array1(&TETRA_ENERGY37,1, temp);

  if(strlen(cp)){
    fprintf(stderr,"rd_Tetra_energy: %s\n", cp);
    exit(1);
  }
  free(cp);
  return;
}
/*------------------------------------------------------------*/
PRIVATE void  rd_Tetra_loop(void)
{
  int    i, c;
  char   *buf;
  
  for(i=0; i<N_TETRALOOPS; i++ ){
    ignore_blanks(fp);
    buf = get_line(fp);
    if(strlen(buf)!=4){
      fprintf(stderr,"rd_Tetra_loop: wrong tetraloop %s\n", buf);
      exit(1);
    }
    if(strcmp(buf,"****")!= 0)
      strcpy(Tetraloops[i],buf);
    free(buf);
  }
  return;
}
/*------------------------------------------------------------*/
PRIVATE float rd_temp(FILE *f1)
{
  char  *line, c;
  float  temp;

  if((c=fgetc(f1)) != '#')
    nrerror("error in parameter file for getting temperature");
  line = get_line(f1);
  sscanf(line,"%f", &temp );
  free(line);

  return temp;
}
/*------------------------------------------------------------*/
PUBLIC void ignore_blanks(FILE *f1)
{
  int c;

  while(isspace((c=fgetc(f1))));
  if(ungetc(c,f1)!=c)
    nrerror("error ocurred while skipping blanks in paramter file");

  return;
}
/*------------------------------------------------------------*/  
PUBLIC void ignore_comment(FILE *f1)
{
  int c, b;

  while((c=fgetc(f1)) == ' ');  /* skip blanks */
  if(c=='/'){
    if((b=fgetc(f1)) == '*'){
      while((b=fgetc(f1))!= '/' && b != EOF && c != '*' ) c = b;
      return;
    }
    else if(b=='/') {
      while(fgetc(fp)!='\n');
      return; 
    }
    else{
      if(ungetc(b,f1)!=b)
	nrerror("error ocurred while skipping comment in paramter file");
    }
  }
  if(ungetc(c,f1)!=c)
    nrerror("error ocurred while skipping comment in paramter file");

  return;
}
/*------------------------------------------------------------*/  
PUBLIC void read_to_eol(FILE *f1)
{
  char *line;

  line = get_line(f1);
  free(line);
  return;
}
/*------------------------------------------------------------*/  
PRIVATE void rescale(int *arr, char *ctrl, int size, float temp)
     /* rescale parameters to 37 C */
{
  short  i;
  int   *p;
  char  *cp;
  float  tempf;

  tempf =  (K0+DEF_TEMP)/(K0+temp);
#ifdef CONTROL
  fprintf(stderr,
	  "\n tempf = (K0+DEF_TEMP)/(K0+temp) = (%g+%g)/(%g+%g) = %g\n",
	  K0,DEF_TEMP,K0,temp,tempf);
#endif
            /* inverse to tempf in scale_parameters!!!  */
  p  = arr;
  cp = ctrl;

 
  for(i=0; i<size; i++, p++, cp++){
    if(*cp){
      switch(*p)
	{
	case    0:
	case  INF:
	case -INF:  /* do nothing */
	  break;
	default:
	  *p  = (int)(*p *(tempf));
	  if(*p<-INF)      *p = -INF;
	  else if(*p>INF) *p =  INF;
	  break;
	}
    }
  }


  return;
}
/*------------------------------------------------------------*/
PRIVATE void extrapolate(int *arr, char *to_be_x, int size)
     /* extrapolate unknown values, from the last known one */
{
  short   i,
          firstx,
          last;
  int     base,
         *p;
  double  no;
  char   *cp,
          flag;

  p       =  arr;
  cp      =  to_be_x;
  flag    =  0;
  last    = -1;
  firstx  = -1;

#ifdef CONTROL
  fprintf(stderr," Extrapolation performed\n");
#endif
  
  /* check for "negative" extrapolation */
  for(i=0;i<size; i++, cp++, p++){
    if(*p== INF || *p == -INF) *cp = 'i';
    else{
      if(*cp && firstx<0)  firstx = i;
      if(!(*cp) && last<0) last = i;
      if(firstx >= 0 && last >= 0){
	p    = arr + firstx;
	cp   = to_be_x + firstx;
	base = *(arr+last);
	no   = (double) last;
	if(firstx < last ){ /* negative extrapolation */
	  for(i=firstx; i<last; i++, cp++, p++){
	    *cp = 0;
	    *p  = base + (int)(lxc37*log((double)(i)/no));
	         /* taken from scale_parameters */
	    if(*p<-INF)     *p = -INF;
	    else if(*p>INF) *p =  INF;
	  }
	  i = size;
	}
	else i = size;  /* continue with normal extrapolation */ 
	 
      }
    }
  }

  /* normal extrapolation */
  p    =  arr;
  cp   =  to_be_x;
  
  for(i=0; i<size; i++, p++, cp++){
    if(*cp==1){
      if     (base== INF) *p =  INF;
      else if(base==-INF) *p = -INF;
      else{
	*p = base+(int)(lxc37*log((double)(i)/no));
	     /* taken from scale_parameters */
	if(*p<-INF)     *p = -INF;
	else if(*p>INF) *p =  INF;
      }
    }
    else if(*cp!='i'){
      base = *p;
      no   = (double)i;
    }
  }

  return;
}
/*------------------------------------------------------------*/ 
