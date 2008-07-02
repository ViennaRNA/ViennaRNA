/* Last changed Time-stamp: <2008-07-02 17:10:04 berni> */
/*
		  Ineractive Access to folding Routines

		  c Ivo L Hofacker
		  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"
#include "PS_dot.h"
#include "energy_const.h"
#include "LPfold.h"

extern float Lfold(char *string, char *structure, int winsize);
extern void  read_parameter_file(const char fname[]);

/*@unused@*/
static char rcsid[] = "$Id: RNAplfold.c,v 1.10 2008/07/02 15:34:24 ivo Exp $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
	"....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);
PRIVATE void putout_pup(double *pup,int length, int winsize, char *name);
PRIVATE void putoutpU_G(double **pU,int length, int ulength, FILE *fp);
int unpaired;
PRIVATE void putoutphakim_u(double **pU,int length, int ulength, FILE *fp);
/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *string, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[80], ffname[100];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, length, l, sym,r;
  int   istty;
  int noconv=0;
  int winsize=70;
  int pairdist=0;
  float cutoff=0.01;
  int tempwin=0;
  int temppair=0;
  int tempunpaired=0;
  FILE *pUfp=NULL;
  FILE *spup=NULL;
  double **pup=NULL; /*prob of being unpaired, lengthwise*/
  int plexoutput=0;
  plist *dpp=NULL;
  plist *pl;
  int simply_putout=0;
  int openenergies=0;
  string=NULL;
  dangles=2;
  unpaired=0;
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-')
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'n':
	  if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	  if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
	  if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	  if ( strcmp(argv[i], "-nsp") ==0) {
	    if (i==argc-1) usage();
	    ns_bases = argv[++i];
	  }
	  if ( strcmp(argv[i], "-noconv")==0) noconv=1;
	  break;
	case '4':
	  tetra_loop=0;
	  break;
	case 'e':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &energy_set);
	  if (!r) usage();
	  break;
	case 'c':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%f", &cutoff);
	  if (!r) usage();
	  break;
/*	case 'C':
	  fold_constrained=1;
	  break;
	case 'S':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%lf", &sfact);
	  if (!r) usage();
	  break;
*/
	case 'd': dangles=0;
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &dangles);
	    if (r!=1) usage();
	    if (!dangles%2) usage(); /*1,3 not possible*/
	  }
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	case 'O':
	  openenergies=1;
	   break;
	case 'W':
	  if (i==argc-1) usage();
	  r=sscanf(argv[++i], "%d", &winsize);
	  if (r!=1) usage();
	  break;
	case 'u':
	  /*default = 31*/
	  if (i==argc-1) {
	    unpaired=31;
	    fprintf(stderr,"Computing unpaired probabilities with default parameters: unpaired length=%d\n", unpaired);
	  }
	  else {
	    r=sscanf(argv[i+1], "%d", &unpaired);
	    if (r!=1){
	      unpaired=31;
	      fprintf(stderr,"Computing unpaired probabilities with default parameters: unpaired length=%d\n", unpaired);
	    }
	    else {
	      i++;
	    }
	  }
	  if (unpaired<0) usage();
	
	  break;
	case 'L':
	  if (i==argc-1) usage();
	  r=sscanf(argv[++i], "%d", &pairdist);
	  if (r!=1) usage();
	  break;
	case 'o':
	  simply_putout=1;
	  break;
	case 'p':
	  plexoutput=1;
	  break;
	default: usage();
	}
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

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
  if ((openenergies) &&(unpaired==0)) unpaired=31;
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
   if (pairdist==0) pairdist=winsize;
  if (pairdist>winsize) {
    fprintf(stderr, "pairdist (-L %d) should be <= winsize (-W %d);"
	    "Setting pairdist=winsize\n",pairdist, winsize);
    pairdist=winsize;
  }
  do {				/* main loop: continue until end of file */
   if (istty) {
      printf("\nInput string (upper or lower case); @ to quit\n");
      printf("%s\n", scale);
    }
    fname[0]='\0';
    if ((line = get_line(stdin))==NULL) break;

    if (tempwin!=0) {
      winsize=tempwin;
      tempwin=0;
    }
    if (temppair!=0) {
      pairdist=temppair;
      temppair=0;
    }
    if (tempunpaired!=0){
      unpaired=tempunpaired;
      tempunpaired=0;
    }
    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
      if (*line=='>')
	(void) sscanf(line, ">%42s", fname);
      printf("%s\n", line);
      free(line);
      if ((line = get_line(stdin))==NULL) break;
    }

    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;

    string = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",string);
    free(line);
    length = (int) strlen(string);
   if (unpaired) {
      pup=(double **)space((length+1)*sizeof(double *));
      pup[0]=(double *)space(sizeof(double)); /*I only need entry 0*/
      pup[0][0]=unpaired;
    }

    structure = (char *) space((unsigned) length+1);

    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    if (istty)
      printf("length = %d\n", length);

    if (length>1000000) {

      if (!simply_putout) {
	printf("Switched to simple output mode!!!\n");
      }
      simply_putout=1;
      
      }
    /* initialize_fold(length); */
    update_fold_params();
    if (length<winsize) {
      fprintf(stderr, "WARN: window size %d larger than sequence length %d\n",
	      winsize, length);
      tempwin=winsize;
	winsize=length;
	if (pairdist>winsize) {	  
	  temppair=pairdist;
	  pairdist=winsize;
	}
	if (unpaired>winsize) {
	  tempunpaired=unpaired;
	  unpaired=winsize;
	}
    }
    if (length >= 5) {
      pf_scale = -1;
      pUfp=NULL;
      if (simply_putout) /**/ {
	char oname[50];
	if (unpaired!=0) {
	  char aname[50];
	  if (fname[0]!='\0') {
	    strcpy(aname, fname);
	    strcat(aname, "_lunp");
	  }
	  else strcpy(aname, "plfold_lunp");
	pUfp=fopen(aname,"w");
	}
	if (fname[0]!='\0') {
	  strcpy(oname, fname);
	  strcat(oname, "_basepairs");
	}
	else strcpy(oname, "plfold_basepairs");
	spup=fopen(oname,"w");
      }
     
      pl=pfl_fold(string, winsize, pairdist, cutoff, pup, &dpp, pUfp,spup);
      if (pUfp!=NULL) fclose(pUfp);
      if (spup!=NULL)  fclose(spup);
     if (!simply_putout) {
	if (fname[0]!='\0') {
	  strcpy(ffname, fname);
	  strcat(ffname, "_dp.ps");
	}
	else strcpy(ffname, "plfold_dp.ps");
	PS_dot_plot_turn(string, pl, ffname, pairdist);
	free(pl);
	if (unpaired){
	  char aname[50];
	  if(plexoutput) {
	    FILE *POP;
	    char oname[50];
	    if (fname[0]!='\0') {
	      strcpy(oname, fname);
	      strcat(oname, "_uplex");
	    } 
	    else strcpy(oname, "plfold_uplex");
	    POP=fopen(oname,"w");
 	    putoutphakim_u(pup,length, unpaired, POP);
	    fclose(POP);
	  }
	  if (openenergies) {
	    if (fname[0]!='\0') {
	      strcpy(aname, fname);
	      strcat(aname, "_openen");
	    }
	    else strcpy(aname, "plfold_openen");
	  }
	  else {
	    if (fname[0]!='\0') {
	      strcpy(aname, fname);
	      strcat(aname, "_lunp");
	    }
	    else strcpy(aname, "plfold_lunp");
	  }
	  pUfp=fopen(aname,"w");
	  putoutpU_prob(pup, length, unpaired, pUfp,openenergies);
	  fclose(pUfp);
	   
	}
	
      }
     else {
       free(pl);
       if (unpaired) {
	 free(pup[0]);
	 free(pup);
       }
      
     }
      if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    }
    free(string);
    free(structure);
  } while (1);
  return 0;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAplfold [-L span] [-W winsize] [-u size of unpaired region]\n"
	  "          [-T temp] [-4] [-o] [-d[0|1|2]] [-noGU] [-noCloseGU]\n"
	  "          [-noLP] [-P paramfile] [-nsp pairs] [-noconv] [-O]\n");
}

PRIVATE void putout_pup(double *pup,int length, int winsize, char *name) {
  int i;
  float factor;
  float tfact;
  FILE *FP;
  
  FP=fopen(name,"w");
  fprintf(FP,"&prob of being unpaired between i-%d and i\n",unpaired);
  fflush(NULL);
  for (i=unpaired; i<=length; i++) {
    factor=0.;
    if (i<winsize) {
      factor=1./(i-unpaired+1);
    }
    if (i>length-winsize+unpaired-1) {
      tfact=1./(length-i+1);
      if (tfact>factor) {
	factor=tfact;
      }
     
    }
    else {
      tfact=1./(winsize-unpaired+1);
      if (tfact>factor) {
	factor=tfact;
      }
    }
    fprintf(FP,"%d %.6f\n",i,pup[i]*factor);
  }
  fclose(FP);

}
PRIVATE void putoutpU_G(double **pU,int length, int ulength, FILE *fp) {
  /*put out unpaireds */
  int i,k;
  fprintf(fp,"#unpaired probabilities\n #i\tl=");
  for (i=1; i<=ulength; i++) {
    fprintf(fp,"%d\t", i);
  }
  fprintf(fp,"\n");
  for (k=1; k<=length; k++){
    fprintf(fp,"%d\t",k);
    for (i=1; i<=ulength; i++) {
      if (k+(i-1)<=length) fprintf(fp,"%.7g\t",pU[k+(i-1)][i]);
    }
    fprintf(fp,"\n");
    free(pU[k]);
  }
  free(pU[0]);
  free(pU);
  fflush(fp);
}
PRIVATE void putoutphakim_u(double **pU,int length, int ulength, FILE *fp) {
  /*put out Fopen in dekacalories per mol, and F(cond,open) also ind dekacal*/
  int k;
  
  float RT = (temperature+K0)*GASCONST;
  float p0;
  float pdep;
  int f0;
  int fdep;
  
  fprintf(fp,"#energy necessary to unpair as well as to unpair if i-1 is unpaired also, if i+1 is unpaired also in dekacal/mol\n");
  for (k=1; k<=length; k++){
    fprintf(fp,"%d\t",k);   
    p0=pU[k][1];
    f0=(int) -RT*log(p0)/10;
    fprintf(fp,"%d\t", f0);
    if (k>1) {
      pdep=pU[k][2]/pU[k-1][1];
      fdep=(int) -RT*log(pdep)/10;
      fprintf(fp,"%d\t",fdep);


    }
    else  fprintf(fp,"-0\t");
    if (k<length) {
      pdep=pU[k+1][2]/pU[k+1][1];
      fdep=(int) -RT*log(pdep)/10;
      fprintf(fp,"%d\t",fdep);
    }
    else  fprintf(fp,"-0\t");
    fprintf(fp,"\n");
  }
 
  fflush(fp);
}

