/* Last changed Time-stamp: <2005-02-16 14:58:00 ivo> */

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
#include "PS_dot.h"
#include "cofold.h"
#include "fold.h"
#include "co_part_func.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"

extern void  read_parameter_file(const char fname[]);

/*@unused@*/
static char rcsid[] = "$Id: RNAcofold.c,v 1.4 2005/02/16 13:59:19 ivo Exp $";

#define PRIVATE static

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))
PRIVATE char *costring(char *string);
PRIVATE char *tokenize(char *line);
PRIVATE void usage(void);
PRIVATE double do_partfunc(char *string, int length, int Switch, struct plist **tpr, struct plist **mf);
PRIVATE void free_franz(char *Astring, char *Bstring, plist *prAB, plist *prAA, plist *prBB, plist *prA, plist *prB, struct plist *mfAB, struct plist *mfAA, struct plist *mfBB, struct plist *mfA, struct plist *mfB);
PRIVATE int read_concentrationfile(char *fname, double *startc);
PRIVATE struct ConcEnt *do_the_concenratinger(char *Conc_file,double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconces);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *string, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[53], ffname[20], gfname[20];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  char *Concfile;
  int   i, length, l, sym, r;
  double energy, min_en;
  double kT, sfact=1.07;
  int   pf=0, istty;
  int noconv=0;
  int doT=0;    /*compute dimere free energies etc.*/
  int doC=0;    /*toggle to compute concentrations*/ 
  int doQ=0;    /*toggle to compute prob of base being paired*/
  int cofi=0;   /*toggle concentrations stdin / file*/
  double FEAB; /*free energy  of AB dimer*/
  double FEAA; /*free energy  of AA dimer*/
  double FEBB; /*free energy  of BB dimer*/
  double FEA;
  double FEB;
  struct plist *prAB;
  struct plist *prAA;   /*pair probabilities of AA dimer*/
  struct plist *prBB;
  struct plist *prA;
  struct plist *prB;
  struct plist *mfAB;
  struct plist *mfAA;   /*pair mfobabilities of AA dimer*/
  struct plist *mfBB;
  struct plist *mfA;
  struct plist *mfB;
  double *ConcAandB;
  
  do_backtrack = 1; 
  string=NULL;
  Concfile=NULL;
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') 
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'p':  pf=1;
	  if (argv[i][2]!='\0')
	    (void) sscanf(argv[i]+2, "%d", &do_backtrack);
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
	case 'C':
	  fold_constrained=1;
	  break;
	case 'S':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%lf", &sfact);
	  if (!r) usage();
	  break;
	case 'd': dangles=0;
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &dangles);
	    if (r!=1) usage();
	  }
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	case 'a':
	  doT=1;
	  pf=1;
	  break;
	case 'c':/*concentrations from stdin*/
	  doC=1; 
	  doT=1;
	  pf=1;
	  break;
	case 'f':/*concentrations in file*/
	  if (i==argc-1) usage();
	  Concfile = argv[++i];
	  doC=1;
	  cofi=1;
	  doT=1;
	  pf=1;
	  break;
	case 'q':
	  pf=1;
	  doQ=1;
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
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
  if ((fold_constrained)&&(istty)) {
    printf("Input constraints using the following notation:\n");
    printf("| : paired with another base\n");
    printf(". : no constraint at all\n");
    printf("x : base must not pair\n");
    printf("< : base i is paired with a base j<i\n");
    printf("> : base i is paired with a base j>i\n");
    printf("matching brackets ( ): base i pairs base j\n");
  } 
	
  do {				/* main loop: continue until end of file */
    cut_point = -1;
    if (doC) {
      ConcAandB=(double *)space(20*sizeof(double));
    }
    if (istty) {
      printf("\nInput string (upper or lower case); @ to quit\n");
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      printf("%s%s\n", scale1, scale2);
    }
    fname[0]='\0';
    if ((line = get_line(stdin))==NULL) break;
    
    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
      if (*line=='>')
	(void) sscanf(line, ">%51s", fname);
      printf("%s\n", line);
      free(line);
      if ((line = get_line(stdin))==NULL) break;
    } 

    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;

    string = tokenize(line); /* frees line */

    length = (int) strlen(string);
    if ((doC)&&!(cofi)) { /*read concentrations out of file*/
      i=0; /*kann i e nehmen??*/
      double tmp1, tmp2;
      printf("Please enter concentrations\n format: ConcA ConcB\n return to end\n");
      while ((line = get_line(stdin))!=NULL){
	if(*line=='\0') break;
	sscanf(line,"%lf %lf",&tmp1,&tmp2);
	ConcAandB[i++]=tmp1;
	ConcAandB[i++]=tmp2;
	
	if (i%20==0) {
	  ConcAandB=(double *)xrealloc(ConcAandB,(i+22)*sizeof(double));
	 
	}
	free(line);
      }
      ConcAandB[i++]=0;
      ConcAandB[i]=0;
    }
    
    structure = (char *) space((unsigned) length+1);
    if (fold_constrained) {
      cstruc = tokenize(get_line(stdin));
      if (cstruc!=NULL)
        strncpy(structure, cstruc, length);
      else
        fprintf(stderr, "constraints missing\n");
    }

    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    if (istty) {
      if (cut_point == -1)
	printf("length = %d\n", length);
      else
	printf("length1 = %d\nlength2 = %d\n",
	       cut_point-1, length-cut_point+1);
    }

    /*compute mfe of AB dimer*/
    min_en = cofold(string, structure);
    mfAB=(struct plist *) space(sizeof(struct plist) * (length+1));
    mfAB=get_mfe_plist(mfAB);
    
/*     for (l=0; mfAB[l].i!=0; l++) { */
/*       if (!SAME_STRAND(mfAB[l].i,mfAB[l].j)) { */
/* 	min_en += 4.1; /\*can i use DuplexInit there? if not: Ivooo??*\/ */
/* 	break; */
/*       } */
/*     } */
    if (cut_point == -1)    printf("%s\n%s", string, structure); /*no cofold*/
      
    else {
      char *pstring, *pstruct;
      pstring = costring(string);
      pstruct = costring(structure);
      printf("%s\n%s", pstring,  pstruct);
      free(pstring);
      free(pstruct);
    }
    if (istty) 
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);
    
      (void) fflush(stdout);
       
    if (fname[0]!='\0') {
      strcpy(ffname, fname);
      strcat(ffname, "_ss.ps");
      strcpy(gfname, fname);
      strcat(gfname, "_ss.g");
    } else {
      strcpy(ffname, "rna.ps");
      strcpy(gfname, "rna.g");
    }
    if (length<2000)
      (void) PS_rna_plot(string, structure, ffname);
    else { 
     fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");
     free_co_arrays();  
    }
    
    /*compute partition function*/
    if (pf) {
      
      if (dangles==1) {
	dangles=2;   /* recompute with dangles as in pf_fold() */
	min_en = energy_of_struct(string, structure);
	dangles=1;
      }
      	 
      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      init_co_pf_fold(length);

      if (cstruc!=NULL)
	strncpy(structure, cstruc, length+1);
      energy = co_pf_fold(string, structure);
      prAB=(struct plist *) space(sizeof(struct plist) * (2*length));
      prAB=get_plist(prAB, length,0.00001); 
      FEAB=energy;
      if (doQ) make_probsum(length,fname); /*compute prob of base paired*/
      free_co_arrays();
      if (doT) { /*cofold of all dimers, monomers*/
	int Blength, Alength;
	char  *Astring, *Bstring;
	char *Newstring;
	char *Newname;
	char *comment;
	struct ConcEnt *Conc;
	if (cut_point<0) {
	  printf("Sorry, i cannot do that with only one molecule, please give me two or leave it\n");
	  free(mfAB);
	  free(prAB);
	  continue;
	}
	if (dangles==1) dangles=2; /*merkmas??*/
	Alength=cut_point-1;    /*length of first molecule*/
	Blength=length-cut_point+1; /*length of 2nd molecule*/
	
	Astring=(char *)space(sizeof(char)*(Alength+1));/*Sequence of first molecule*/
	Bstring=(char *)space(sizeof(char)*(Blength+1));/*Sequence of second molecule*/
	strncat(Astring,string,Alength);      
	strncat(Bstring,string+Alength,Blength);
	
	
	
	/*compute AA dimer*/
	prAA=(struct plist *) space(sizeof(struct plist) * (4*Alength));
	mfAA=(struct plist *) space(sizeof(struct plist) * (Alength+1));
	FEAA=do_partfunc(Astring, Alength, 2, &prAA, &mfAA);
		/*compute BB dimer*/
	prBB=(struct plist *) space(sizeof(struct plist) * (4*Blength));
	mfBB=(struct plist *) space(sizeof(struct plist) * (Blength+1));
	FEBB=do_partfunc(Bstring, Blength, 2, &prBB, &mfBB);
	/*free_co_pf_arrays();*/
	
	/*compute A monomer*/
	prA=(struct plist *) space(sizeof(struct plist) * (2*Alength));
	mfA=(struct plist *) space(sizeof(struct plist) * (Alength+1));
	if (Alength>4) {  /*only if sec_struc is possible*/
	   FEA=do_partfunc(Astring, Alength, 1, &prA, &mfA);
	   /*	  free_pf_arrays();*/

	}
	else { /*no secondary structure*/
	  FEA=0.;
	  prA=(struct plist *)xrealloc(prA,sizeof(struct plist));
	  mfA=(struct plist *)xrealloc(mfA,sizeof(struct plist));
	  prA[0].i=mfA[0].i=0;
	  prA[0].j=mfA[0].j=0;
	  prA[0].p=mfA[0].p=0.;
	}
	/*compute B monomer*/
	prB=(struct plist *) space(sizeof(struct plist) * (2*Blength));
	mfB=(struct plist *) space(sizeof(struct plist) * (Blength+1));
	
	if (Blength>4) {
	  
	  FEB=do_partfunc(Bstring, Blength, 1, &prB, &mfB); 
	  /*	  free_pf_arrays();*/
	  
	  
	}
	else {
	  FEB=0.;
	  prB=(struct plist *)xrealloc(prB,sizeof(struct plist));
	  mfB=(struct plist *)xrealloc(mfB,sizeof(struct plist));
	  prB[0].i=mfB[0].i=0;
	  prB[0].j=mfB[0].j=0;
	  prB[0].p=mfB[0].p=0.;
	 
	} 

	compute_probabilities(&FEAB,&FEAA,&FEBB,&FEA,&FEB,
				      prAB,prAA,prBB,prA,prB,
				      Alength,Blength);
	printf("\nFree Energies:\nblubAB\t\tAA\t\tBB\t\tA\t\tB\n%.6f\t%6f\t%6f\t%6f\t%6f\n",FEAB,FEAA,FEBB,FEA,FEB);
	
	if (doC) {
	  Conc=do_the_concenratinger(Concfile,FEAB, FEAA, FEBB, FEA, FEB, ConcAandB);
	  free(Conc);/*freeen*/
	}
	
	if (fname[0]!='\0') {
	  strcpy(ffname, fname);
	  strcat(ffname, "_dp5.ps");
	} else strcpy(ffname, "dot5.ps");
	/*output of the 5 dot plots*/
	
	/*AB dot_plot*/
	/*write Free Energy into comment*/ 
	comment=(char *)space(80*sizeof(char));
	sprintf(comment,"\n%%FreeEnergy= %.9f\n",FEAB);
	/*reset cut_point*/
	cut_point=Alength+1;
	/*write New name*/
	Newname=(char*)space((strlen(ffname)+3)*sizeof(char));
	sprintf(Newname,"AB");
	strcat(Newname,ffname);
	(void)PS_dot_plot_list(string, Newname, prAB, mfAB, comment);
	free(comment);
	free(Newname);
	
	/*AA dot_plot*/
	comment=(char *)space(80*sizeof(char));
	sprintf(comment,"\n%%FreeEnergy= %.9f\n",FEAA);
	/*write New name*/
	Newname=(char*)space((strlen(ffname)+3)*sizeof(char));
	sprintf(Newname,"AA");
	strcat(Newname,ffname);
	/*write AA sequence*/
	Newstring=(char*)space((2*Alength+1)*sizeof(char));
	strcat(Newstring,Astring);
	strcat(Newstring,Astring);
	(void)PS_dot_plot_list(Newstring, Newname, prAA, mfAA, comment);
	free(comment);
	free(Newname);
	free(Newstring);

	/*BB dot_plot*/
	comment=(char *)space(80*sizeof(char));
	sprintf(comment,"\n%%FreeEnergy= %.9f\n",FEBB);
	/*write New name*/
	Newname=(char*)space((strlen(ffname)+3)*sizeof(char));
	sprintf(Newname,"BB");
	strcat(Newname,ffname);
	/*write BB sequence*/
	Newstring=(char*)space((2*Blength+1)*sizeof(char));
	strcat(Newstring,Bstring);
	strcat(Newstring,Bstring);
	/*reset cut_point*/
	cut_point=Blength+1;
	(void)PS_dot_plot_list(Newstring, Newname, prBB, mfBB, comment);
	free(comment);
	free(Newname);
	free(Newstring);

	/*A dot plot*/
	/*reset cut_point*/
	cut_point=-1;
	comment=(char *)space(80*sizeof(char));
	sprintf(comment,"\n%%FreeEnergy= %.9f\n",FEA);
	/*write New name*/
	Newname=(char*)space((strlen(ffname)+2)*sizeof(char));
	sprintf(Newname,"A");
	strcat(Newname,ffname);
	/*write BB sequence*/
	(void)PS_dot_plot_list(Astring, Newname, prA, mfA, comment);
	free(comment);
	free(Newname);

	/*B monomer dot plot*/
	comment=(char *)space(80*sizeof(char));
	sprintf(comment,"\n%%FreeEnergy= %.9f\n",FEB);
	/*write New name*/
	Newname=(char*)space((strlen(ffname)+2)*sizeof(char));
	sprintf(Newname,"B");
	strcat(Newname,ffname);
	/*write BB sequence*/
	(void)PS_dot_plot_list(Bstring, Newname, prB, mfB, comment);
	free(comment);
	free(Newname);
		
	free_franz(Astring, Bstring,  prAB, prAA, prBB, prA, prB,  mfAB, mfAA, mfBB, mfA, mfB);
	
	
      } /*end if(doT)*/
      if (do_backtrack) {
	printf("%s", structure);
	if (!istty) printf(" [%6.2f]\n", energy);
	else printf("\n");/*8.6.04*/
      }
      if ((istty)||(!do_backtrack)) 
	printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
      printf(" frequency of mfe structure in ensemble %g\n",
	     exp((energy-min_en)/kT));
     
	
    }/*end if(pf)*/
      
        
    if (do_backtrack) {
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_dp.ps");
      } else strcpy(ffname, "dot.ps");
      
      if (!doT) {
	if (pf) {	  (void) PS_dot_plot_list(string, ffname, prAB, mfAB, "doof");
	free(prAB);}
	free(mfAB);
      }
    }
    if (!doT) free_co_pf_arrays();
    
      
    if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    free(string);
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

PRIVATE char *costring(char *string)
{
  char *ctmp;
  int len;
  
  len = strlen(string);
  ctmp = (char *)space((len+2) * sizeof(char));
  /* first sequence */
  (void) strncpy(ctmp, string, cut_point-1);
  /* spacer */
  ctmp[cut_point-1] = '&';
  /* second sequence */
  (void) strcat(ctmp, string+cut_point-1);

  return ctmp;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAfold [-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n" 
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] "
	  "[-noconv] [-a] [-c] [-f concfile]\n");
}


PRIVATE double do_partfunc(char *string, int length, int Switch, struct plist **tpr, struct plist **mfpl) {
  /*compute mfe and partition function of dimere or  monomer*/  
  double En;
  char *Newstring;
  char *tempstruc;
  double min_en;
  double sfact=1.07;
  double kT;
  kT = (temperature+273.15)*1.98717/1000.;
  switch (Switch)
    {
    case 1: /*monomer*/
      cut_point=-1;
      tempstruc = (char *) space((unsigned)length+1);
      min_en = fold(string, tempstruc);
      pf_scale = exp(-(sfact*min_en)/kT/(length));
      *mfpl=get_mfe_plist(*mfpl);
      free_arrays();
      /*En=pf_fold(string, tempstruc);*/
      init_co_pf_fold(length);
      En=co_pf_fold(string, tempstruc);
      
      *tpr=get_plist(*tpr, length,0.00001);
      free_co_pf_arrays();
      free(tempstruc);
      break;

    case 2: /*dimer*/
      Newstring=(char *)space(sizeof(char)*(length*2+1)); 
      strcat(Newstring,string);
      strcat(Newstring,string);
      cut_point=length+1;
      tempstruc = (char *) space((unsigned)length*2+1);
      min_en = cofold(Newstring, tempstruc);
      pf_scale =exp(-(sfact*min_en)/kT/(2*length));
      *mfpl=get_mfe_plist(*mfpl);
      free_co_arrays(); 
      init_co_pf_fold(2*length);
      En=co_pf_fold(Newstring, tempstruc);
      *tpr=get_plist(*tpr, 2*length,0.00001);
      free_co_pf_arrays();
      free(Newstring);
      free(tempstruc);
      break;

    default:
      printf("Error in get_partfunc\n, computing neither mono- nor dimere!\n");
      exit (42);
      
    }
  return En;
}

PRIVATE void free_franz(char *Astring, char *Bstring, struct plist *prAB, struct plist *prAA, struct plist *prBB, struct plist *prA, struct plist *prB, struct plist *mfAB, struct plist *mfAA, struct plist *mfBB, struct plist *mfA, struct plist *mfB) {
  /*free arrays for dimer/monomer computations*/
  free(Astring);
  free(Bstring);
  free(prAB);
  free(prAA);
  free(prBB);
  free(prB);
  free(prA);
  free(mfAB);
  free(mfAA);
  free(mfBB);
  free(mfA);
  free(mfB);
  return;
}

PRIVATE struct ConcEnt *do_the_concenratinger(char *Conc_file,double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconces) {
  /*compute concentrations out of  free energies, calls get_concentrations*/
  struct ConcEnt *result;
  int i;
  char *line;
  double temp;
  int n=1;
  i=0;
  if (Conc_file!=NULL) i=read_concentrationfile(Conc_file,startconces);/*??*/
  else if (startconces[0]==0) {
    
    printf("Please enter concentrations, alternating A,B, '*'to stop\n");
    do {
      if ((line = get_line(stdin))==NULL) break;
      if ((line ==NULL) || (strcmp(line, "*") == 0)) break;
      sscanf(line,"%lf", &temp);
      startconces[i++]=temp;
      free(line);
      if (i==n*20-1) {
	n*=2;
	startconces=(double *)xrealloc(startconces,(2+n*20)*sizeof(double));
      }
    }while (1);
    
    if (!i%2) {
      printf("Warning! number of concentrations is not a multiple of 2!!\n");
    }
    startconces[i++]=0;
    startconces[i]=0;
  }
   
  result=get_concentrations(FEAB, FEAA, FEBB, FEA, FEB, startconces);
  
  free(startconces);
  return result;
}
PRIVATE int read_concentrationfile(char *fname, double *startc) {
  /*reads file of concentrations*/
  FILE *fp;
  char *line;
  int i,n;
  double tmp1, tmp2;
  n=1;
  i=0;
  if(!(fp = fopen(fname,"r"))){
  
    printf("Sorry, th file you specified could not be opend!!\n Please try again!!\n");
    return 0;
  }
  line=get_line(fp);
  while(line!=NULL) {
    sscanf(line,"%lf %lf",&tmp1,&tmp2);
    startc[i++]=tmp1;
    startc[i++]=tmp2;
    free(line);
    if (i==n*20) {
      n*=2;
      startc=(double *)xrealloc(startc,(2+20*n)*sizeof(double));
    }
    line=get_line(fp);
  }
  free(line);
  startc[i++]=0;
  startc[i--]=0;
  return --i;
}

