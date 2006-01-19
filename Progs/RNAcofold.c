/* Last changed Time-stamp: <2006-01-19 11:49:34 ivo> */
/*
		  c Ivo L Hofacker, Vienna RNA package
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
#include "part_func_co.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"

extern void  read_parameter_file(const char fname[]);

/*@unused@*/
static char rcsid[] = "$Id: RNAcofold.c,v 1.6 2006/01/19 11:30:04 ivo Exp $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
		       "....,....5....,....6....,....7....,....8";

PRIVATE char *costring(char *string);
PRIVATE char *tokenize(char *line);
PRIVATE void usage(void);
PRIVATE cofoldF do_partfunc(char *string, int length, int Switch, struct plist **tpr, struct plist **mf);
PRIVATE double *read_concentrations(FILE *fp);
PRIVATE void do_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconces);
PRIVATE struct plist *get_mfe_plist(struct plist *pl);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *string, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[53], ffname[20];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  char *Concfile;
  int   i, length, l, sym, r;
  double min_en;
  double kT, sfact=1.07;
  int   pf=0, istty;
  int noconv=0;
  int doT=0;    /*compute dimere free energies etc.*/
  int doC=0;    /*toggle to compute concentrations*/
  int doQ=0;    /*toggle to compute prob of base being paired*/
  int cofi=0;   /*toggle concentrations stdin / file*/
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
    if (istty) {
      printf("\nInput sequence(s); @ to quit\n");
      printf("Use '&' as spearator between 2 sequences that shall form a complex.\n");
      printf("%s\n", scale);
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
    if (doC) {
      FILE *fp;
      if (cofi) { /* read from file */
	fp = fopen(Concfile, "r");
	if (fp==NULL) {
	  fprintf(stderr, "could not open concentration file %s", Concfile);
	  nrerror("\n");
	}
	ConcAandB = read_concentrations(fp);
	fclose(fp);
      } else {
	printf("Please enter concentrations [mol/l]\n format: ConcA ConcB\n return to end\n");
	ConcAandB = read_concentrations(stdin);
      }
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
    } else {
      strcpy(ffname, "rna.ps");
    }
    if (length<2000)
      (void) PS_rna_plot(string, structure, ffname);
    else {
     fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");
     free_co_arrays();
    }

    /*compute partition function*/
    if (pf) {
      cofoldF AB, AA, BB;
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
      AB = co_pf_fold(string, structure);

      if (do_backtrack) {
	char *costruc;
	costruc = (char *) space(sizeof(char)*(strlen(structure)+2));
	if (cut_point<0) printf("%s", structure);
	else {
	  strncpy(costruc, structure, cut_point-1);
	  strcat(costruc, "&");
	  strcat(costruc, structure+cut_point-1);
	  printf("%s", costruc);
	}
	if (!istty) printf(" [%6.2f]\n", AB.FAB);
	else printf("\n");/*8.6.04*/
      }
      if ((istty)||(!do_backtrack))
	printf(" free energy of ensemble = %6.2f kcal/mol\n", AB.FAB);
      printf(" frequency of mfe structure in ensemble %g",
	     exp((AB.FAB-min_en)/kT));

      printf(" , delta G binding=%6.2f\n", AB.FcAB - AB.FA - AB.FB);

      prAB=(struct plist *) space(sizeof(struct plist) * (2*length));
      prAB=get_plist(prAB, length,0.00001);

      /* if (doQ) make_probsum(length,fname); */ /*compute prob of base paired*/
      free_co_arrays();
      if (doT) { /* cofold of all dimers, monomers */
	int Blength, Alength;
	char  *Astring, *Bstring;
	char *Newstring;
	char Newname[30];
	char comment[80];
	if (cut_point<0) {
	  printf("Sorry, i cannot do that with only one molecule, please give me two or leave it\n");
	  free(mfAB);
	  free(prAB);
	  continue;
	}
	if (dangles==1) dangles=2;
	Alength=cut_point-1;        /*length of first molecule*/
	Blength=length-cut_point+1; /*length of 2nd molecule*/

	Astring=(char *)space(sizeof(char)*(Alength+1));/*Sequence of first molecule*/
	Bstring=(char *)space(sizeof(char)*(Blength+1));/*Sequence of second molecule*/
	strncat(Astring,string,Alength);
	strncat(Bstring,string+Alength,Blength);

	/* compute AA dimer */
	prAA=(struct plist *) space(sizeof(struct plist) * (4*Alength));
	mfAA=(struct plist *) space(sizeof(struct plist) * (Alength+1));
	AA=do_partfunc(Astring, Alength, 2, &prAA, &mfAA);
	/* compute BB dimer */
	prBB=(struct plist *) space(sizeof(struct plist) * (4*Blength));
	mfBB=(struct plist *) space(sizeof(struct plist) * (Blength+1));
	BB=do_partfunc(Bstring, Blength, 2, &prBB, &mfBB);
	/*free_co_pf_arrays();*/

	/* compute A monomer */
	prA=(struct plist *) space(sizeof(struct plist) * (2*Alength));
	mfA=(struct plist *) space(sizeof(struct plist) * (Alength+1));
	do_partfunc(Astring, Alength, 1, &prA, &mfA);

	/* compute B monomer */
	prB=(struct plist *) space(sizeof(struct plist) * (2*Blength));
	mfB=(struct plist *) space(sizeof(struct plist) * (Blength+1));
	do_partfunc(Bstring, Blength, 1, &prB, &mfB);

	compute_probabilities(AB.F0AB, AB.FA, AB.FB, prAB, prA, prB, Alength);
	compute_probabilities(AA.F0AB, AA.FA, AA.FA, prAA, prA, prA, Alength);
	compute_probabilities(BB.F0AB, BB.FA, BB.FA, prBB, prA, prB, Blength);
	printf("Free Energies:\nAB\t\tAA\t\tBB\t\tA\t\tB\n%.6f\t%6f\t%6f\t%6f\t%6f\n",
	       AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB);

	if (doC) {
	  do_concentrations(AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB, ConcAandB);
	  free(ConcAandB);/*freeen*/
	}

	if (fname[0]!='\0') {
	  strcpy(ffname, fname);
	  strcat(ffname, "_dp5.ps");
	} else strcpy(ffname, "dot5.ps");
	/*output of the 5 dot plots*/

	/*AB dot_plot*/
	/*write Free Energy into comment*/
	sprintf(comment,"\n%%Heterodimer AB FreeEnergy= %.9f\n", AB.FcAB);
	/*reset cut_point*/
	cut_point=Alength+1;
	/*write New name*/
	strcpy(Newname,"AB");
	strcat(Newname,ffname);
	(void)PS_dot_plot_list(string, Newname, prAB, mfAB, comment);

	/*AA dot_plot*/
	sprintf(comment,"\n%%Homodimer AA FreeEnergy= %.9f\n",AA.FcAB);
	/*write New name*/
	strcpy(Newname,"AA");
	strcat(Newname,ffname);
	/*write AA sequence*/
	Newstring=(char*)space((2*Alength+1)*sizeof(char));
	strcpy(Newstring,Astring);
	strcat(Newstring,Astring);
	(void)PS_dot_plot_list(Newstring, Newname, prAA, mfAA, comment);
	free(Newstring);

	/*BB dot_plot*/
	sprintf(comment,"\n%%Homodimer BB FreeEnergy= %.9f\n",BB.FcAB);
	/*write New name*/
	strcpy(Newname,"BB");
	strcat(Newname,ffname);
	/*write BB sequence*/
	Newstring=(char*)space((2*Blength+1)*sizeof(char));
	strcpy(Newstring,Bstring);
	strcat(Newstring,Bstring);
	/*reset cut_point*/
	cut_point=Blength+1;
	(void)PS_dot_plot_list(Newstring, Newname, prBB, mfBB, comment);
	free(Newstring);

	/*A dot plot*/
	/*reset cut_point*/
	cut_point=-1;
	sprintf(comment,"\n%%Monomer A FreeEnergy= %.9f\n",AB.FA);
	/*write New name*/
	strcpy(Newname,"A");
	strcat(Newname,ffname);
	/*write BB sequence*/
	(void)PS_dot_plot_list(Astring, Newname, prA, mfA, comment);

	/*B monomer dot plot*/
	sprintf(comment,"\n%%Monomer B FreeEnergy= %.9f\n",AB.FB);
	/*write New name*/
	strcpy(Newname,"B");
	strcat(Newname,ffname);
	/*write BB sequence*/
	(void)PS_dot_plot_list(Bstring, Newname, prB, mfB, comment);
	free(Astring); free(Bstring); free(prAB); free(prAA); free(prBB); free(prA); free(prB);
	free(mfAB); free(mfAA); free(mfBB); free(mfA); free(mfB);

      } /*end if(doT)*/

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

PRIVATE cofoldF do_partfunc(char *string, int length, int Switch, struct plist **tpr, struct plist **mfpl) {
  /*compute mfe and partition function of dimere or  monomer*/
  char *Newstring;
  char *tempstruc;
  double min_en;
  double sfact=1.07;
  double kT;
  cofoldF X;
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
      X=co_pf_fold(string, tempstruc);

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
      X=co_pf_fold(Newstring, tempstruc);
      *tpr=get_plist(*tpr, 2*length,0.00001);
      free_co_pf_arrays();
      free(Newstring);
      free(tempstruc);
      break;

    default:
      printf("Error in get_partfunc\n, computing neither mono- nor dimere!\n");
      exit (42);

    }
  return X;
}


PRIVATE void do_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc) {
  /* compute and print concentrations out of free energies, calls get_concentrations */
  struct ConcEnt *result;
  int i, n;

  result=get_concentrations(FEAB, FEAA, FEBB, FEA, FEB, startconc);

  printf("Initial concentrations\t\trelative Equilibrium concentrations\n");
  printf("A\t\t B\t\t AB\t\t AA\t\t BB\t\t A\t\t B\n");
  for (n=0; (startconc[2*n]>0) || (startconc[2*n+1]>0); n++); /* count */
  for (i=0; i<n ;i++) {
    double tot = result[i].A0+result[i].B0;
    printf("%-10g\t%-10g\t%.5f \t%.5f \t%.5f \t%.5f \t%.5f\n", result[i].A0, result[i].B0,
	   result[i].ABc/tot, result[i].AAc/tot, result[i].BBc/tot ,result[i].Ac/tot, result[i].Bc/tot);
  }
  free(result);
}


PRIVATE double *read_concentrations(FILE *fp) {
  /* reads concentrations, returns list of double, -1. marks end */
  char *line;
  double *startc;
  int i=0, n=2;

  startc = (double *) space((2*n+1)*sizeof(double));

  while ((line=get_line(fp))!=NULL) {
    int c;
    if (i+4>=2*n) {
      n*=2;
      startc=(double *)xrealloc(startc,(2*n+1)*sizeof(double));
    }
    c = sscanf(line,"%lf %lf", &startc[i], &startc[i+1]);
    free(line);
    if (c<2) break;
    i+=2;
  }
  startc[i]=startc[i+1]=0;
  return startc;
}

static struct plist *get_mfe_plist(struct plist *pl) {
  /*get list of mfe structure out of base_pairs array*/
  int l;
  for(l=1; l<=base_pair[0].i; l++)
    {
      pl[l-1].i=base_pair[l].i;
      pl[l-1].j=base_pair[l].j;
      pl[l-1].p=0.95;
    }
  pl[l-1].i=0;
  pl[l-1].j=0;
  pl[l-1].p=0.;
  pl=(struct plist *)xrealloc(pl,l*sizeof(struct plist));
  return pl;
}


PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAcofold [-a] [-c] [-f concfile]\n"
	  "[-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n"
	  "[-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] "
	  "[-noconv]\n");
}
