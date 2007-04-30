/*
  Last changed Time-stamp: <2007-04-30 15:49:01 ulim>
  $Id: RNAup.c,v 1.1 2007/04/30 14:57:19 ivo Exp $
  
  Ineractive Access to cofolding routines
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
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
#include "pair_mat.h"
#include "part_func.h"
#include "part_func_up.h"


extern void  read_parameter_file(const char fname[]);


/*@unused@*/
static char rcsid[] = "$Id: RNAup.c,v 1.1 2007/04/30 14:57:19 ivo Exp $";

#define PRIVATE static

PRIVATE void tokenize(char *line, char **seq1, char **seq2);

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);


/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  char *string1=NULL, *string2=NULL, *temp, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[53], my_contrib[10], *up_out;
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, length1,length2,length, l, sym, r;
  double energy, min_en;
  double kT, sfact=1.07;
  int   pf, istty;
  int noconv=0;
  double Zu, Zup;
  /* variables for output */
  pu_contrib *unstr_out, *unstr_short;
  FLT_OR_DBL **inter_out;
  char *title;
  /* commandline parameters */
  int w;       /* length of region of interaction */
  int incr3;   /* add x unpaired bases after 3'end of short RNA*/
  int incr5;   /* add x unpaired bases after 5'end of short RNA*/
  int  unstr;  /* length of unpaired region for output*/
  int  upmode; /* output mode for pf_unpaired and pf_up()*/
  upmode = 0;
  unstr = 4; 
  incr3=0;
  incr5=0;
  w=25;
  do_backtrack = 1;
  pf=1; /* partition function has to be calculated */
  length1=length2=0;
  up_out=NULL;
  title=NULL;
  unstr_out=NULL;
  inter_out=NULL;
  my_contrib[0] = 'S';
  my_contrib[1] = '\0';
  
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') 
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'w':
	  /* -w maximal length of unstructured region */  
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &w);
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
	case 'o': upmode=1;
	  /* output mode 0: non, 1:only pr_unpaired, 2: pr_unpaired + pr_up */
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &upmode);
	    if (r!=1) usage();
	  }
	  break;
	case 'u':
	  /* -u length of unstructured region in pr_unpaired output
	     makes only sense in combination with -o1 or -o2 */  
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &unstr);
	  if (!r) usage();
	  break;
	  /* incr5 and incr3 are only for the longer (target) sequence */
	  /* increments w (length of the unpaired region) to incr5+w+incr3*/
	  /* the longer sequence is given in 5'(= position 1) to */
	  /* 3' (=position n) direction */
	  /* incr5 adds incr5 residues to the 5' end of w */
	case '5':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &incr5);
	  if (!r) usage();
	  break; 
	  /* incr3 adds incr3 residues to the 3' end of w */
	case '3':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &incr3);
	  if (!r) usage();
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	case 'x':  
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%s", my_contrib);
	  if (!r) usage();
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
    cut_point=-1;
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
      free(line);
      if ((line = get_line(stdin))==NULL) break;
    } 
    if ((line == NULL) || (strcmp(line, "@") == 0)) break;

    tokenize(line,&string1,&string2);
    
    if(upmode != 0){
      if(cut_point == -1 && upmode == 2) {
	  nrerror("only one sequence - can not cofold one sequence!");
      }
    } else {
      if(cut_point == -1){
	upmode=1;
      } else {
	upmode=2;
      }
    }
    
    if(string1 != NULL)
      length1 = (int) strlen(string1);
    if(string2 != NULL) 
      length2 = (int) strlen(string2);
    else
      length2=0;    

    /* write longer seq in string1 and and shorter one in string2 */ 
    if(length1 < length2)
      {
	length=length1; length1=length2; length2=length;
	
	temp=(char *) space(strlen(string1)+1);
	(void) sscanf(string1,"%s",temp);
	string1 = (char *) xrealloc (string1,sizeof(char)*length1+1);
	(void) sscanf(string2,"%s",string1);
	string2 = (char *) xrealloc(string2,sizeof(char)*length2+1);
	(void) sscanf(temp,"%s",string2);
	free(temp);
      }
   
    structure = (char *) space((unsigned) length1+1);
    if (fold_constrained) {
      cstruc = get_line(stdin);
      if (cstruc!=NULL) 
	strncpy(structure, cstruc, length1);
      else
	fprintf(stderr, "constraints missing\n");
    }
    for (l = 0; l < length1; l++) {
      string1[l] = toupper(string1[l]);
      if (!noconv && string1[l] == 'T') string1[l] = 'U';
    }
    for (l = 0; l < length2; l++) {
      string2[l] = toupper(string2[l]);
      if (!noconv && string2[l] == 'T') string2[l] = 'U';
    }

    if (istty)
      printf("length1 = %d\n", length1);
    
    /* initialize_fold(length); */
    update_fold_params();
    printf("\n%s", string1);
    min_en = fold(string1, structure);
    
    if (istty)
      {
	printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
      }
    else
      printf(" (%6.2f)\n", min_en);
    
    (void) fflush(stdout);
    
    /* parse cml parameters for the filename*/
    if(upmode > 0) {
      char wuadd[10];
      up_out = (char*) space(sizeof(char)*53);
      /* create the name of the output file */
      if(fname[0]!='\0' && up_out[0] =='\0' ){
	if(strlen(fname)< 30){
	  strcpy(up_out, fname);
	} else {  
	  strncpy(up_out, fname,30);
	}
      }
      else if(fname[0]=='\0' && up_out[0] == '\0'){
	char defaultn[10] = "RNA";
	sprintf(up_out,"%s",defaultn);
      }
	
      sprintf(wuadd,"%d",w);
      strcat(up_out, "_w");
      strcat(up_out, wuadd);
      strcat(up_out, "u");
      sprintf(wuadd,"%d",unstr);
      strcat(up_out, wuadd);
      strcat(up_out, "_up.out");
      printf("RNAup output in file: %s\n",up_out);
	    
      /* create the title for the output file */      
      if (title == NULL) {
	char wuadd[10];
	title = (char*) space(sizeof(char)*60);
	if(fname[0]!='\0'){
	  if(strlen(fname)< 30){
	    strcpy(title, fname);
	  } else {  
	    strncpy(title, fname,30);
	  }
	}
	else if (fname[0]=='\0'){
	  char defaultn[10]= "RNAup";
	  sprintf(title,"%s",defaultn);
	}
	sprintf(wuadd,"%d",unstr);
	strcat(title," u=");
	strcat(title, wuadd);
	sprintf(wuadd,"%d",w);
	strcat(title," w=");
	strcat(title, wuadd);
	sprintf(wuadd,"%d",length1);
	strcat(title," n=");
	strcat(title, wuadd);
      }
    } else {
      nrerror("no output format given: use [-o[1|2]] to select output format");
    }
    
    
    if (pf) {
      
      if (dangles==1) {
	dangles=2;   /* recompute with dangles as in pf_fold() */
	min_en = energy_of_struct(string1, structure);
	dangles=1;
      }
	 
      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      
      if(upmode != 0){
	int wplus;
	wplus=w+incr3+incr5;
	/* calculate prob. unstructured for the shorter seq */
	if(upmode == 3) {
	  min_en = fold(string2, structure);
	  pf_scale = exp(-(sfact*min_en)/kT/length2);
	  if (length2>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
	  init_pf_fold(length2);
	  if (cstruc!=NULL)
	    strncpy(structure, cstruc, length2+1);
	  energy = pf_fold(string2, structure);
	  if(wplus > length2){ wplus = length2;} /* for the shorter seq */
	  unstr_short = pf_unstru(string2, structure, wplus);
	  free_pf_unstru();
	  free_pf_arrays(); /* for arrays for pf_fold(...) */
	}

	/* calculate prob. unstructured for the longer seq */
	wplus=w+incr3+incr5; 
	min_en = fold(string1, structure);
	pf_scale = exp(-(sfact*min_en)/kT/length1);
	if (length1>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
	init_pf_fold(length1);
	if (cstruc!=NULL)
	  strncpy(structure, cstruc, length1+1);
	energy = pf_fold(string1, structure);
	unstr_out = pf_unstru(string1, structure, wplus);
	free_pf_unstru();
	free_pf_arrays(); /* for arrays for pf_fold(...) */
	/* calculate the interaction between the two sequences */
	if(upmode > 1 && cut_point > -1){
	  inter_out = pf_interact(string1,string2,unstr_out,w, incr3, incr5);
	  if(Up_plot(unstr_out,inter_out,length1,up_out,unstr,my_contrib)==0){
	    nrerror("Up_plot: no output values assigned");
	  }
	} else if(cut_point == -1 && upmode > 1) { /* no second seq given */
	  nrerror("only one sequence given - cannot cofold one sequence!");
	} else { /* plot only the results for prob unstructured */
	  if(Up_plot(unstr_out,NULL,length1,up_out,unstr,my_contrib)==0){
	    nrerror("Up_plot: no output values assigned");
	  }
	}	
      } else {
	nrerror("no output format given: use [-o[1|2]] to select output format");
      }
       
      if (do_backtrack) {
	printf("%s", structure);
	if (!istty) printf(" [%6.2f]\n", energy);
	else printf("\n");
      }
      if ((istty)||(!do_backtrack)) 
	printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
      energy = pf_fold(string1, structure);
      printf(" frequency of mfe structure in ensemble %g; "
	     "ensemble diversity %-6.2f\n", exp((energy-min_en)/kT),
	     mean_bp_dist(length1));
      free_pf_arrays();
    }
    if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    if (string1!=NULL) free(string1);
    if (string2!=NULL) free(string2);
    free(structure);
    if(up_out != NULL) free(up_out);
    up_out=NULL;
    if(title != NULL) free(title);
    title=NULL;
    if(upmode == 1) free_pf_two(unstr_out,NULL);
    if(upmode > 1) free_pf_two(unstr_out,inter_out);
    if(upmode == 3)free_pf_two(unstr_short,NULL);
    free_arrays(); /* for arrays for fold(...) */
    
  } while (1);
  return 0;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAup [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n"  
	  "      [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]\n"
          "      [-w max.lenght of unpaired region] [-5 length] [-3 length]\n"
          "      [-o[1|2]] [-u length of unpaired region in output]\n"
	  "      [-x string determining p_unpaired output][-noconv]\n");
}

/* call:  tokenize(line,&seq1,&seq2); the sequence string is spli at the "&"
   and the first seq is written in seq1, the second into seq2 */
void tokenize(char *line, char **seq1, char **seq2) {
  char *pos;
  int cut = -1;
  int i;
  pos = strchr(line, '&');
  if (pos) {
    cut = (int) (pos-line)+1;
    (*seq1) = (char *) space((cut+1)*sizeof(char));
    (*seq2)= (char *) space(((strlen(line)-cut)+2)*sizeof(char));
  
    if (cut >= strlen(line)) cut = -1;
    if (strchr(pos+1, '&')) nrerror("more than one cut-point in input");
    strncpy((*seq1),line,cut-1);
    (*seq1)[cut-1]='\0';
    for (i=0;i<=strlen(line)-cut;i++) {
      (*seq2)[i]=line[i+cut];   
    } /* splice out the & */
    (*seq2)[i]='\0';
      
  } else {
    (*seq1) = (char *) space((strlen(line)+1)*sizeof(char));
    (*seq2)=NULL;
    strcpy((*seq1),line);
  }
    
  if (cut > -1) {
    if (cut_point==-1) cut_point = cut;
    else if (cut_point != cut) {
      fprintf(stderr,"cut_point = %d cut = %d\n", cut_point, cut);
      nrerror("Sequence and Structure have different cut points.");
    }
  }
  free(line);
  line=NULL;
  return;
}
