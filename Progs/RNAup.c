/*
  Last changed Time-stamp: <2007-12-20 17:52:27 ulim>
  $Id: RNAup.c,v 1.3 2007/12/22 09:03:01 ivo Exp $
  
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
#include <float.h> 
#include "fold.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
#include "part_func.h"
#include "part_func_up.h"
#include "duplex.h"
#include "energy_const.h"

extern void  read_parameter_file(const char fname[]);


/*@unused@*/
static char rcsid[] = "$Id: RNAup.c,v 1.3 2007/12/22 09:03:01 ivo Exp $";

#define PRIVATE static
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define EQUAL(A,B) (fabs((A)-(B)) < 1000*DBL_EPSILON)
PRIVATE void tokenize(char *line, char **seq1, char **seq2);
PRIVATE char *tokenize_one(char *line);
PRIVATE int comp_nums(const int *num1, const int *num2);
PRIVATE int get_u_values(char unstrs[], int **u_vals, int l1);
PRIVATE void seperate_bp(char **inter,int len1,char **intra_l,char **intra_s);
PRIVATE void print_interaction(interact *Int, char *s1, char *s2, pu_contrib *p_c, pu_contrib *p_c2, int w, int incr3, int incr5);
PRIVATE void print_unstru(pu_contrib *p_c, int w);

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);
/* defaults for -u and -w */
PRIVATE int default_u; /* -u options for plotting: plot pr_unpaired for 4 nucleotides */
PRIVATE int default_w; /* -w option for interaction: maximal region of interaction is 25 nucleotides */

/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  char *string1=NULL, *string2=NULL, *dummy=NULL, *temp=NULL, *line=NULL;
  char *structure=NULL, *cstruc=NULL, *cstruc_l=NULL, *cstruc_s=NULL;
  char fname[53], ffname[53], temp_name[53], first_name[53], my_contrib[10];
  char up_out[60], unstrs[100];
  char *ParamFile=NULL;
  char *ns_bases=NULL, *c;
  int  i, length1,length2,length, l, sym, r, *u_vals;
  double energy, min_en;
  double kT, sfact=1.07;
  int   pf, istty;
  int noconv=0;
  /* variables for output */
  pu_contrib *unstr_out, *unstr_short;
  interact *inter_out;
  char *title;
  /* commandline parameters */
  int w;       /* length of region of interaction */
  int incr3;   /* add x unpaired bases after 3'end of short RNA*/
  int incr5;   /* add x unpaired bases after 5'end of short RNA*/
  int unstr;   /* length of unpaired region for output*/
  int upmode ;  /* output mode for pf_unpaired and pf_up()*/
  int task;    /* input mode for calculation of interaction */
  /* default settings for RNAup */
  task=0;
  upmode = 1; /* default is one sequence, option -X[p|f] has to be set
		 for the calculation of an interaction, if no "&" is in
		 the sequence string  */
  unstrs[0]='\0';
  default_u = 4;
  unstr=default_u;
  default_w = 25;
  w=default_w;
  u_vals=NULL;
  incr3=0;
  incr5=0;
  do_backtrack = 1;
  pf=1; /* partition function has to be calculated */
  length1=length2=0;
  title=NULL;
  unstr_out=NULL;
  inter_out=NULL;
  my_contrib[0] = 'S';
  my_contrib[1] = '\0';
  first_name[0] = '\0';
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
	case 'b':
	  upmode=3;
	  break;
	case 'X':
	  /* interaction mode invoked */
	  if(upmode == 1) upmode=2;
	  switch (argv[i][2]) { /* now determine which sequences interact */
	  case 'p': task=1; break; /* pairwise interaction */
	  case 'f': task=2; break; /* first one interacts with all others */
	  }
	  break;
	case 'u':
	  /* -u length of unstructured region in pr_unpaired output */  
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%s", unstrs);
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
	case 'c':  
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
    printf(". : no constraint at all\n");
    printf("x : base must not pair\n");
    printf("matching brackets ( ): base i pairs base j\n");
    printf("constraints for unimolecular folding only:\n"); 
    printf("< : base i is paired with a base j<i\n");
    printf("> : base i is paired with a base j>i\n");    
    printf("constraints for cofolding only:\n");
    printf("| : paired with another base intermolecularly\n");        
  } 
	
  do {	/* main loop: continue until end of file */
    cut_point=-1;
    if (istty) {
      if(upmode == 1) {
	printf("\nInput string (upper or lower case); @ to quit\n");
	printf("%s%s\n", scale1, scale2);
      }
      else if (upmode > 1) {
	if(task == 1 || (task == 0 && upmode == 3)) {
	  printf("\nUse either '&' to connect the 2 sequences or give each sequence on an extra line.\n"); 
	  printf("%s%s\n", scale1, scale2);
	}
	else if(task == 2) {
	  printf("\nGive each sequence on an extra line. The first seq. is stored, every other seq. is compared to the first one.\n"); 
	  printf("%s%s\n", scale1, scale2);
	}
      }
    }
    fname[0]='\0';
    ffname[0]='\0';
    /* read the first sequence */
    if ((line = get_line(stdin))==NULL) break;

    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
      if (*line=='>')
	(void) sscanf(line, ">%51s", fname);
      free(line);
      line=NULL;
      if ((line = get_line(stdin))==NULL) break;
    } 
    if ((line == NULL) || (strcmp(line, "@") == 0)) break;

    if(first_name[0] == '\0' && fname[0] !='\0' && task == 2) {
	strncpy(first_name,fname,30);
      }
    /* if upmode == 2: check if the sequences are seperated via "&" (cut_point > -1) or given on extra lines */
    if(task < 3) {
      tokenize(line,&string1,&string2);
      if(task == 2 && cut_point != -1) task = 3;
      /* two sequences with & are given: calcuklate interaction */
      if(task == 0 && cut_point != -1) {
	task = 1;
	if(upmode == 1) upmode = 2;
      }
    }
    else if(task == 3) { /* option -Xf*/
      strncpy(ffname,fname,30);
      strncpy(fname,first_name,30);  /* first_name: name of first seq */
      if(temp != NULL) { /*strings have been switched - write temp to string1*/
	string1 = (char *) xrealloc (string1,sizeof(char)*strlen(temp)+1);
	(void) sscanf(temp,"%s",string1);
	free(temp);temp=NULL;
	
      }	
      tokenize(line,&string2,&dummy); /*compare every seq to first one given */
      free(dummy);dummy=NULL;
      if(cut_point != -1) {
	nrerror(
	   "After the first sequence pair: Input a single sequence (no &)!\n"
	   "Each input seq. is compared to the first seq. given.\n");
      }
    }
    /* interaction mode -> get the second seq. if seq are on seperate lines*/
    if(upmode > 1){ /* interaction mode */
      if(cut_point == -1 && task < 3) { /* seqs are given on seperate lines */
	/* read the second sequence */
	if(task == 2) task = 3;
	if ((line = get_line(stdin))==NULL) {
	  nrerror("only one sequence - can not cofold one sequence!");
	}
	/* skip comment lines and get filenames */
	while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	  if (*line=='>')
	    (void) sscanf(line, ">%51s", ffname); /* name of the 2nd seq */
	  free(line);
	  line=NULL;
	  if ((line = get_line(stdin))==NULL) break;
	} 
	if ((line ==NULL) || (strcmp(line, "@") == 0)) break;
	free(string2); /* string2 has been allocated in tokenize() */
    
	string2 = (char *) space(strlen(line)+1);
	(void) sscanf(line,"%s",string2); free(line);line=NULL;
      }
    } else { /* default mode pr_unpaired for ONE seq */
      /* if a second sequence is give, cofold the sequences*/
      if(cut_point != -1){
	upmode = 2;	
      }
    }

    /* parse cml parameters for output filename*/    
    /* create the name of the output file */
    if(fname[0]!='\0') {
      strncpy(up_out,fname,30);
      if(upmode >1 && ffname[0] != '\0') {
	strncpy(up_out,fname,15);
	strcat(up_out, "_");
	strncat(up_out,ffname,15);
      }	
    } else {
      char defaultn[10] = "RNA";
      sprintf(up_out,"%s",defaultn);
    }
    if(upmode >1) {
      sprintf(temp_name,"%d",w);
      strcat(up_out, "_w");
      strcat(up_out, temp_name);
    }

    if(string1 != NULL){length1 = (int) strlen(string1);}
    else {nrerror("sequence is NULL, check your input.");}
	    
    /* do this only when -X[p|f] is used or if two sequences seperated by & are given */
    if(upmode > 1) {
      if(task == 3) strncpy(temp_name,fname,30);
      
      if(string2 != NULL) {length2 = (int) strlen(string2);}
      else{nrerror("one of the sequences is NULL, check your input.");}

      /* write longer seq in string1 and and shorter one in string2 */ 
      if(length1 < length2) {
	strncpy(temp_name,fname,30);
	strncpy(fname,ffname,30);
	strncpy(ffname,temp_name,30);
	  
	length=length1; length1=length2; length2=length;
	
	temp=(char *) space(sizeof(char)*strlen(string1)+1);
	(void) sscanf(string1,"%s",temp);
	string1 = (char *) xrealloc (string1,sizeof(char)*length1+1);
	(void) sscanf(string2,"%s",string1);
	string2 = (char *) xrealloc(string2,sizeof(char)*length2+1);
	(void) sscanf(temp,"%s",string2);
	if(task == 1) {
	  free(temp);
	  temp = NULL;
	}
      }
    }
    
    /* get values for -u */
    if( ! get_u_values(unstrs,&u_vals,length1)) {
      nrerror("option -u: length value exceeds sequence length\n");
    }
      
    
    for (l = 0; l < length1; l++) {
      string1[l] = toupper(string1[l]);
      if (!noconv && string1[l] == 'T') string1[l] = 'U';
    }
    for (l = 0; l < length2; l++) {
      string2[l] = toupper(string2[l]);
      if (!noconv && string2[l] == 'T') string2[l] = 'U';
    }
    
    if (fold_constrained) {
      char *temp_cstruc=NULL;
      int old_cut;
      temp_cstruc = get_line(stdin);
      old_cut = cut_point;
      cut_point=-1;
      /* get contrained string without & */
      cstruc = tokenize_one(temp_cstruc);
      /* free(temp_cstruc); */
      /* only one seq, cstruc should not have an & */
      if(upmode == 1 && cut_point == -1) {
	if(strlen(cstruc) == length1) {
	  cstruc_l=(char*)space(sizeof(char)*(length1+1));
	  strncpy(cstruc_l,cstruc,length1);
	}else{
	  fprintf(stderr, "%s\n%s\n",string1,cstruc);
	  nrerror("RNAup -C: constrain string and structure have unequal length");
	}
      }	else if(upmode == 1 && cut_point != -1) {
	fprintf(stderr, "%s\n%s\n",string1,cstruc);
	nrerror("RNAup -C: only one sequence but constrain structure for cofolding");
      }
      /* constrain string is for both seqs */
      else if(upmode > 1 && cut_point != -1) {
	if(old_cut != cut_point) {
	  nrerror("RNAup -C: different cut points in sequence und constrain string");
	}
	seperate_bp(&cstruc,length1,&cstruc_l,&cstruc_s);
	if(strlen(cstruc) != (length1+length2)) {
	  fprintf(stderr, "%s&%s\n%s\n",string1,string2,cstruc);
	  nrerror("RNAup -C: constrain string and structure have unequal length");
	}
	if(strlen(cstruc_l) != (length1)) {
	  fprintf(stderr, "%s\n%s\n",string1,cstruc_l);
	  nrerror("RNAup -C: constrain string and structure have unequal length");
	}
	if(strlen(cstruc_s) != (length2)) {
	  fprintf(stderr, "%s\n%s\n",string2,cstruc_s);
	  nrerror("RNAup -C: constrain string and structure have unequal length");
	} 
      } else {
	fprintf(stderr, "%s&%s\n%s\n",string1,string2,cstruc);
	nrerror("RNAup -C: no cutpoint in constrain string");
      }      
    }
    
    structure = (char *) space(sizeof(char)*(length1+1));
    update_fold_params();
    if (cstruc_s != NULL)
      strncpy(structure, cstruc_s, length2+1);
    min_en = fold(string1, structure);    
    (void) fflush(stdout);

    if (pf) {      
      if (dangles==1) {
	dangles=2;   /* recompute with dangles as in pf_fold() */
	min_en = energy_of_struct(string1, structure);
	dangles=1;
      }	 
      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      if(upmode != 0){
	int wplus,w_sh;
	if(upmode == 3) { /* calculate prob. unstruct. for shorter seq */  
	  w_sh = w;
	  /* len of unstructured region has to be <= len shorter seq. */
	  if(w > length2) w_sh = length2;
	  if (cstruc_s != NULL)
	    strncpy(structure, cstruc_s, length2+1);
	  min_en = fold(string2, structure);	  
	  pf_scale = exp(-(sfact*min_en)/kT/length2);
	  if (length2>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
	  init_pf_fold(length2);
	  if (cstruc_s != NULL)
	    strncpy(structure, cstruc_s, length2+1);
	  energy = pf_fold(string2, structure);
	  unstr_short = pf_unstru(string2, w_sh);
	  free_pf_arrays(); /* for arrays for pf_fold(...) */
	}

	/* calculate prob. unstructured for longer seq */
	wplus=w+incr3+incr5;
	/* calculate prob. unpaired for the maximal length of -u */
	if(u_vals[u_vals[0]] > wplus) wplus=u_vals[u_vals[0]];
	/* length of the unstructured region has to be <= len longer seq. */
	if(wplus > length1) wplus=length1;
	if (cstruc_l !=NULL)
	  strncpy(structure, cstruc_l, length1+1);
	min_en = fold(string1, structure);
	pf_scale = exp(-(sfact*min_en)/kT/length1);
	if (length1>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
	init_pf_fold(length1);
	if (cstruc_l !=NULL)
	  strncpy(structure, cstruc_l, length1+1);
	energy = pf_fold(string1, structure);
	if(upmode > 1) {
	  unstr_out = pf_unstru(string1, wplus);
	} else {
	  unstr_out = pf_unstru(string1, u_vals[u_vals[0]]);
	}
	free_pf_arrays(); /* for arrays for pf_fold(...) */

	if(upmode > 1){/* calculate interaction between two sequences */
	  int count;
	  if(upmode == 2) {
	    inter_out = pf_interact(string1,string2,unstr_out,NULL,w,cstruc,incr3,incr5);
	    print_interaction(inter_out,string1,string2,unstr_out,NULL,w,incr3,incr5);
	  } else if(upmode == 3){
	    inter_out = pf_interact(string1,string2,unstr_out,unstr_short,w,cstruc,incr3,incr5);
	    print_interaction(inter_out,string1,string2,unstr_out,unstr_short,w,incr3,incr5);
	  }
	  if(u_vals[0] == 1) { printf("RNAup output in file: ");}
	  else { printf("RNAup output in files: ");}
	  /* plot for all -u values */
	  for(count = 1; count <= u_vals[0]; count++) {
	    char name[60];
	    unstr = u_vals[count];

	    strcpy(name,up_out);
	    strcat(name, "_u");
	    sprintf(temp_name,"%d",unstr);
	    strcat(name, temp_name);
	    strcat(name, "_up.out");
	    printf("%s",name);
	    if(count < u_vals[0]) printf(", ");
	    if(upmode == 2) {
	      Up_plot(unstr_out,NULL,inter_out,length1,name,unstr,my_contrib);
	    }
	    else if (upmode == 3){
	      Up_plot(unstr_out,unstr_short,inter_out,length1,name,unstr,my_contrib);	      
	    }
	  }
	  printf("\n");
	} else { /* one sequence:  plot only results for prob unstructured */
	  int count;
	  char collect_out[1000];
	  collect_out[0]='\0';
	  for(count = 1; count <= u_vals[0]; count++) {
	    unstr = u_vals[count];
	    print_unstru(unstr_out,unstr);
	  }
	  if(u_vals[0] == 1) { printf("RNAup output in file: ");}
	  else { printf("RNAup output in files: ");}
	  for(count = 1; count <= u_vals[0]; count++) {
	    char name[60];
	    unstr = u_vals[count];

	    strcpy(name,up_out);
	    strcat(name, "_u");
	    sprintf(temp_name,"%d",unstr);
	    strcat(name, temp_name);
	    strcat(name, "_up.outTESTV");
	    printf("%s",name);
	    if(count < u_vals[0]) printf(", ");
	    Up_plot(unstr_out,NULL,NULL,length1,name,unstr,my_contrib);
	  }
	  printf("\n");
	}	
      } else {
	nrerror("no output format given\n");
      }
    }
    
    free(structure);
    if(title != NULL) free(title);
    title=NULL;
    if(u_vals != NULL) free(u_vals);
    u_vals=NULL;
    if(upmode == 1) free_pu_contrib(unstr_out);
    if(upmode > 1) {
      free_pu_contrib(unstr_out);
      free_interact(inter_out);
    }
    if(upmode == 3)free_pu_contrib(unstr_short);
    free_arrays(); /* for arrays for fold(...) */   
    if (cstruc!=NULL) free(cstruc);
    cstruc=NULL;
    if (cstruc_l!=NULL) free(cstruc_l);
    cstruc_l=NULL;
    if (cstruc_s!=NULL) free(cstruc_s);
    cstruc_s=NULL;
    (void) fflush(stdout);
    if (string1!=NULL && task != 3) {
      free(string1);
      string1 = NULL;
    }
    if (string2!=NULL) free(string2);
    string2 = NULL;
    
  } while (1);
  if (string1!=NULL) free(string1);
  if (string2!=NULL) free(string2);
  if (cstruc!=NULL) free(cstruc);
  if (cstruc_l!=NULL) free(cstruc_l);
  if (cstruc_s!=NULL) free(cstruc_s);  
  return 0;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAup [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n"  
	  "      [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]\n"
          "      [-noconv] [-u len] [-c \"SHIME\"] [-X[p|f] interaction]\n"
	  "      [-w len] [-5 length] [-3 length] [-b]\n");
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
/* remove the & from a string for two sequences */
PRIVATE char *tokenize_one(char *line)
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
  line=NULL;
  return copy;
}

int comp_nums(const int *num1, const int *num2) {
  if (*num1 <  *num2) return -1;
  if (*num1 == *num2) return  0;
  if (*num1 >  *num2) return  1;
  return 0;
}

/* get the values for option -u, write them in u_vals*/
/* max. length of the unstructured region has to be <= len longer seq.!!!*/
/* u_vals[u_vals[0]] contains the largest -u value <= len longer seq. */
int get_u_values(char unstrs[], int **u_vals, int l1) {
  int min,max,tmp,uc,count,uunstr;
  char *token, *cp;

  if((*u_vals) != NULL) free((*u_vals));
  (*u_vals) = (int*) space (100*sizeof(int));
  if(unstrs[0] != '\0' && strchr(unstrs,'-')) {/*range contains symbol "-"*/
    const char delimiters[] = " -";

    if(strchr(unstrs,','))
       nrerror("option -u : enter either a range using \"-\" or a comma seperated list\n");
    cp = strdup(unstrs);                
    token = strtok(cp,delimiters);     
    min = atoi(token);
    token = strtok (NULL,delimiters);
    max = atoi(token);
    free(cp);
    if(min > max) {
      tmp = min;
      min = max;
      max = tmp;
    } else if (min == max) {
      nrerror("option -u : you enterd a range where min = max, use min < max to define a range");
    }
    if(max - min > 20) {
      fprintf(stderr, "only the first 20 length value are used\n");
    }
    
    (*u_vals)[0] = (max - min+1) <= 20 ? (max - min+1) : 20;
    uc = 0;
    max = max < min+19 ? max : min+19;
    for(tmp = min;tmp <= max; tmp++) {
      if(tmp <= l1) { (*u_vals)[++uc]=tmp;
	/* printf("%d,",tmp); */
      } else {
	fprintf(stderr, "option -u: length %d is longer than length of longer sequence. Only values <= length of longer sequence are allowed.\n",tmp);
	break;
      }
    }
    (*u_vals)[0]=uc;
    if(uc < 1) return(0);
    return(1);
    /* comma seperated list of values, symbol "," */
  } else if (unstrs[0] != '\0' && strchr(unstrs,',')) {
    const char delimiters[] = " ,";
    if(strchr(unstrs,'-'))
       nrerror("option -u : enter either a range using \"-\" or a comma seperated list\n");
    
    cp = strdup(unstrs);               
    token = strtok (cp,delimiters);
    uc = 1;
    (*u_vals)[1] = atoi(token);
    while((token=strtok(NULL,delimiters)) && uc<20 )
      (*u_vals)[++uc] = atoi(token);
    if((token=strtok(NULL,delimiters))) {
      fprintf(stderr,"the first 20 length value are used\n");
    }
    free(cp);
    (*u_vals)[0] = 0;
    uunstr = (uc) <= 20 ? uc : 20;
    qsort((*u_vals),(uunstr+1),sizeof(int),(void *)comp_nums );
    for(count = 0; count < uunstr+1; count++) {
      if((*u_vals)[count] > l1) {
	fprintf(stderr, "option -u: length %d is longer than length of longer sequence. Only values <= length of longer sequence are allowed.\n",(*u_vals)[count]);
	break;
      }      
    }
    (*u_vals)[0] = count - 1;
    if((count -1) < 1) return(0);   
    return(1);
  } else if (unstrs[0] != '\0') {
    (*u_vals)[0] = 1;
    uunstr = atoi(unstrs);
    if(uunstr > l1) return(0);
    (*u_vals)[1] = uunstr;
    return(1);
  } else { /* default value */
    (*u_vals)[0] = 1;
    if(default_u > l1) {
      (*u_vals)[1] = l1;
      fprintf(stderr, "option -u = %d exceeds length of longer sequence, %d. -u is set length of longer sequence.\n",default_u,l1);
    }
    (*u_vals)[1] = default_u;
  }
  return 1;
}
/* dived the constrains string in intermolecular constrains (inter)
   and intramolecular constrains within both sequences */
/* len1 is the length of the LONGER input seq ! */
void seperate_bp(char **inter, int len1, char **intra_l, char **intra_s) {
  int i,j,len;
  short *pt=NULL;
  char *temp_inter, *pt_inter;

  len=strlen((*inter));  
  /* printf("inter\n%s\n",(*inter)); */
  i = len+1;
  temp_inter=(char*)space(sizeof(char)*i);
  /* to make a pair_table convert <|> to (|) */
  pt_inter=(char*)space(sizeof(char)*i); 
  /* if shorter seq is first seq in constrained string, write the
     longer one as the first one */				
  temp_inter[strlen((*inter))] = '\0';
  pt_inter[strlen((*inter))] = '\0';
  if(cut_point < len1) {
    /* write the constrain for the longer seq first */
    for(j=0,i=cut_point-1;i<len;i++,j++) {
      switch ((*inter)[i]){
      case '(':
	temp_inter[j] = ')';
	pt_inter[j] = ')';
	break;
      case ')':
	temp_inter[j] = '(';
	pt_inter[j] = '(';
	break;
      default:
	temp_inter[j] = (*inter)[i];
	pt_inter[j] = '.';
      }
    }
    /* then add the constrain for the shorter seq */
    for(i=0;i< cut_point-1;i++,j++) {
      switch ((*inter)[i]){
      case '(':
	temp_inter[j] = ')';
	pt_inter[j] = ')';
	break;
      case ')':
	temp_inter[j] = '(';
	pt_inter[j] = '(';
	break;
      default:
	temp_inter[j] = (*inter)[i];
	pt_inter[j] = '.';
      }
    }
    cut_point = len1+1;
    strcpy((*inter),temp_inter);
  } else {
    for(i=0;i<strlen((*inter));i++) {
      switch ((*inter)[i]){
      case '(':
	pt_inter[i] = '(';
	break;
      case ')':
	pt_inter[i] = ')';
	break;
      default:
	pt_inter[i] = '.';
      }
    }
  }    
	
  pt = make_pair_table(pt_inter);
  /* printf("temp_inter\n%s\n%s\t%d\n",temp_inter,pt_inter,cut_point); */
  /* for(i=0;i<=pt[0];i++) { printf("%d   %d\n",i,pt[i]);} */
/*   printf("\n"); */
  /* intramolecular structure in longer (_l) and shorter (_s) seq */
  (*intra_l)=(char*)space(sizeof(char)*(len1+1));
  (*intra_s)=(char*)space(sizeof(char)*(strlen((*inter))-len1+2));
  (*intra_l)[len1] = '\0';
  (*intra_s)[strlen((*inter))-len1+1] = '\0';
  /* now seperate intermolecular from intramolecular bp */
  for(i=1;i<=pt[0];i++) {
    if(pt[i] == 0) {
      temp_inter[i-1] = (*inter)[i-1];
      if(i<cut_point) {
	(*intra_l)[i-1] = (*inter)[i-1];
	if((*inter)[i-1] == '|')
	  (*intra_l)[i-1] = '.';
      } else {
	(*intra_s)[i-cut_point] = (*inter)[i-1];
	if((*inter)[i-1] == '|')
	  (*intra_s)[i-cut_point] = '.';
      }
    } else {
      if(i<cut_point) {
	/* intermolekular bp */
	if(pt[i]>=cut_point){
	  temp_inter[i-1] = (*inter)[i-1];
	  (*intra_l)[i-1] = '.';
	  (*intra_s)[pt[i]-cut_point] = '.';
	} else { /* intramolekular bp */
	  (*intra_l)[i-1] = (*inter)[i-1];
	  temp_inter[i-1] = '.';
	}
      } else { /* i>=cut_point */
	/* intermolekular bp */
	if(pt[i] < cut_point){
	  temp_inter[i-1] = (*inter)[i-1];
	  /* (*intra_s)[i-1] = '.'; */
	} else { /* intramolekular bp */
	  (*intra_s)[i-cut_point] = (*inter)[i-1];
	  temp_inter[i-1] = '.';
	}
      }
    }
  }
  
  /* printf("%s -1\n%s -2\n%s -3\n%s -4\n",(*inter),temp_inter,(*intra_l),(*intra_s)); */
  strcpy((*inter),temp_inter);
  free(temp_inter);
  free(pt_inter);
  free(pt);
}

PRIVATE void print_interaction(interact *Int, char *s1, char *s2, pu_contrib *p_c, pu_contrib *p_c2, int w, int incr3, int incr5) {
  char *i_long,*i_short;
  int len, l_l, l_s, pc_size, len1, i, j, n1, n2, end5, end3;
  double **p_c_S, **p_c2_S;
  double G_min,Gi_min,Gul;
  duplexT mfe;

  
  G_min = Int->Gikjl;
  Gi_min = Int->Gikjl_wo;
  len1 = Int->length;
  len=strlen(s1)+strlen(s2);
  n1 = (int) strlen(s1); 
  n2 = (int) strlen(s2);
  p_c_S = (double **) space (sizeof(double *)*(n1+1));
  
  for (i=1; i<=n1; i++) {
    pc_size = MIN((w+incr3+incr5),n1);
    /* printf("pc_size = %d\n",pc_size); */
    p_c_S[i] = (double *) space (sizeof(double)*(pc_size+1));    
    for (j=0; j < (pc_size); j++) {
      p_c_S[i][j] = p_c->H[i][j]+p_c->I[i][j]+p_c->M[i][j]+p_c->E[i][j];
      
    }
  }
  if(p_c2 != NULL) {
    p_c2_S = (double **) space (sizeof(double *)*(n2+1));
    for (i=1; i<=n2; i++) {
      pc_size = MIN(w,n2);
      p_c2_S[i] = (double *) space (sizeof(double)*(pc_size+2));    
      for (j=0; j < (pc_size); j++) {
	p_c2_S[i][j] = p_c2->H[i][j]+p_c2->I[i][j]+p_c2->M[i][j]+p_c2->E[i][j];
      }
    }
  }
  
  /* use duplexfold() to fold the interaction site */
  l_l = (Int->i-Int->k+1);
  i_long = (char*)space(sizeof(char)*(l_l+1));
  l_s = (Int->l-Int->j+1);
  i_short = (char*)space(sizeof(char)*(l_s+1));
  
  strncpy(i_long,&s1[Int->k-1],(Int->i-Int->k+1));
  i_long[Int->i-Int->k+1] = '\0';
  strncpy(i_short,&s2[Int->j-1],(Int->l-Int->j+1));
  i_short[Int->l-Int->j+1] = '\0';  
  
  mfe = duplexfold(i_long,i_short);
  

  free(i_long);
  free(i_short);

  end5=Int->k-incr5 > 0 ? Int->k-incr5 : 1;
  pc_size = MIN((w+incr3+incr5),len1);
  end3= Int->i-Int->k+incr3 < pc_size ? Int->i-Int->k+incr3 : (pc_size-1);
  Gul = -log(p_c_S[end5][end3])*((temperature+K0)*GASCONST/1000.0);

  if(p_c2 == NULL) {
    double G_sum,diff;
    G_sum =  Gi_min + Gul;
    
    if(!EQUAL(G_min,G_sum)) {
      printf("ERROR\n");
      diff = fabs((G_min)-(G_sum));
      printf("diff %.18f\n",diff);
    }
    printf("dG = dGint + dGu_l\n");
    printf("%s %3d,%-3d : %3d,%-3d (%.2f = %.2f + %.2f)\n", mfe.structure, Int->k, Int->i, Int->j, Int->l, G_min, Gi_min, Gul);
  } else {    
     /* get values when calculating k<i and j<l */
    double Gus,G_sum,diff;
    Gus = -log(p_c2_S[Int->j][(Int->l)-(Int->j)])*((temperature+K0)*GASCONST/1000.0);
    G_sum = Gi_min+Gul+Gus;
    if(!EQUAL(G_min,G_sum)) {
      printf("ERROR\n");
      diff = fabs((G_min)-(G_sum));
      printf("diff %.18f\n",diff);
    }     
    printf("dG = dGint + dGu_l + dGu_s\n");
    printf("%s %3d,%-3d : %3d,%-3d (%.2f = %.2f + %.2f + %.2f)\n", mfe.structure, Int->k, Int->i, Int->j, Int->l, G_min, Gi_min, Gul, Gus);
  }

  free(mfe.structure);
  for (i=1; i<=n1; i++) {
    free(p_c_S[i]);      
  }
  free(p_c_S);
  if(p_c2 != NULL) {
    for (i=1; i<=n2; i++) {
      free(p_c2_S[i]);      
    }
    free(p_c2_S);
  }
}
/* print the coordinates and the free energy for the region of highest accessibility */ 
PRIVATE void print_unstru(pu_contrib *p_c, int w) {
  int i,j,len,min_i,min_j;
  double dG_u, min_gu;
  
  if(p_c != NULL) {
    min_gu = 1000.0;
    len = p_c->length;

    for (i=1; i<=len; i++) {
      for (j=i; j < MIN((i+w),len+1); j++) {
	double blubb;
	if( (j-i+1) == w && i+w-1 <= len) {
	  blubb = p_c->H[i][j-i]+p_c->I[i][j-i]+p_c->M[i][j-i]+p_c->E[i][j-i];
	  
	  dG_u = -log(blubb)*(temperature+K0)*GASCONST/1000.0;
	  if(dG_u < min_gu ) {
	    min_gu = dG_u;
	    min_i=i;
	    min_j=i+w-1;
	  }
	}
      }
    }
    printf("%d,%d \t (%.3f) \t for u=%d\n",min_i,min_j,min_gu,w);
  } else {
    nrerror("error with prob unpaired");
  }
}
