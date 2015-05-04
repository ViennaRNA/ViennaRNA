/*
                Distances of Secondary Structure Ensembles
          Peter F Stadler, Ivo L Hofacker, Sebastian Bonhoeffer
                        Vienna RNA Package
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include "part_func.h"
#include "fold_vars.h"
#include "dist_vars.h"
#include "utils.h"
#include "PS_dot.h"
#include "read_epars.h"
#include "profiledist.h"
#include "RNApdist_cmdl.h"


#define MAXLENGTH  10000
#define MAXSEQ      1000
/*@unused@*/
static char rcsid[] = "$Id: RNApdist.c,v 1.8 2002/11/07 12:19:41 ivo Exp $";

PRIVATE void command_line(int argc, char *argv[]);
PRIVATE void usage(void);
PRIVATE void print_aligned_lines(FILE *somewhere);

PRIVATE char task;
PRIVATE char outfile[FILENAME_MAX_LENGTH];
PRIVATE char  ruler[] ="....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";
static int noconv = 0;

int main(int argc, char *argv[])

{
  float     *T[MAXSEQ];
  int        i,j, istty, n=0;
  int        type, length, taxa_list=0;
  float      dist;
  FILE      *somewhere=NULL;
  char      *structure;
  char      *line=NULL, fname[FILENAME_MAX_LENGTH], *list_title=NULL;
  plist     *pr_pl, *mfe_pl;

  pr_pl = mfe_pl = NULL;

  command_line(argc, argv);

  if((outfile[0]=='\0')&&(task=='m')&&(edit_backtrack))
    strcpy(outfile,"backtrack.file");
  if (outfile[0]!='\0') somewhere = fopen(outfile,"w");
  if (somewhere==NULL) somewhere = stdout;
  istty   = (isatty(fileno(stdout))&&isatty(fileno(stdin)));

  while (1) {
    if ((istty)&&(n==0)) {
      printf("\nInput sequence;  @ to quit\n");
      printf("%s\n", ruler);
    }

    type = 0;
    do {  /* get sequence to fold */
      if (line!=NULL) free(line);
      *fname='\0';
      if ((line=get_line(stdin))==NULL) {type = 999; break;}
      if (line[0]=='@') type = 999;
      if (line[0]=='*') {
        if (taxa_list==0) {
          if (task=='m') taxa_list=1;
          printf("%s\n", line);
          type = 0;
        } else {
          list_title = strdup(line);
          type = 888;
        }
      }
      if (line[0]=='>') {
        if (sscanf(line,">%" XSTR(FILENAME_ID_LENGTH) "s", fname)!=0)
          strcat(fname, "_dp.ps");
        if (taxa_list)
          printf("%d : %s\n", n+1, line+1);
        else printf("%s\n",line);
        type = 0;
      }
      if (isalpha(line[0]))  {
        char *cp;
        cp =strchr(line,' ');
        if (cp) *cp='\0';
        type = 1;
      }
    } while(type==0);

    if( (task == 'm')&&(type>800) ) {
      if (taxa_list)
        printf("* END of taxa list\n");
      printf("> p %d (pdist)\n",n);
      for (i=1; i<n; i++) {
        for (j=0; j<i; j++) {
          printf("%g ",profile_edit_distance(T[i], T[j]));
          if(edit_backtrack) fprintf(somewhere,"> %d %d\n",i+1,j+1);
          print_aligned_lines(somewhere);
        }
        printf("\n");
      }
      if (type==888) {  /* do another distance matrix */
        n = 0;
        printf("%s\n", list_title);
        free(list_title);
      }
    }

    if(type>800) {
      for (i=0; i<n; i++)
        free_profile(T[i]);
      if (type == 888) continue;
      if (outfile[0]!='\0') (void) fclose(somewhere);
      if (line!= NULL) free(line);
      return 0; /* finito */
    }

    length = (int) strlen(line);
    for (i=0; i<length; i++) {
      line[i]=toupper(line[i]);
      if (!noconv && line[i] == 'T') line[i] = 'U';
    }

    /* init_pf_fold(length); <- obsolete */
    structure = (char *) space((length+1)*sizeof(char));
    (void) pf_fold(line,structure);

    if (*fname=='\0')
      sprintf(fname, "%d_dp.ps", n+1);

    /* PS_dot_plot(line, fname); <- NOT THREADSAFE and obsolete function! */

    /* get pairlist of probability matrix */
    assign_plist_from_pr(&pr_pl, pr, length, 1e-5);
    /* no previous mfe call thus no mfe structure information known */
    mfe_pl = (plist *)space(sizeof(plist));
    mfe_pl[0].i = mfe_pl[0].j = 0;

    /* call threadsafe dot plot printing function */
    PS_dot_plot_list(line, fname, pr_pl, mfe_pl, "");

    T[n] = Make_bp_profile_bppm(pr, length);
    if((istty)&&(task=='m')) printf("%s\n",structure);
    free(structure);
    free(mfe_pl);
    free(pr_pl);
    free_pf_arrays();

    n++;
    switch (task) {
    case 'p' :
      if (n==2) {
        dist = profile_edit_distance(T[0],T[1]);
        printf("%g\n",dist);
        print_aligned_lines(somewhere);
        free_profile(T[0]);
        free_profile(T[1]);
        n=0;
      }
      break;
    case 'f' :
      if (n>1) {
        dist = profile_edit_distance(T[1], T[0]);
        printf("%g\n",dist);
        print_aligned_lines(somewhere);
        free_profile(T[1]);
        n=1;
      }
      break;
    case 'c' :
      if (n>1) {
        dist = profile_edit_distance(T[1], T[0]);
        printf("%g\n",dist);
        print_aligned_lines(somewhere);
        free_profile(T[0]);
        T[0] = T[1];
        n=1;
      }
      break;

    case 'm' :
      break;

    default :
      nrerror("This can't happen.");
    }    /* END switch task */
    (void) fflush(stdout);
  }    /* END while */
  if (line !=NULL) free(line);
  return 0;
}

/* ----------------------------------------------------------------- */

PRIVATE void command_line(int argc, char *argv[])
{

  int i, sym;
  char *ns_bases=NULL, *c;
  char *ParamFile=NULL;
  struct  RNApdist_args_info args_info;

  task = 'p';

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNApdist_cmdline_parser (argc, argv, &args_info) != 0) exit(1);

  /* temperature */
  if(args_info.temp_given)
    temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)
    tetra_loop=0;

  /* set dangle model */
  if(args_info.dangles_given){
    dangles = args_info.dangles_arg;
    if(dangles) dangles = 2;
  }

  /* set energy model */
  if(args_info.energyModel_given)
    energy_set = args_info.energyModel_arg;

  /* do not allow weak pairs */
  if(args_info.noLP_given)
    noLonelyPairs = 1;

  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)
    noGU = 1;

  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given)
    no_closingGU = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)
    noconv = 1;

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* take another energy parameter set */
  if(args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  if(args_info.compare_given){
    switch(args_info.compare_arg[0]){
      case 'p':
      case 'm':
      case 'f':
      case 'c': task=args_info.compare_arg[0];
                break;
      default:  RNApdist_cmdline_parser_print_help();
                exit(EXIT_FAILURE);
    }
  }

  if(args_info.backtrack_given){
    if(strcmp(args_info.backtrack_arg, "none")){
      strncpy(outfile, args_info.backtrack_arg, FILENAME_MAX_LENGTH-1);
    }
    edit_backtrack = 1;
  }

  /* free allocated memory of command line data structure */
  RNApdist_cmdline_parser_free (&args_info);

  /* do some preparations */
  if (ParamFile!=NULL)
    read_parameter_file(ParamFile);

  if (ns_bases!=NULL) {
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

}

/* ------------------------------------------------------------------------- */

PRIVATE void usage(void)
{
  nrerror("usage: RNApdist [-Xpmfc] [-B [file]] [-T temp] [-4] [-d] [-noGU]\n"
          "                [-noCloseGU] [-noLP] [-e e_set] [-P paramfile] [-nsp pairs]");
}

/*--------------------------------------------------------------------------*/

PRIVATE void print_aligned_lines(FILE *somewhere)
{
  if (edit_backtrack)
    fprintf(somewhere, "%s\n%s\n", aligned_line[0], aligned_line[1]);
}

/*--------------------------------------------------------------------------*/
