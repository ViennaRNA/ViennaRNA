/* Last changed Time-stamp: <2005-07-23 16:50:24 ivo> */
/*                
             Compute duplex structure of two RNA strands

                           c Ivo L Hofacker
                          Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "read_epars.h"
#include "subopt.h"
#include "duplex.h"

PRIVATE void  print_struc(duplexT const *dup);
/*@unused@*/
static char rcsid[] = "$Id: RNAduplex.c,v 1.5 2007/08/26 09:41:12 ivo Exp $";

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

static void usage(void);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[]){
  char *s1, *s2, *line;
  char  fname[13], *input_string;
  unsigned int  input_type;
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, l, sym, r;
  double deltaf;
  double kT, sfact=1.07;
  int   pf=0, istty, delta=-1;
  int noconv=0;
   
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') 
      switch ( argv[i][1] )
        {
        case 'T':  if (argv[i][2]!='\0') usage();
          if(i==argc-1) usage();
          r=sscanf(argv[++i], "%lf", &temperature);
          if (!r) usage();
          break;
        case 'p':
          fprintf(stderr, "partition function folding not yet implemented\n");
          usage();
          pf=1;
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
          if (i>=argc-1) usage();
          r=sscanf(argv[++i], "%lf", &deltaf);
          if (r!=1) usage();
          delta = (int) (0.1+deltaf*100);
          break;
        case 's': subopt_sorted=1;
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
        default: usage();
        } 
  }

  /*
  #############################################
  # begin initializing
  #############################################
  */
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

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do{
    duplexT mfe, *subopt;
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty) print_tty_input_seq_str("Input two sequences (one line each)");

    /* extract filename from fasta header if available */
    fname[0] = '\0';
    while((input_type = get_input_line(&input_string, (istty) ? VRNA_INPUT_NOPRINT : 0)) == VRNA_INPUT_FASTA_HEADER){
      (void) sscanf(input_string, "%13s", fname);
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      s1 = strdup(input_string);
      free(input_string);
    }

    /* get second sequence */
    input_type = get_input_line(&input_string, (istty) ? VRNA_INPUT_NOPRINT : 0);
    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      s2 = strdup(input_string);
      free(input_string);
    }

    if(noconv){
      str_RNA2RNA(s1);
      str_RNA2RNA(s2);
    }
    else{
      str_DNA2RNA(s1);
      str_DNA2RNA(s2);
    }

    if (istty) printf("lengths = %d,%d\n", strlen(s1), strlen(s2));

    /*
    ########################################################
    # done with 'stdin' handling, now init everything properly
    ########################################################
    */

     /*
    ########################################################
    # begin actual computations
    ########################################################
    */
    update_fold_params();

    if (delta>=0) {
      duplexT *sub;
      subopt = duplex_subopt(s1, s2, delta, 5);
      for (sub=subopt; sub->i >0; sub++) {
        print_struc(sub);
        free(sub->structure);
      }
      free(subopt);
    }
    else {
      mfe = duplexfold(s1, s2);
      print_struc(&mfe);
      free(mfe.structure);
    }
    (void) fflush(stdout);
    free(s1);
    free(s2);
  } while (1);
  return 0;
}

static void print_struc(duplexT const *dup) {
  int l1;
  l1 = strchr(dup->structure, '&')-dup->structure;
  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, dup->i+1-l1,
         dup->i, dup->j, dup->j+strlen(dup->structure)-l1-2, dup->energy);
}
    
static void usage(void)
{
  nrerror("usage:\n"
          "RNAduplex [-e range] [-s]\n"
          "          [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n" 
          "          [-noLP] [-P paramfile] [-nsp pairs] [-noconv]\n");
}
