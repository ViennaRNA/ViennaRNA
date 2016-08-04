
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "split.h"
#include "cluster.h"
#include "distance_matrix.h"
#include "treeplot.h"
#include "ViennaRNA/utils.h"


#define PUBLIC
#define PRIVATE   static

PRIVATE void usage(void);

int main(int argc, char *argv[])
{
   int     i,j;
   float **dm;
   Split  *S;
   Union  *U;
   char    type[5];

   short   Do_Split=1, Do_Wards=0, Do_Nj=0;

   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') {
	 switch ( argv[i][1] ) {
	  case 'X':  if (argv[i][2]=='\0') { Do_Split = 1 ; break; }
	    Do_Split = 0;
	    Do_Wards = 0;
	    Do_Nj    = 0;
	    for(j=2;j<strlen(argv[i]);j++) {
	       switch(argv[i][j]) {
		case 's' :  Do_Split = 1; 
		  break;
		case 'w' :  Do_Wards = 1;
		  break;
		case 'n' :  Do_Nj    = 1;
		  break;
		  default :
		  usage();
	       }
	    }
	    break;
	    default : 
	    usage();
         }
      }
   }

   while ((dm=read_distance_matrix(type))!=NULL) {

      printf_taxa_list();
      printf("> %s\n",type);
      
      if(Do_Split) {
         S = split_decomposition(dm);
         sort_Split(S);
         print_Split(S);
         free_Split(S);
      }
      if(Do_Wards) {
         U = wards_cluster(dm);

         printf_phylogeny(U,"W");
         PSplot_phylogeny(U,"wards.ps","Ward's Method");
         free(U);
      }
      if(Do_Nj) {
         U = neighbour_joining(dm);
         printf_phylogeny(U,"Nj");
         PSplot_phylogeny(U,"nj.ps","Neighbor Joining");
         free(U);
      }
      free_distance_matrix(dm);
   }
   return 0;
}


PRIVATE void usage(void)
{
   vrna_message_error("usage: AnalyseDist [-X[swn]]");
   exit(0);
}
