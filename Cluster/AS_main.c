#include <stdio.h>
#include <string.h>
#include "distance_matrix.h"
#include "statgeom.h"
#include "split.h"
#include "cluster.h"
#include "treeplot.h"
#include "utils.h"


#define PUBLIC
#define PRIVATE   static

PRIVATE char  scale1[] = "....,....1....,....2....,....3....,....4";
PRIVATE char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);

main(int argc, char *argv[])
{
   int      n,i,j,l;
   int      intty;
   int      outtty;
   char    *mask, junk[20];
   char   **s;
   char   **ss[4];
   float   *B;
   float  **dm;
   Split   *S;
   Union   *U;
   char     DistAlgorithm='H';
   int      nn[4];
   short    Do_Split=0, Do_Wards=0, Do_Stg=1, Do_4_Stg=0, Do_Nj=0, Do_Mat=0;
   float    per_digit, per_gap;
   
   mask   = space(sizeof(char)*54);
   strcpy (mask,"%ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') {
	 switch ( argv[i][1] ) {
	  case 'X':  
	     if (argv[i][2]=='\0') { Do_Stg = 1 ; break; }
             Do_Split = 0;
             Do_Wards = 0;
             Do_Stg   = 0;
             Do_Nj    = 0;
	     for(j=2;j<strlen(argv[i]);j++) {
                switch(argv[i][j]) {
                 case 's' :  Do_Split = 1; 
                    break;
                 case 'w' :  Do_Wards = 1;
                    break;
                 case 'b' :  Do_Stg   = 1;
                    break;
                 case 'n' :  Do_Nj    = 1;
                    break;
		 case 'm' :  Do_Mat   = 1;
		   break;
                 default :
                    usage();
                }
             }
             break;
          case 'Q': 
             Do_4_Stg = 1;
             break;
          case 'M':
             if(mask) { free(mask); mask = NULL; }
             switch (argv[i][2] ) {
              case '\0' :
                usage();
                break;
              case 'a' : 
                mask   = space(sizeof(char)*54);
                strcpy(mask,
		"%ABCDEFGHIJKLMNOPQRSTUVWabcdefghijklmnopqrstuvwxyz");
                if(argv[i][3]=='+') mask[0] = '~';  /* make case sensitive */
                break;
              case 'u' :
                mask   = space(sizeof(char)*28);
                strcpy(mask,"~ABCDEFGHIJKLMNOPQRSTUVW");
                break;
              case 'l' :
                mask   = space(sizeof(char)*28);
                strcpy(mask,"~abcdefghijklmnopqrstuvwxyz");
                break;
              case 'c' :
                mask   = space(sizeof(char)*12);
                strcpy(mask,"~1234567890");
                break;
              case 'n' :
                mask   = space(sizeof(char)*64);
                strcpy (mask,
                "%ABCDEFGHIJKLMNOPQRSTUVWabcdefghijklmnopqrstuvwxyz1234567890");
                if(argv[i][3]=='+') mask[0] = '~';  /* make case sensitive */
                break;
              case 'R' :     /* RNA */
                mask   = space(sizeof(char)*10);
                strcpy(mask,"%GCAUgcau");
                if(argv[i][3]=='+') mask[0] = '~';  /* make case sensitive */
                break;
              case 'D' :     /* DNA */
                mask   = space(sizeof(char)*10);
                strcpy (mask,"%GCATgcat");
                if(argv[i][3]=='+') mask[0] = '~';  /* make case sensitive */
                break;
              case 'A' :    /*  AMINOACIDS */
                mask   = space(sizeof(char)*42);
                strcpy(mask,"%ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy");
                if(argv[i][3]=='+') mask[0] = '~';  /* make case sensitive */
                break;
              case 'S' :    /* SECONDARY STRUCTURES */
                mask   = space(sizeof(char)*6);
                strcpy(mask,"~().^");
                break;
              case '%' :    /* ARBITRARY ALPHABETS */
                l = strlen(argv[i]);
                if(argv[i][l] == '+'){
                   mask =   space(sizeof(char)*(l-2));
                   mask[0] = '~';
                   for(j=1;j<=l-4;j++) mask[j] = argv[i][j+2];
		   mask[l-3]='\0';
                }                
                if(argv[i][l] == '!'){
                   mask =   space(sizeof(char)*(l-2));
                   mask[0] = '!';
                   for(j=1;j<=l-4;j++) mask[j] = argv[i][j+2];
		   mask[l-3]='\0';
                }
                else { 
                   mask =   space(sizeof(char)*(l-1));
                   mask[0] = '%';
                   for(j=1;j<=l-3;j++) mask[j] = argv[i][j+2];
		   mask[l-2]='\0';
                }
                break;
              default : 
                usage();
            }   
            break;
          case 'D' :               /* choose algorithm */
            switch(argv[i][2]) {
              case '\0' : 
                usage();
                break;
              case 'H' :
                DistAlgorithm = 'H';
                break;
              case 'A' :
                DistAlgorithm = 'A';
                if(argv[i][3]==',') {
                   per_digit=-1.;
                   sscanf(argv[i],"%[^,],%g",junk,&per_digit);
                   if(per_digit<0.) usage();
                   per_gap = per_digit;
                   Set_StrEdit_GapCosts(per_digit,per_gap);
                }
                break;
              case 'G' :
                DistAlgorithm = 'G';
                if(argv[i][3]==',') {
                   per_digit=-1.;
                   per_gap  =-1.;
                   sscanf(argv[i]+4,"%f,%f",&per_digit,&per_gap);
                   if((per_digit<0.)||(per_gap<0)) usage(); 
                   Set_StrEdit_GapCosts(per_digit,per_gap);
                }
                break;
              default:
                usage();
             }
             break;
           case 'd' :               /* choose distance matrix */  
             switch(argv[i][2]) {
               case 'D' :    /* Dayhoff Distances */
	       case 'A' :    /* Aminoacid Distance (Hofacker & Borstnik) */
               case 'B' :    /* RY distances for nucleotides */
               case 'H' :    /* Hogeweg's Distance for Secondary Structures */
               case 'S' :    /* Simple Distance (superfluous option) */
                 Set_StrEdit_CostMatrix(argv[i][2]);
                 break;
               default :
                 usage();
             }
             break;
          default : 
             usage();
         }
      }
   }

   /* END PARSING OF COMMAND LINE */

   intty = isatty(fileno(stdin));
   outtty= isatty(fileno(stdout));
   
   if(intty){
      if(outtty) {
         printf("Input sequences; @ to mark end of input\n");
         printf("%s%s\n", scale1, scale2);
      }
      else {
         fprintf(stderr,"Input sequences; @ to mark end of input\n");
         fprintf(stderr,"%s%s\n", scale1, scale2);
      }
   }
   
   while ((s=read_sequence_list(&n,mask))!=NULL) {
      dm = NULL;
      if(Do_4_Stg) {
         ss[0] = s;
         nn[0] = n;
         for(i=1;i<4;i++) {
            ss[i] = read_sequence_list(&n,mask);
            if(ss[i]==NULL) nrerror("read_sequences: wrong or insufficient data.");
            nn[i] = n;
         }
         printf_taxa_list();
         B = statgeom4(ss,nn);
         printf_stg(B);
         SimplifiedBox(B,"box.ps");    /* This is preliminary !!! */ 
         free(B);
         
         for(i=0;i<4;i++){
	    for(j=0;j<nn[i];j++) free(ss[i][j]);
	    free(ss[i]);
         }
	 free(ss);
      }
      else {
         printf_taxa_list();
         if(Do_Stg) {
            B = statgeom(s,n);
	    if (B) {
	       printf_stg(B);
	       SimplifiedBox(B,"box.ps");
	       free(B);
	    }
         }
         if((Do_Split)||(Do_Wards)||(Do_Nj)||(Do_Mat)) {
            switch(DistAlgorithm) {
             case 'H' :
		dm = Hamming_Distance_Matrix(s,n);
		printf("> %s\n","H (Hamming Distance)");
		break;
	      case 'A' :
		dm = StrEdit_SimpleDistMatrix(s,n);
		printf("> %s\n","A (Needleman-Wunsch Distance)");
		break;
	      case 'G' :
		 dm = StrEdit_GotohDistMatrix(s,n);
		printf("> %s\n","G (Gotoh Distance)");
		break;
	      default:
		nrerror("This can't happen.");
	     }
         }
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
	 if(Do_Mat) printf_distance_matrix(dm);
      }

      if (dm!=NULL) free_distance_matrix(dm);
      for(i=0;i<n;i++) free(s[i]);
      free(s);
   }
   return 0;
}

/* -----------------------------------------------------------------------*/   

PRIVATE void usage(void)
{
   nrerror("usage: AnalyseSeqs [-X[bswnm]] [-Q] [-M{mask}] \n"
   "                   [-D{H|A[,cost]|G[,cost1,cost2]}] [-d{D|B|H|S}]");
}
