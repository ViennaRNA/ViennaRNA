/* Last changed Time-stamp: <95/07/12 19:49:53 ivo>*/
/*
	       Split Decomposition of Distance Matrices
	     as described by H.J.Bandelt and A.W.M.Dress
		       Adv Math, 92 (1992) p47
		 c Peter Stadler and Ivo L. Hofacker
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

#define  PUBLIC
#define  PRIVATE      static
#define  DEBUG        0    
#define  SHORT_OUTPUT 1    

#define  DINFTY       1.e32    
#define  ZERO         1.e-10


typedef struct {
   short   *splitlist[2];
   int      splitsize;
   double   isolation_index; } Split;
  
PUBLIC Split *split_decomposition(float **dist);
PUBLIC void free_Split(Split *x);
PUBLIC void print_Split(Split *x);
PUBLIC void sort_Split(Split *x);
 
PUBLIC Split *split_decomposition(float **dist)
{

   int      elm, n_of_splits;
   short   *full_slots;
   int      i,j,sp,new_sp,spp;
   int      i2,j1,j2,maxlen;
   int      x,y,z;
   double   alpha,beta,tmp;
   double   test1,test2;
   Split   *SD, *S;
   int number_of_points;
   
   number_of_points = (int) dist[0][0];
   
   /* Initialize */ 
   elm = 2;
   n_of_splits = 1;
   maxlen = (number_of_points+1)*(number_of_points+1);
   full_slots = (short *) space(2*(maxlen+1)*sizeof(short));
   full_slots[1]=1;
   
   SD = space((maxlen+1)*sizeof(Split));
   
   SD[0].splitlist[0]    = NULL;
   SD[0].splitlist[1]    = NULL; 
   SD[0].splitsize       = 1;
   SD[0].isolation_index = 0.0;
   /**/ 
   SD[1].splitlist[0]    = (short *) space((number_of_points+1)*sizeof(short));
   SD[1].splitlist[0][0] = number_of_points;
   SD[1].splitlist[0][1] = 1;
   SD[1].splitlist[1]    = (short *) space((number_of_points+1)*sizeof(short));
   SD[1].splitlist[1][1] = 2;
   SD[1].splitsize       = 1;
   SD[1].isolation_index = dist[1][2];


   /* Iteration */

   for( elm=3; elm <= number_of_points; elm++){

      for (sp=1, spp=n_of_splits;sp<=n_of_splits; sp++){
           
	 SD[sp].splitlist[0][0] = elm;

	 /* 1.  Split = {old_A | old_B+{elm}}  */

	 alpha = 2.0*SD[sp].isolation_index;
	 for(i2=1; i2<=(elm - SD[sp].splitsize); i2++){
	    if(i2==(elm-SD[sp].splitsize)) {
	       x = elm;
	    }
	    else{
	       x= SD[sp].splitlist[1][i2];
	    }
	    for(j1=1; j1 <= SD[sp].splitsize; j1++){
	       y= SD[sp].splitlist[0][j1];
	       for(j2=1; j2 <= (SD[sp].splitsize); j2++){
		  z= SD[sp].splitlist[0][j2];

		  /* calculate the value beta = beta(elm,x; y,z) */ 

		  beta = dist[elm][y] + dist[x][z];
		  tmp  = dist[elm][z] + dist[x][y];
		  if(tmp>beta) beta=tmp;
		  tmp  = dist[elm][x] + dist[y][z];
		  if(tmp>beta) beta=tmp;
		  beta -= ( dist[elm][x] + dist[y][z] );

		  if(beta<alpha) alpha=beta;
	       }
	    }
	 }
	 alpha/=2.;
           
	 if (alpha > ZERO){
	    /* add {old_A | old_B + {elm}} to the split-list */
	    spp++;
	    full_slots[spp] =1;
	    SD[spp].splitsize       = SD[sp].splitsize;
	    SD[spp].isolation_index = alpha;
	    SD[spp].splitlist[0]    = (short *) space((number_of_points+1)*sizeof(short));
	    SD[spp].splitlist[0][0] = elm;
	    SD[spp].splitlist[1]    = (short *) space((number_of_points+1)*sizeof(short));
	    for(i=1; i<=SD[spp].splitsize; i++)
	       SD[spp].splitlist[0][i] = SD[sp].splitlist[0][i];
	    for(i=1; i<=(elm-1-SD[spp].splitsize); i++)
	       SD[spp].splitlist[1][i] = SD[sp].splitlist[1][i];
	    SD[spp].splitlist[1][elm-SD[spp].splitsize] = elm;
	 }
    
	 /* 2. Split = {old_A+{elm} | old_B}   */

	 alpha = 2.0*SD[sp].isolation_index;
	 for(i2=1; i2<= (SD[sp].splitsize+1); i2++){
	    if(i2==(SD[sp].splitsize+1)){
	       x= elm;
	    }
	    else
	       {
		  x= SD[sp].splitlist[0][i2];
	       }
	    for(j1=1; j1 <= (elm-1-SD[sp].splitsize); j1++){
	       y= SD[sp].splitlist[1][j1];
	       for(j2=1; j2 <= (elm-1-SD[sp].splitsize); j2++){
		  z= SD[sp].splitlist[1][j2];

		  /* calculate the value beta = beta(elm,x; y,z) */ 

		  beta = dist[elm][y] + dist[x][z];
		  tmp  = dist[elm][z] + dist[x][y];
		  if(tmp>beta) beta=tmp;
		  tmp  = dist[elm][x] + dist[y][z];
		  if(tmp>beta) beta=tmp;
		  beta -= ( dist[elm][x] + dist[y][z] );

		  if(beta<alpha) alpha=beta;
	       }
	    }
	 }
	 alpha/=2.;

	 if (alpha > ZERO){
	    /* replace {old_A | old_B} by {old_A+{elm} | old_B} in the splitlist */
	    SD[sp].splitsize++;
	    SD[sp].splitlist[0][SD[sp].splitsize] = elm;
	    SD[sp].isolation_index = alpha;
	 }
	 else{
	    /* remove {old_A | old_B} from the split list */
	    full_slots[sp] = 0;
	    free(SD[sp].splitlist[0]);
	    SD[sp].splitlist[0]    = NULL;
	    free(SD[sp].splitlist[1]);
	    SD[sp].splitlist[1]    = NULL;
	    SD[sp].splitsize       = 0;
	    SD[sp].isolation_index = 0.0;
	 }
      }

      /* 3.  Split =  { {elm} | {1,...,elm-1} } */

      alpha=DINFTY;
      for(i=1;i<=elm-1;i++){
	 for(j=1; j<=elm-1;j++){
	    tmp = dist[elm][i]+dist[elm][j] - dist[i][j];
	    if( tmp < alpha) alpha = tmp;   
	 }
      }
      alpha/=2.;
      if (alpha > ZERO){
	 spp++;
	 full_slots[spp] = 1;
	 SD[spp].splitsize         = 1;
	 SD[spp].isolation_index   = alpha;
	 SD[spp].splitlist[0]      = (short *) space((number_of_points+1)*sizeof(short));
	 SD[spp].splitlist[0][0]   = elm;
	 SD[spp].splitlist[1]      = (short *) space((number_of_points+1)*sizeof(short));
	 SD[spp].splitlist[0][1]   = elm;
	 for (i=1; i<= elm-1; i++)
	    SD[spp].splitlist[1][i] = i;
      }
      /* note that spp now points to the last entry in the splitlist !!! */

      /* garbage collection .... unfortunately a bit boring */

      for ( sp=1, new_sp=0; sp<= spp; sp++ ){
	 if (full_slots[sp]){
	    new_sp++;
	    if(sp != new_sp){
	       /* copy all the junk from sp to new_sp */
	       full_slots[new_sp]=1; full_slots[sp]=0;
	       SD[new_sp].splitsize        = SD[sp].splitsize;
	       SD[new_sp].isolation_index  = SD[sp].isolation_index;
	       if(SD[new_sp].splitlist[0]==NULL) 
		  SD[new_sp].splitlist[0] = (short *) space((number_of_points+1)*sizeof(short));
	       if(SD[new_sp].splitlist[1]==NULL) 
		  SD[new_sp].splitlist[1] = (short *) space((number_of_points+1)*sizeof(short));
	       for(i=0;i<=elm;i++){
		  SD[new_sp].splitlist[0][i]  = SD[sp].splitlist[0][i];
		  SD[new_sp].splitlist[1][i]  = SD[sp].splitlist[1][i];
	       }
	       free(SD[sp].splitlist[0]);
	       SD[sp].splitlist[0] = NULL;
	       free(SD[sp].splitlist[1]);
	       SD[sp].splitlist[1] = NULL;
	       SD[sp].splitsize       = 0;
	       SD[sp].isolation_index = 0.0;
	    }
	 }
      }
      n_of_splits = new_sp;
      SD[0].splitsize = n_of_splits;
#if DEBUG
      for(test1=0, i=2; i<=elm; i++) for( j=1; j<i; j++) test1+=dist[i][j];
      for(test2=0, i=1; i<= n_of_splits; i++)
	 test2 += (elm-SD[i].splitsize)*SD[i].splitsize*SD[i].isolation_index;
      SD[0].isolation_index = (test1 - test2)/test1;
      print_Split(SD);
#endif
   } /* End of iteration */


   /* Calculate fraction of split-prime part for the full matrix
      if not already done */

   for(test1=0.0, i=2; i<=number_of_points; i++) 
      for( j=1; j<i; j++) test1+=dist[i][j];
   for(test2=0.0, i=1; i<= n_of_splits; i++) 
      test2 += (number_of_points-SD[i].splitsize)*
	 SD[i].splitsize*SD[i].isolation_index;
   SD[0].isolation_index = (test1 - test2)/test1;
   
#if DEBUG
   print_Split(SD);
#endif
   /* free superfluous space */ 

   S = space((n_of_splits+1)*sizeof(Split));
   for (i=0; i<= n_of_splits; i++) S[i] = SD[i];
   free(SD);

   if(DEBUG)
      print_Split(S);

   return(S);
}  

/* -------------------------------------------------------------------------- */

PUBLIC void free_Split(Split *x)
{
   int i;
   for (i=0; i<=x[0].splitsize; i++) {
      if( x[i].splitlist[0] != NULL ) free( x[i].splitlist[0] ); 
      if( x[i].splitlist[1] != NULL ) free( x[i].splitlist[1] );
   }
   free(x);
}

/* -------------------------------------------------------------------------- */

PUBLIC void print_Split(Split *x) 
{
   int j,k;
   printf("> %d Split Decomposition",x[0].splitsize);
   for (k=1; k<=x[0].splitsize; k++) {
      if( (x[k].splitlist[0])&&(x[k].splitlist[1]) ) {
         printf("\n%3d   %8.4f  : {",k,x[k].isolation_index);
         for(j=1;j<=x[k].splitsize;j++)
            printf(" %3d",x[k].splitlist[0][j]);
         printf("   |");
#if SHORT_OUTPUT
	 printf(" ...");
#else
	 for(j=1;j<=(x[k].splitlist[0][0]-x[k].splitsize);j++)
	    printf(" %3d",x[k].splitlist[1][j]);
	 printf(" } ");
#endif
      }
   }
   printf("\n      %8.4f  : { [Split prime fraction] }\n", x[0].isolation_index);
}

    
/* -------------------------------------------------------------------------- */
           
int CompareSplit(const void *v1 , const void *v2)
{
   Split *S1, *S2;

   S1 = (Split *) v1; S2 = (Split *) v2;
   if(S1->isolation_index < S2->isolation_index) return +1;
   if(S2->isolation_index < S1->isolation_index) return -1;
   return 0;
}
   
/* -------------------------------------------------------------------------- */

PUBLIC void sort_Split(Split *x) 
{
   int i,j;
   int t1;

   qsort((void *) &(x[1]), x[0].splitsize ,sizeof(Split),CompareSplit );
  
   for(i=1; i<=x[0].splitsize; i++) {
      if(x[i].splitsize > x[i].splitlist[0][0]-x[i].splitsize) {
         for(j=1; j<=x[i].splitsize; j++) {
            t1 = x[i].splitlist[0][j];
            x[i].splitlist[0][j] = x[i].splitlist[1][j];
            x[i].splitlist[1][j] = t1;
         }
         x[i].splitsize = x[i].splitlist[0][0]-x[i].splitsize;
      }
   }
}


/* -------------------------------------------------------------------------- */
