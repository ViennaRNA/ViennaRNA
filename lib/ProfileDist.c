/*
	  Functions for handling the Base Pair Probability Matrix
		      Peter F Stadler, Ivo L Hofacker
			    Vienna RNA Package
*/ 
/* Last changed Time-stamp: <2002-11-07 12:46:16 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include  "dist_vars.h"
#include  "fold_vars.h"
#include  "part_func.h"
#include  "utils.h"
/*@unused@*/
static char rcsid[] = "$Id: ProfileDist.c,v 1.6 2002/11/07 11:49:59 ivo Exp $";

#define PUBLIC
#define PRIVATE        static
#define MIN(x,y)       (((x)<(y)) ? (x) : (y))      
#define MAX(x,y)       (((x)>(y)) ? (x) : (y))
#define MIN3(x,y,z)    (MIN(  (MIN((x),(y))) ,(z)))

PRIVATE int *alignment[2];

PUBLIC float    profile_edit_distance(const float *T1, const float *T2);
PUBLIC float    *Make_bp_profile(int length); 
PUBLIC void     free_profile(float *T);
PRIVATE void    sprint_aligned_bppm(const float *T1, const float *T2);
PRIVATE double  PrfEditCost(int i, int j, const float *T1, const float *T2);
PRIVATE double  average(double x, double y);

/*---------------------------------------------------------------------------*/

PUBLIC float profile_edit_distance(const float *T1, const float *T2)
{
  /* align the 2 probability profiles T1, T2 */
  /* This is like a Needleman-Wunsch alignment,
     we should really use affine gap-costs ala Gotoh */

  float    **distance;
  short    **i_point, **j_point;
  
  int           i, j, i1, j1, pos, length1,length2;
  float         minus, plus, change, temp;
  
  length1 = (int) T1[0];
  length2 = (int) T2[0];
  distance = (float **)     space((length1 +1)*sizeof(float *));
  if(edit_backtrack){
    i_point  = (short **)  space((length1 +1)*sizeof(short *));
    j_point  = (short **)  space((length1 +1)*sizeof(short *));
  }
  for(i=0; i<= length1; i++){
    distance[i] = (float *) space( (length2+1)*sizeof(float));
    if(edit_backtrack){
      i_point[i]  = (short *) space( (length2+1)*sizeof(short));
      j_point[i]  = (short *) space( (length2+1)*sizeof(short));
    }
  }
  
  for(i = 1; i <= length1; i++) {
    distance[i][0] = distance[i-1][0]+PrfEditCost(i,0,T1,T2);
    if(edit_backtrack){ i_point[i][0] = (short) i-1; j_point[i][0] = 0;   }
  }
  for(j = 1; j <= length2; j++) {
    distance[0][j] = distance[0][j-1]+PrfEditCost(0,j,T1,T2);
    if(edit_backtrack){ i_point[0][j] = 0;   j_point[0][j] = (short) j-1; }
  }    
  for (i = 1; i <= length1; i++) {
    for (j = 1; j <= length2 ; j++) {
      minus  = distance[i-1][j]  + PrfEditCost(i,0,T1,T2);
      plus   = distance[i][j-1]  + PrfEditCost(0,j,T1,T2);
      change = distance[i-1][j-1]+ PrfEditCost(i,j,T1,T2);
      
      distance[i][j] = MIN3(minus, plus, change);
      /* printf("%g ", distance[i][j]); */
      
      if(edit_backtrack){
	if(distance[i][j] == change) {
	  i_point[i][j]= (short)i-1; j_point[i][j]= (short) j-1;  }
	else if(distance[i][j] == plus) {
	  i_point[i][j]= (short)i  ; j_point[i][j]= (short)j-1;  }
	else {
	  i_point[i][j]= (short)i-1; j_point[i][j]= (short)j  ;  }
      }
    } 
    /* printf("\n"); */
  }
  /* printf("\n"); */
  temp = distance[length1][length2];
  for(i=0;i<=length1;i++) 
    free(distance[i]);
  free(distance);
  
  if(edit_backtrack){
    alignment[0] = (int *) space((length1+length2+1)*sizeof(int));
    alignment[1] = (int *) space((length1+length2+1)*sizeof(int));
    
    pos = length1+length2;
    i   = length1;
    j   = length2;
    while( (i>0)||(j>0) ) {
      i1 = i_point[i][j];
      j1 = j_point[i][j];
      if( ((i-i1)==1)&&((j-j1)==1) )  {  /* substitution    */ 
	alignment[0][pos] = i;
	alignment[1][pos] = j;
      }
      if( ((i-i1)==1)&&(j==j1) )      {  /* Deletion in [1] */
	alignment[0][pos] = i;
	alignment[1][pos] = 0;
      }
      if( (i==i1)&&((j-j1)==1)  )      {  /* Deletion in [0] */
	alignment[0][pos] = 0;
	alignment[1][pos] = j;
      }
      pos--;
      i = i1;
      j = j1;
    }
    for(i=pos+1; i<=length1+length2; i++){
      alignment[0][i-pos] = alignment[0][i];
      alignment[1][i-pos] = alignment[1][i];
    }
    alignment[0][0] = length1+length2-pos;   /* length of alignment */
    
    for(i=0; i<=length1; i++){
      free(i_point[i]); free(j_point[i]);
    }
    free(i_point); free(j_point);
    sprint_aligned_bppm(T1,T2);
    free(alignment[0]);
    free(alignment[1]);       
  }
  
  return temp;
}


/*---------------------------------------------------------------------------*/

PRIVATE double PrfEditCost(int i, int j, const float *T1, const float *T2)
{
  double  dist;
  int    k, kmax;

  kmax = (int) T1[1];
  if ((int) T2[1] != kmax) nrerror("inconsistent Profiles in PrfEditCost");
  
  if(i==0) {
    for(dist = 0. ,k=0 ; k<kmax ; k++) 
      dist += T2[j*kmax+k];
  }
  if(j==0) {
    for(dist = 0. ,k=0 ; k<kmax ; k++) 
      dist += T1[i*kmax+k];
  }
  if((i>0)&&(j>0)) {
    for(dist = 2.,k=0; k<kmax; k++)
      dist -= 2.*average(T1[i*kmax+k],T2[j*kmax+k]);
  }
  return dist;
}

/*---------------------------------------------------------------------------*/

PRIVATE double average(double x, double y)

/* can be essentially anything that fulfils :
   1.)     a(x,y)  =  a(y,x)
   2.)     a(x,y) >=  0       for 0<= x,y <= 1
   3.)     a(x,y) <=  (x+y)/2
   4.)     a(x,x) >=  a(x,y)  for 0<= x,y <= 1
   As in Bonhoeffer et al (1993) 'RNA Multi Structure Landscapes',
   Eur. Biophys. J. 22: 13-24 we have chosen  the geometric mean.
*/ 

{
    float a;
    a = (float) sqrt(x*y);    
    return a;
}

/*---------------------------------------------------------------------------*/

PUBLIC float *Make_bp_profile(int length)
{
   int i,j;
   int L=3;
   float *P; /* P[i*3+0] unpaired, P[i*3+1] upstream, P[i*3+2] downstream p */
   
   P =  (float *) space((length+1)*3*sizeof(float));
   /* indices start at 1 use first entries to store length and dimension */
   P[0] = (float) length;
   P[1] = (float) L;

   for( i=1; i<length; i++)
     for( j=i+1; j<=length; j++ ) {
       P[i*L+1] += pr[iindx[i]-j];
       P[j*L+2] += pr[iindx[i]-j];
     } 
   for( i=1; i<=length; i++)
     P[i*3+0] = 1 - P[i*3+1] - P[i*3+2];
   return (float *) P;
}
 
/*---------------------------------------------------------------------------*/
         
PRIVATE void sprint_aligned_bppm(const float *T1, const float *T2)
{
   int     i, length; 
   length = alignment[0][0];
   aligned_line[0] = (char *) space((length+1)*sizeof(char));
   aligned_line[1] = (char *) space((length+1)*sizeof(char));
   for(i=1; i<=length; i++){
      if(alignment[0][i] ==0) aligned_line[0][i-1] = '_';
      else { aligned_line[0][i-1] = bppm_symbol(T1+alignment[0][i]*3); }
      if(alignment[1][i] ==0) aligned_line[1][i-1] = '_';
      else { aligned_line[1][i-1] = bppm_symbol(T2+alignment[1][i]*3); }
   }
}

/*---------------------------------------------------------------------------*/

PUBLIC void print_bppm(const float *T)
{
   int i;
   for(i=1; i<=( (int)T[0]); i++)
      printf("%c",bppm_symbol(T+i*3));
   printf("\n");
}

/*---------------------------------------------------------------------------*/

PUBLIC void     free_profile(float *T)
{
   free(T);
}
   
