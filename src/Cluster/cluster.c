/*
		 Cluster Analysis using Ward's Method
		Ward J Amer Stat Ass, 58 (1963), p236
		   c Peter Stadler and Ivo Hofacker
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ViennaRNA/utils.h"

#define PUBLIC
#define PRIVATE static

#define INFINITY   1000000

typedef struct{
        int   set1;
        int   set2;
        float distance;
        float distance2;
        } Union;

typedef struct {
                 int  type; 
                 int  weight;
                 int  father;
                 int  sons;
                 int  leftmostleaf;                 
               } Postorder_list;

PUBLIC Union *wards_cluster(float **clmat);
PUBLIC Union *neighbour_joining(float **clmat);
PUBLIC void   printf_phylogeny(Union *tree, char *type);
       

/*--------------------------------------------------------------------*/

PUBLIC Union *wards_cluster(float **clmat)
{
   float  **d;
   int      *indic;
   int      *size;
   float    *help;
   Union    *tree;

   float    min,deno,xa,xb,x;
   int      i,j,step,s=0,t=0,n;

   n= (int)(clmat[0][0]);

   size  = (int *)     vrna_alloc((n+1)*sizeof(int));
   d     = (float **)  vrna_alloc((n+1)*sizeof(float *));
   for(i=0;i<=n;i++)
      d[i] = (float *) vrna_alloc((n+1)*sizeof(float));
   indic = (int *)     vrna_alloc((n+1)*sizeof(int));
   help  = (float *)   vrna_alloc((n+1)*sizeof(float));
   tree  = (Union *)   vrna_alloc((n+1)*sizeof(Union));

   tree[0].set1      = n;
   tree[0].set2      = 0;
   tree[0].distance  = 0.0;   
   tree[0].distance2 = 0.0;    

   for (i=1;i<=n;i++) size[i]=1;
   for (i=1; i<=n; i++){
      for(j=1; j<=n; j++){
         d[i][j] = clmat[i][j];
      }
   }

   /* look for the indices [s,t]  with minimum  d[s][t]*/
   for(step=1;step<n; step++){
      min = INFINITY;
      for (i=1; i<=n; i++){
         if (indic[i]==0){
            for (j=1; j<=n; j++){
               if(j!=i){
                  if (indic[j]==0){
                     if(d[i][j] < min) {
                        min = d[i][j];
                        s = i;
                        t = j;
                     }
                  }
               }
            }
         }  
      }

      /* now we have to join the clusters s and t and update the working array*/

      tree[step].set1     = s;
      tree[step].set2     = t;
      tree[step].distance = min;
      indic[t] =1;
   
      for (i=1; i<=n; i++){
         if (indic[i]==0){
            deno = (float) (size[i]+size[s]+size[t]);
            xa = ((float) (size[i]+size[s]))/deno; 
            xb = ((float) (size[i]+size[t]))/deno;
             x = ((float) size[i])/deno;
            help[i] = xa*d[i][s] + xb*d[i][t] - x*d[s][t];
         }
      }
      for (i=1; i<=n; i++){
         if (indic[i]==0){
            d[s][i] = help[i];
            d[i][s] = help[i];
         }
      }
      size[s] += size[t];
   }

   free(help);
   free(indic);
   for(i=0;i<=n;i++) free(d[i]);
   free(d);
   free(size);
 
   return tree;
}        

/*--------------------------------------------------------------------*/

PUBLIC Union *neighbour_joining(float **clmat)
{            
  int n,i,j,k,l,step,ll[3];
  float b1,b2,b3,nn,tot,tmin,d1,d2;
  int mini=0, minj=0;
  int    *indic;
  float  *av, *temp;
  float **d;
  Union   *tree;

  n = (int) (clmat[0][0]);

  tree = (Union *) vrna_alloc((n+1)*sizeof(Union));
  indic = (int   *) vrna_alloc((n+1)*sizeof(int)   );
  av   = (float *) vrna_alloc((n+1)*sizeof(float) );
  temp = (float *) vrna_alloc((n+1)*sizeof(float) );
  d    = (float**) vrna_alloc((n+1)*sizeof(float*));
  for(i=0;i<=n;i++) d[i] = (float*) vrna_alloc((n+1)*sizeof(float));


  for (i=1; i<=n; i++){
     for(j=1; j<=n; j++){
        d[i][j] = clmat[i][j];
     }
     d[i][i]=0.0;
  }

  tree[0].set1      = n;
  tree[0].set2      = 0;
  tree[0].distance  = 0.0; 
  tree[0].distance2 = 0.0;
  
  nn = (float) n;

  for(step=1;step<=n-3;step++) {
/*    for(j=2;j<=n;j++) for(i=1;i<j;i++) d[j][i]=d[i][j]; */
     tmin=99999.9;
     for(l=2; l<=n; l++){
        if(!indic[l]) {
           for(k=1;k<l;k++){                          
              if(!indic[k]) {                                       
                 d1=0.0; d2=0.0;
                 for(i=1; i<=n;i++) {
                    d1 += d[i][k];
                    d2 += d[i][l];
                 }
                 tot=(nn-2.0)*d[k][l]-d1-d2;
                 if(tot<tmin){
                    tmin=tot;
                    mini=k;
                    minj=l;
                 }
              }
           }
        }
     }

     d1=0.0; d2=0.0;                                                                
     for(i=1;i<=n;i++) {
        d1 +=d[i][mini];
        d2 +=d[i][minj];
     }
     d1 = (d1-d[mini][minj])/(nn-2.0);
     d2 = (d2-d[mini][minj])/(nn-2.0);

     tree[step].set1      = mini;
     tree[step].distance  = (d[mini][minj]+d1-d2)*0.5-av[mini];
     tree[step].set2      = minj; 
     tree[step].distance2 = d[mini][minj]-(d[mini][minj]+d1-d2)*0.5-av[minj];

     av[mini]=d[mini][minj]*0.5;

     nn=nn-1.0;
     indic[minj]=1;
     for(j=1;j<=n;j++) { 
        if(!indic[j]) 
           temp[j]=(d[mini][j]+d[minj][j])*0.5;
     }
     for(j=1;j<=n;j++) {
        if(!indic[j]) {
           d[mini][j] = temp[j];
           d[j][mini] = temp[j];
        }
     }                               
     for(j=1;j<=n;j++){
         d[minj][j] = 0.0;
         d[j][minj] = 0.0;
         d[j][j]    = 0.0;
     }
  }  
                                            
  j=0;   
  for(i=1;i<=n;i++) {
     if(!indic[i]){
        ll[j]=i;
        j++;
     }
  }          
  b1=(d[ll[0]][ll[1]]+d[ll[0]][ll[2]]-d[ll[1]][ll[2]])*0.5;
  b2=d[ll[0]][ll[1]]-b1;
  b3=d[ll[0]][ll[2]]-b1;
  b1 -= av[ll[0]];
  b2 -= av[ll[1]];
  b3 -= av[ll[2]];
  tree[step].set1      = ll[1];
  tree[step].distance  = b2;
  tree[step].set2      = ll[2];
  tree[step].distance2 = b3;
  step++;
  tree[step].set1      = ll[0];
  tree[step].distance  = 0.0;
  tree[step].set2      = ll[1];
  tree[step].distance2 = b1;

  for(i=0;i<=n;i++) free(d[i]); free(d); 
  free(temp);
  free(av);
  free(indic);

  return tree;
}

/*--------------------------------------------------------------------*/


PUBLIC void   printf_phylogeny(Union *tree, char *type)
{
   int i,n;
   n = tree[0].set1;
   
   printf("> %d %s ( Phylogeny using ",n, type);
   switch(type[0]){
    case 'W'  : 
      printf("Ward's Method )\n");
      printf("> Nodes      Variance\n");
      for(i=1; i<n; i++)
         printf("%3d %3d    %9.4f \n", 
         tree[i].set1, tree[i].set2, tree[i].distance);
      break;
    case 'N' :
      printf("Saitou's Neighbour Joining Method )\n");
      printf("> Nodes      Branch Length in Tree\n");
      for(i=1; i<n; i++)
         printf("%3d %3d    %9.4f   %9.4f\n", 
         tree[i].set1, tree[i].set2, tree[i].distance, tree[i].distance2);
      break;
    default:
      printf(" non-identified Method )\n");
      for(i=1; i<n; i++)
         printf("%3d %3d    %9.4f   %9.4f\n", 
         tree[i].set1, tree[i].set2, tree[i].distance, tree[i].distance2);
   }
}

/*--------------------------------------------------------------------*/
