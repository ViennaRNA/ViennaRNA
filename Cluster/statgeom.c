/*      Quadruple Statistics on an arbitrary alphabet    

              Vienna RNA Package   ---  Peter F Stadler   1993

*/

#include <stdio.h>
#include <strings.h>
#include <ctype.h>
#include "utils.h"
#include "PS3D.h"
#include "distance_matrix.h"

#define PUBLIC    
#define PRIVATE    static

#define MIN2(A, B)         ((A) < (B) ? (A) : (B))
#define MIN4(A, B, C, D)   MIN2( (MIN2((A),(B))), (MIN2((C),(D))) )
#define MAX2(A, B)         ((A) > (B) ? (A) : (B))
#define MAX4(A, B, C, D)   MAX2( (MAX2((A),(B))), (MAX2((C),(D))) )

PRIVATE int IBox[16];


PUBLIC float *statgeom(char **seqs, int n_of_seqs);
PUBLIC float *statgeom4(char **ss[4], int nn[4]);
PUBLIC void   printf_stg(float *B);

PRIVATE void SingleBox(char *x1, char *x2, char *x3, char *x4);
PRIVATE void SortSingleBox(void);

/* ----------------------------------------------------------------------- */

PUBLIC float *statgeom(char **seqs, int n_of_seqs)
{
   int i,j,k,l;
   int i1;

   float *B;
 float  temp;
   
   if(n_of_seqs < 4) {
      fprintf(stderr,"Less than 4 sequences for statistical geometry.\n");
      return NULL;
   }

   B = (float *) space(16*sizeof(float));

   for(i=3; i<n_of_seqs; i++) {
    for(j=2; j<i; j++) {
     for(k=1; k<j; k++) {
      for(l=0; l<k; l++) {
          SingleBox(seqs[i], seqs[j], seqs[k], seqs[l]);
          SortSingleBox();
          B[0] = (float) IBox[0];      /* transfer length */
          for(i1=1;i1<=15;i1++) B[i1] += ( (float) IBox[i1] );
      }
     }
    }
   }

   temp = (float) n_of_seqs;
   temp = temp*(temp-1.)*(temp-2.)*(temp-3.)/24.;
   for(i1=1;i1<=15;i1++) B[i1] /= (temp*B[0]);

   return B;
}
   
/* ----------------------------------------------------------------------- */

PUBLIC float *statgeom4(char **ss[4], int  nn[4])
{
   int i,j,k,l;
   int i1;

   float *B;
   float  temp;
   
   B = (float *) space(16*sizeof(float));

   for(i=0; i<nn[0]; i++) {
    for(j=0; j<nn[1]; j++) {
     for(k=0; k<nn[2]; k++) {
      for(l=0; l<nn[3]; l++) {
          SingleBox(ss[0][i], ss[1][j], ss[2][k], ss[3][l]);
          B[0] = (float) IBox[0];      /* transfer length */
          for(i1=1;i1<=15;i1++) B[i1] += ( (float) IBox[i1] );
      }
     }
    }
   }

   for (temp = 1, i=0;i<4;i++) temp *= ((float) nn[i]) ;
   for(i1=1;i1<=15;i1++) B[i1] /= (temp*B[0]);

   return B;
}
/* ------------------------------------------------------------------------- */

PUBLIC void printf_stg(float *B)
{
   printf("> Statistical Geometry.\n");
   printf("> %d (sequence length)\n", (int) B[0]);
   printf("> AAAA\n");
   printf("  %7.5f\n", B[1]);
   printf("> BAAA        ABAA        AABA        AAAB\n");
   printf("  %7.5f    %7.5f     %7.5f     %7.5f\n", B[2],B[3],B[4],B[5]); 
   printf("> AABB        ABAB        ABBA\n");
   printf("  %7.5f    %7.5f     %7.5f\n", B[6],B[7],B[8]);
   printf("> AABC        ABAC        ABCA        BAAC        BACA        BCAA\n");
   printf("  %7.5f    %7.5f     %7.5f     %7.5f    %7.5f    %7.5f\n"
            ,B[9],B[10],B[11],B[12],B[13],B[14]);
   printf("> ABCD\n");
   printf("  %7.5f\n",B[15]);
}

/* ------------------------------------------------------------------------- */


PUBLIC void SimplifiedBox(float *B, char *filename)
{
   char *tmp, temp[50];
   float X,Y,Z;
   float T,P;
   float t1, t2, t3, t4;
   float p1, p2, p3, p4;
   float x1, x2, y1, y2,z1,z2;
   FILE *fp;
   float view[3] = {0.3, 1.5, 0.1};
   float axis[3] = {1.0, 0.0, 0.0};

   t1 = B[9]+B[10]+B[11]+B[15];
   t2 = B[9]+B[12]+B[13]+B[15];
   t3 = B[10]+B[12]+B[14]+B[15];
   t4 = B[11]+B[13]+B[14]+B[15];

   p1 = B[2]; p2 = B[3]; p3 = B[4]; p4 = B[5];

   x1 = B[6] + B[14];
   x2 = B[6] + B[9];
   y1 = B[7] + B[13];
   y2 = B[7] + B[10];
   z1 = B[8] + B[12];
   z2 = B[8] + B[11];

   X = (x1+x2)/2.;  Y = (y1+y2)/2.;  Z = (z1+z2)/2;
   T = (t1+t2+t3+t4)/4.;
   P = (p1+p2+p3+p4)/4.;
   
   printf("> X = %7.5f \n",X);
   printf("> Y = %7.5f \n",Y);
   printf("> Z = %7.5f \n",Z);
   printf("> T = %7.5f \n",T);
   printf("> P = %7.5f \n",P);

   tmp = get_taxon_label(-1);      /* retrieve the dataset identifier */
   
   temp[0]='\0';
   if(tmp) { strcat(temp,tmp);strcat(temp,"_"); free(tmp); }
   strcat(temp,filename);
   
   fp = fopen(temp,"w");
   if (fp!=NULL) {
      ps3d_Preambel(fp, view, axis, "N");
      PS_DrawSimplifiedBox(X,Y,Z,T,P, fp);
      fclose(fp);
   } else fprintf(stderr,"couldn't open %s -- not drawing box\n", temp);
}


      
/* ------------------------------------------------------------------------- */
      
PRIVATE void SingleBox(char *x1, char *x2, char *x3, char *x4)
{
   int   len,i1,j1,k1,i,M,m;
   int   d[4];
   char  t[4];

   len=strlen(x1);
   if(strlen(x2)!=len) nrerror("Sequences of unequal length in 'SingleBox'");
   if(strlen(x3)!=len) nrerror("Sequences of unequal length in 'SingleBox'");
   if(strlen(x4)!=len) nrerror("Sequences of unequal length in 'SingleBox'");

   IBox[0] = len;
   for(i=1; i<=15; i++) IBox[i] = 0;

   for(i1=0;i1<len;i1++) {
      t[0] = x1[i1];
      t[1] = x2[i1];
      t[2] = x3[i1];
      t[3] = x4[i1];
      for(j1=0;j1<4;j1++) {
         d[j1]=4;
         for(k1=0;k1<4;k1++) d[j1] -=(t[j1]!=t[k1]);
      }
      M=MAX4(d[0],d[1],d[2],d[3]);
      switch(M) {
       case 4 :     /* Four of a kind */
         IBox[1]++;                            /*  A A A A */
         break;
       case 3 :     /* Three of a kind */
         if(d[0] == 1)      IBox[2]++;         /*  B A A A */
         else if(d[1] == 1) IBox[3]++;         /*  A B A A */
         else if(d[2] == 1) IBox[4]++;         /*  A A B A */
         else               IBox[5]++;         /*  A A A B */
         break;
       case 2 : 
         m=MIN4(d[0],d[1],d[2],d[3]);
         if(m==2){   /* Two Pairs */
            if(t[1]==t[0])       IBox[6]++;    /*  A A B B */
            else if(t[2]==t[0])  IBox[7]++;    /*  A B A B */
            else                 IBox[8]++;    /*  A B B A */
         }
         else {      /* One Pair */  
            if(d[0]==2){  /* 0 is in the pair */
               if(t[1]==t[0])      IBox[9]++;  /*  A A B C */
               else if(t[2]==t[0]) IBox[10]++; /*  A B A C */
               else                IBox[11]++ ;/*  A B C A */
            }
            else if(d[1]==2){  /* 1 is in the pair */
               if(t[2]==t[1])      IBox[12]++; /*  B A A C */
               else                IBox[13]++; /*  B A C A */
            }
            else                   IBox[14]++; /*  B C A A */
         }
         break;
       case 1 :      /* No Pair */
         IBox[15]++;                           /*  A B C D */
         break;
       default:
         nrerror("This can't happen.");
      } 
   }
}

/* ----------------------------------------------------------------------- */

PRIVATE void SortSingleBox(void)
{
   int i;
   int M; 
   int IBB[16];
   int s[4];

   M = MAX2(MAX2( IBox[6], IBox[7] ), IBox[8] );
   
   if( M== IBox[6] ) {                   /* 12|34       */
      IBB[6] = IBox[6];
      if( IBox[9] >= IBox[14] ) {        /* 1,2 > 3,4   */
         IBB[9]  = IBox[9]; 
         IBB[14] = IBox[14];
         if( IBox[7] >= IBox[8] ) {      /*    13|24    */
            IBB[7] = IBox[7];
            IBB[8] = IBox[8];
            if( IBox[10] >= IBox[13] ) { /* 1>2>3>4     */
               IBB[10] = IBox[10];
               IBB[13] = IBox[13];
               IBB[11] = IBox[11];
               IBB[12] = IBox[12];
               s[0]=0; s[1]=1; s[2]=2; s[3]=3;            /*  1 2 3 4 */
            }
            else {                       /* 2>1>4>3     */
               IBB[10] = IBox[13];
               IBB[13] = IBox[10];
               IBB[11] = IBox[12];
               IBB[12] = IBox[11];
               s[0]=1; s[1]=0; s[2]=3; s[3]=2;            /*  2 1 4 3 */
            }
         }
         else {                          /*    14|23    */
            IBB[7] = IBox[8];
            IBB[8] = IBox[7];
            if( IBox[11] >= IBox[12]) {   /* 1>2>4>3     */
               IBB[10] = IBox[11];
               IBB[13] = IBox[12];
               IBB[11] = IBox[10];
               IBB[12] = IBox[13];
               s[0]=0; s[1]=1; s[2]=3; s[3]=2;            /*  1 2 4 3 */
            }
            else {                       /* 2>1>3>4     */
               IBB[10] = IBox[12]; 
               IBB[13] = IBox[11];
               IBB[11] = IBox[13];
               IBB[12] = IBox[10];
               s[0]=1; s[1]=0; s[2]=2; s[3]=3;            /*  2 1 3 4 */
            } 
         }
      }
      else {                             /* 3,4 > 1,2   */
         IBB[9] = IBox[14];
         IBB[14]= IBox[9];
         if( IBox[7] >= IBox[8] ) {      /*    31|42    */
            IBB[7] = IBox[7];
            IBB[8] = IBox[8];
            if(IBox[10] >= IBox[13]) {   /* 3>4>1>2     */
               IBB[10] = IBox[10];
               IBB[13] = IBox[13];
               IBB[11] = IBox[12];
               IBB[12] = IBox[11];
               s[0]=2; s[1]=3; s[2]=0; s[3]=1;            /*  3 4 1 2 */
            }
            else {                       /*            */
               IBB[10] = IBox[13];
               IBB[13] = IBox[10];
               IBB[11] = IBox[11];
               IBB[12] = IBox[12];
               s[0]=3; s[1]=2; s[2]=1; s[3]=0;            /*  4 3 2 1 */
            }
         }
         else {                          /*    32|14    */
            IBB[7] = IBox[8];
            IBB[8] = IBox[7];
            if( IBox[11] >= IBox[12]) {   /* 3>4>2>1     */
               IBB[10] = IBox[12];
               IBB[13] = IBox[11];
               IBB[11] = IBox[10];
               IBB[12] = IBox[13];
               s[0]=2; s[1]=3; s[2]=1; s[3]=0;            /*  3 4 2 1 */
            }
            else {                       /* 2>1>3>4     */
               IBB[10] = IBox[10];
               IBB[13] = IBox[13];
               IBB[11] = IBox[12];
               IBB[12] = IBox[11];
               s[0]=2; s[1]=3; s[2]=0; s[3]=1;            /*  3 4 1 2 */
            } 
         }
      }
   }
   else if (M == IBox[7] ) {             /* 13|24       */
      IBB[6] = IBox[7];
      if( IBox[10] >= IBox[13] ) {       /* 1,3 > 2,4   */
         IBB[9]  = IBox[10]; 
         IBB[14] = IBox[13];         
         if( IBox[6] >= IBox[8]) {       /*    12|34    */
            IBB[7] = IBox[6];
            IBB[8] = IBox[8];
            if( IBox[9] >= IBox[14] ) {  /*    1,2>3,4  */
               IBB[10] = IBox[9];
               IBB[13] = IBox[14];
               IBB[11] = IBox[11];
               IBB[12] = IBox[12];
               s[0]=0; s[1]=2; s[2]=1; s[3]=3;            /* 1 3 2 4 */
            }
            else {
               IBB[10] = IBox[14];
               IBB[13] = IBox[9];
               IBB[11] = IBox[12];
               IBB[12] = IBox[11];
               s[0]=2; s[1]=0; s[2]=3; s[3]=1;            /* 3 1 4 2 */
            }    
         }
         else {                          /*    14|23   */
            IBB[7] = IBox[8];
            IBB[8] = IBox[6];
            if( IBox[11] >= IBox[12] ){  /*    1,4 > 2,3 */
               IBB[10] = IBox[11];
               IBB[13] = IBox[12];
               IBB[11] = IBox[9];
               IBB[12] = IBox[14];
               s[0]=0; s[1]=2; s[2]=3; s[3]=1;            /* 1 3 4 2 */
            }
            else {
               IBB[10] = IBox[12];
               IBB[13] = IBox[11];
               IBB[11] = IBox[14];
               IBB[12] = IBox[9];
               s[0]=2; s[1]=0; s[2]=1; s[3]=3;            /* 3 1 2 4 */
            }
         }
      }
      else {                             /* 2,4 > 1,3   */
         IBB[9]  = IBox[13]; 
         IBB[14] = IBox[10];         
         if( IBox[6] >= IBox[8]) {       /*    21|43    */
            IBB[7] = IBox[6];
            IBB[8] = IBox[8];
            if( IBox[9] >= IBox[14] ) {  /*    2,1>4,3  */
               IBB[10] = IBox[9];
               IBB[13] = IBox[14];
               IBB[11] = IBox[12];
               IBB[12] = IBox[11];
               s[0]=1; s[1]=3; s[2]=0; s[3]=2;            /* 2 4 1 3 */
            }
            else {
               IBB[10] = IBox[14];
               IBB[13] = IBox[9];
               IBB[11] = IBox[11];
               IBB[12] = IBox[12];
               s[0]=3; s[1]=1; s[2]=2; s[3]=0;            /* 4 2 3 1 */
            }    
         }
         else {                          /*    14|23   */
            IBB[7] = IBox[8];
            IBB[8] = IBox[6];
            if( IBox[11] >= IBox[12] ){  /*    1,4 > 2,3 */
               IBB[10] = IBox[11];
               IBB[13] = IBox[12];
               IBB[11] = IBox[14];
               IBB[12] = IBox[9];
               s[0]=3; s[1]=1; s[2]=0; s[3]=2;            /* 4 2 1 3 */
            }
            else {
               IBB[10] = IBox[12];
               IBB[13] = IBox[11];
               IBB[11] = IBox[9];
               IBB[12] = IBox[14];
               s[0]=1; s[1]=3; s[2]=2; s[3]=0;            /* 2 4 3 1 */
            }
         }
      }
   }
   else {                                /* 14 | 23   */
      IBB[6] = IBox[8];
      if( IBox[11] >= IBox[12] ) {       /* 1,4 > 2,3 */
         IBB[9]  = IBox[11]; 
         IBB[14] = IBox[12];   
         if(IBox[6] >= IBox[7]) {        /*    12|34  */
            IBB[7] = IBox[6];
            IBB[8] = IBox[7];             
            if(IBox[9] >= IBox[14]) {   /*    1,2>3,4 */
               IBB[10] = IBox[9];
               IBB[13] = IBox[14];
               IBB[11] = IBox[10];
               IBB[12] = IBox[13];
               s[0]=0; s[1]=3; s[2]=1; s[3]=2;          /* 1 4 2 3 */
            }
            else {                      /*    4,3>2,1 */
               IBB[10] = IBox[14];
               IBB[13] = IBox[9];
               IBB[11] = IBox[13];
               IBB[12] = IBox[10];
               s[0]=3; s[1]=0; s[2]=2; s[3]=1;          /* 4 1 3 2 */
            }
         }
         else {                         /*    13|24  */
            IBB[7] = IBox[7];
            IBB[8] = IBox[6];       
            if( IBox[10] >= IBox[13]) { /*    1,3 > 2,4 */
               IBB[10] = IBox[10];
               IBB[13] = IBox[13];
               IBB[11] = IBox[9];
               IBB[12] = IBox[14];
               s[0]=0; s[1]=3; s[2]=2; s[3]=1;          /* 1 4 3 2 */
            }
            else {
               IBB[10] = IBox[13];
               IBB[13] = IBox[10];
               IBB[11] = IBox[14];
               IBB[12] = IBox[9];
               s[0]=3; s[1]=0; s[2]=1; s[3]=2;          /* 4 1 2 3 */
            }
         }
      }
      else {                             /* 2,3 > 1,4 */
         IBB[9]  = IBox[12]; 
         IBB[14] = IBox[11];   
         if(IBox[6] >= IBox[7]) {       /*    34|12  */
            IBB[7] = IBox[6];
            IBB[8] = IBox[7];
            if(IBox[9] >= IBox[14]) {   /*    2,1>4,3 */
               IBB[10] = IBox[9];
               IBB[13] = IBox[14];
               IBB[11] = IBox[13];
               IBB[12] = IBox[10];
               s[0]=1; s[1]=2; s[2]=0; s[3]=3;          /* 2 3 1 4 */
            }
            else {                      /*    4,3>2,1 */
               IBB[10] = IBox[14];
               IBB[13] = IBox[9];
               IBB[11] = IBox[10];
               IBB[12] = IBox[13];
               s[0]=2; s[1]=1; s[2]=3; s[3]=0;          /* 3 2 4 1 */
            }
         }
         else {                         /*    31|42  */
            IBB[7] = IBox[7];
            IBB[8] = IBox[6];
            if( IBox[10] >= IBox[13]) { /*    1,3 > 2,4 */
               IBB[10] = IBox[10];
               IBB[13] = IBox[13];
               IBB[11] = IBox[14];
               IBB[12] = IBox[9];
               s[0]=2; s[1]=1; s[2]=0; s[3]=4;          /* 3 2 1 4 */
            }
            else {
               IBB[10] = IBox[13];
               IBB[13] = IBox[10];
               IBB[11] = IBox[9];
               IBB[12] = IBox[14];
               s[0]=1; s[1]=2; s[2]=3; s[3]=0;          /* 2 3 4 1 */
            }
         }
      }      
   }
   /* HEUREKA */
   IBB[0] = IBox[0];
   IBB[1] = IBox[1];
   IBB[2] = IBox[2+s[0]];
   IBB[3] = IBox[2+s[1]];
   IBB[4] = IBox[2+s[2]];
   IBB[5] = IBox[2+s[3]];
   /* 
   IBB[6...14] see above :-) 
   */
   IBB[15] = IBox[15];
   
   for(i=0;i<=15;i++) IBox[i] = IBB[i];
}

/* ----------------------------------------------------------------------- */         
   
