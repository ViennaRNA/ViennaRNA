#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "utils.h"
#include "StrEdit_CostMatrix.h"

#define  PUBLIC
#define  PRIVATE         static
#define  MAXSEQS         1000
#define  MIN(A, B)       ((A) < (B) ? (A) : (B))   
#define  MIN3(A,B,C)     (MIN((MIN((A),(B))),(C)))
#define  MAX(A, B)       ((A) > (B) ? (A) : (B))

PUBLIC   float **read_distance_matrix(char type[]);
PUBLIC   char  **read_sequence_list(int *n_of_seqs, char *mask);
PUBLIC   float **Hamming_Distance_Matrix(char **seqs, int n_of_seqs);
PUBLIC   float **StrEdit_SimpleDistMatrix(char **seqs, int n_of_seqs);
PUBLIC   float **StrEdit_GotohDistMatrix(char **seqs, int n_of_seqs);
PUBLIC   void    free_distance_matrix(float **x);
PUBLIC   void    printf_distance_matrix(float **x);
PUBLIC   void    printf_taxa_list(void);
PUBLIC   char   *get_taxon_label(int whoami);
PUBLIC   float   StrEdit_SimpleDist(char *str1, char *str2);
PUBLIC   float   StrEdit_GotohDist(char *str1, char *str2);
PUBLIC   void    Set_StrEdit_CostMatrix(char type);
PUBLIC   void    Set_StrEdit_GapCosts(float per_digit, float per_gap);

/* NOTE:   x[0][0] = (float)size_of_matrix;    */

PRIVATE  void    read_taxa_list(void);
PRIVATE  int     string_consists_of(char line[],char *mask);
PRIVATE  float   StrEditCost( int i, int j, char *T1, char *T2);
PRIVATE  int     decode(char id);

PRIVATE  char    Taxa_List[MAXSEQS][50];
PRIVATE  int     Taxa_Numbers[MAXSEQS];
PRIVATE  int     N_of_named_taxa=0;    
PRIVATE  char   *file_name;
PRIVATE  char    N_of_infiles=0;
PRIVATE  float **StrEdit_CostMatrix;
PRIVATE  char   *StrEdit_ValidAlphabet;
PRIVATE  float   StrEdit_GapCost   = 1.;
PRIVATE  float   StrEdit_GotohAlpha = 1.;
PRIVATE  float   StrEdit_GotohBeta  = 1.;



PUBLIC float **read_distance_matrix(char type[])
{
   char   *line;
   float **D;
   float   tmp;
   int     i,j,size;
   
   while(1) {
     type[0]= '\0';
     size   =    0;
     D      = NULL;
     if ((line = get_line(stdin))==NULL) return NULL;
     if (*line =='@') return NULL;
     if (*line =='*') {
       N_of_infiles++;
       if(file_name) free(file_name);
       
       if(strlen(line)>1) {
	 file_name = (char *) space(sizeof(char)*strlen(line));
	 sscanf(line,"*%s",file_name);
       } else {
	 file_name = (char *) space(10);
	 sprintf(file_name,"%d",N_of_infiles);
       }
       read_taxa_list();
     } 
     else if (*line=='>') {
       int r;
       size = 0; 
       r = sscanf(line,"> %1s%*[ ] %d", type, &size);
       fprintf(stderr, "%d ", r);
       if (r==EOF) return NULL;
       if((r==2)&&(size>1)) {
	 D=(float **)space((size+1)*sizeof(float *));
	 for(i=0; i<=size; i++)
	   D[i] = (float *)space((size+1)*sizeof(float));
	 D[0][0] = (float)size;
	 D[1][1] = 0.0;
	 for(i=2; i<= size; i++) {
	   D[i][i] = 0.0;
	   for(j=1; j<i; j++) {
	     if (scanf("%f", &tmp)!=1) {
	       for(i=0;i<=size;i++) free(D[i]);
	       free(D);
	       return NULL;
	     }
	     D[i][j] = tmp;
	     D[j][i] = tmp;
	   }
	 }
	 return D;
       }
       else printf("%s\n",line);
     }
     else printf(" %s\n", line);
     free(line);
   }
}

/* ------------------------------------------------------------------------- */

PUBLIC char **read_sequence_list(int *n_of_seqs, char *mask)
{
   int     i;
   char   *line;
   char   *tt[MAXSEQS];
   char  **sl;
   int     len;
   
   (*n_of_seqs) = 0;
   while(1) {
      if ((line = get_line(stdin))==NULL) break;
      
      if(line[0]=='\0') break;
      if(line[0]=='%')  break;
      if(line[0]=='#')  break;
      if(line[0]=='@')  break;
      if(line[0]=='*') {
         N_of_infiles++;
         if(file_name) free(file_name);
         if(strlen(line)>1) {
	    file_name = (char *) space(sizeof(char)*strlen(line));
            sscanf(line,"*%s",file_name);
	 } else {
	    file_name = (char *) space(10);
            sprintf(file_name,"%d",N_of_infiles);
	 }
         read_taxa_list();
	 free(line);
	 continue;
      }

      len = strlen(line);
      if(string_consists_of(line,mask)){
	 if(mask[0]=='%') {
	    for(i=0;i<len;i++){
	       if(isalpha(line[i])) line[i]=toupper(line[i]);
	    }
	    if(mask[0]=='!') {
	       for(i=0;i<=len;i++){
		  switch(toupper(line[i])){
		   case 'G' :  case 'A' : case 'X' :  
		      line[i] = 'R'; break;
		    case 'U' : case 'C' : case 'K' : case 'T' :
		       line[i] = 'Y'; break;
		    default:
		      line[i] = '*';
		   }
	       }
	    }
	    tt[*n_of_seqs] = (char *)space((len+1)*sizeof(char));
	    sscanf(line,"%s",tt[*n_of_seqs]);
	    (*n_of_seqs)++;
	 }
      }
      free(line);
   }
   if(*n_of_seqs == 0) return NULL;
   else {
     sl = (char **) space((*n_of_seqs)*sizeof(char *));
     for(i=0;i<*n_of_seqs; i++) sl[i] = tt[i]; 
     return sl;
   }
}

/* -------------------------------------------------------------------------- */

PRIVATE int string_consists_of(char *line, char *mask)
{
   int i,j,len;
   if((len=strlen(line))==0) return 0; /* empty string ! */
   for(j=0;j<len;j++){
      for(i=1;i<strlen(mask);i++){     /* mask[0] marks toupper/not_toupper ! */
         if(mask[i]==line[j]) break;
      }
      if(i==(strlen(mask))) break;      /* letter j no in mask */
   }
   if(j!=len) return 0;       /* loop left before last char -> false */
   return 1;                           /* loop left after last char -> true */
}
/* -------------------------------------------------------------------------- */

PUBLIC void free_distance_matrix(float **x)
{
   int i,n;
   n=(int) x[0][0];
   for(i=0;i<=n;i++) free(x[i]);
   free(x);
   x=NULL;
}

/* -------------------------------------------------------------------------- */

PUBLIC void printf_distance_matrix(float **x)
{
   int i,j,n;
   n=(int) x[0][0];
   printf("> X  %d\n",n);
   if(n>1){
      for(i=2;i<=n;i++) {
         for(j=1;j<i;j++) printf("%g ",x[i][j]);
         printf("\n");
      }
   }
}

/* -------------------------------------------------------------------------- */

PUBLIC float **Hamming_Distance_Matrix(char **seqs, int n_of_seqs)
{
   int i,j,k;
   float **D;
   D = (float **) space((n_of_seqs+1)*sizeof(float *));
   for(i=0;i<=n_of_seqs;i++)
      D[i] = (float *) space((n_of_seqs+1)*sizeof(float));
   D[0][0] = (float) n_of_seqs;
   
   for(i=1; i<n_of_seqs; i++) {
      D[i][i] = 0.;
      for(j=0;j<i;j++){
         if(strlen(seqs[i])!=strlen(seqs[j])) 
            nrerror("Unequal Seqence Length for Hamming Distance.");
         D[i+1][j+1] = 0.0;
         for(k=0;k<strlen(seqs[i]);k++)
            D[i+1][j+1] += StrEditCost(k+1,k+1,seqs[i],seqs[j]);
	 /* was :  (float)(seqs[i][k]!=seqs[j][k]); */
         D[j+1][i+1] = D[i+1][j+1];
      }
      D[n_of_seqs][n_of_seqs] = 0.;
   }
   return D;
}

/* -------------------------------------------------------------------------- */

PUBLIC float **StrEdit_SimpleDistMatrix(char **seqs, int n_of_seqs)
{
   int i,j;
   float **D;
   D = (float **) space((n_of_seqs+1)*sizeof(float *));
   for(i=0;i<=n_of_seqs;i++)
      D[i] = (float *) space((n_of_seqs+1)*sizeof(float));
   D[0][0] = (float) n_of_seqs;
   
   for(i=1; i<n_of_seqs; i++) {
      D[i][i] = 0.;
      for(j=0;j<i;j++){
         D[i+1][j+1] = StrEdit_SimpleDist(seqs[i],seqs[j]);
         D[j+1][i+1] = D[i+1][j+1];
      }
      D[n_of_seqs][n_of_seqs] = 0.;
   }
   return D;
}

/* -------------------------------------------------------------------------- */

PUBLIC float **StrEdit_GotohDistMatrix(char **seqs, int n_of_seqs)
{
   int i,j;
   float **D;
   D = (float **) space((n_of_seqs+1)*sizeof(float *));
   for(i=0;i<=n_of_seqs;i++)
      D[i] = (float *) space((n_of_seqs+1)*sizeof(float));
   D[0][0] = (float) n_of_seqs;
   
   for(i=1; i<n_of_seqs; i++) {
      D[i][i] = 0.;
      for(j=0;j<i;j++){
         D[i+1][j+1] = StrEdit_GotohDist(seqs[i],seqs[j]);
         D[j+1][i+1] = D[i+1][j+1];
      }
      D[n_of_seqs][n_of_seqs] = 0.;
   }
   return D;
}

/* -------------------------------------------------------------------------- */

PRIVATE void read_taxa_list(void)
{
   char *line;
   int i,add_it;

   add_it=0;
   if ((line = get_line(stdin))==NULL) return;
   if(line[0]=='#') {
      i=0;
      sscanf(line,"#%d", &i);
      if(i<=1) N_of_named_taxa = 0;
      add_it = i*100000;
      free(line);
      if ((line = get_line(stdin))==NULL) return;
   } else N_of_named_taxa = 0;
   
   do {
      if(line[0]=='\0') break;
      if(line[0]=='%')  break;
      if(line[0]=='#')  break;
      if(line[0]=='@')  break;
      if(line[0]=='*')  break;
      *Taxa_List[N_of_named_taxa]='\0';
      sscanf(line,"%d :%49s", &i, Taxa_List[N_of_named_taxa]);
      if(*Taxa_List[N_of_named_taxa]) { 
         Taxa_Numbers[N_of_named_taxa]=i+add_it;
	 N_of_named_taxa++;
      }
      free(line);
   } while ((line = get_line(stdin))!=NULL);
   if (line!=NULL) free(line);
   return;
}

/* -------------------------------------------------------------------------- */

PUBLIC void printf_taxa_list(void)
{
   int i;
   if(N_of_named_taxa>0){
      printf("* List of Taxa: %s\n", file_name);
      for(i=0;i<N_of_named_taxa;i++)
         printf("%3d : %s\n",Taxa_Numbers[i],Taxa_List[i]);   
      printf("* End of Taxa List\n");
   }
}

/* -------------------------------------------------------------------------- */

PUBLIC char *get_taxon_label(int whoami)
{
   char *label;
   char tmp[20];
   int i;

   if(whoami<0) {    /* negative arguments return the identifier of the data set */
      if(!file_name) return NULL;
      label = (char *) space(sizeof(char)*(strlen(file_name)+1));
      strcpy(label,file_name);
      return label;
   }
   for(i=0;i<N_of_named_taxa;i++) {
      if(whoami==Taxa_Numbers[i]) {
         label = (char *) space(sizeof(char)*(strlen(Taxa_List[i])+1));
         strcpy(label,Taxa_List[i]);
         return label;
      }    
   }
   sprintf(tmp,"%d",whoami);
   
   label = (char *) space(sizeof(char)*(strlen(tmp)+1));
   strcpy(label,tmp);
   return label;
}

/* -------------------------------------------------------------------------- */

PUBLIC float StrEdit_SimpleDist(char *str1, char *str2 )

{
   float  **distance;

   int           i, j, length1,length2;
   float         minus, plus, change, temp;
    
   length1 = strlen(str1);
   length2 = strlen(str2);

   distance = (float **)  space((length1 +1)*sizeof(float *));
   for(i=0; i<= length1; i++)
      distance[i] = (float *) space( (length2+1)*sizeof(float));

   for(i = 1; i <= length1; i++) 
      distance[i][0] = distance[i-1][0]+StrEditCost(i,0,str1,str2);
   for(j = 1; j <= length2; j++) 
      distance[0][j] = distance[0][j-1]+StrEditCost(0,j,str1,str2);
    
   for (i = 1; i <= length1; i++) {
      for (j = 1; j <= length2 ; j++) {
         minus  = distance[i-1][j]  + StrEditCost(i,0,str1,str2);
         plus   = distance[i][j-1]  + StrEditCost(0,j,str1,str2);
         change = distance[i-1][j-1]+ StrEditCost(i,j,str1,str2);
            
         distance[i][j] = MIN3(minus, plus, change);  
      } 
   }
   temp = distance[length1][length2];
   for(i=0;i<=length1;i++) free(distance[i]);
   free(distance);

   return temp;
}

/* -------------------------------------------------------------------------- */

PUBLIC float StrEdit_GotohDist(char *str1, char *str2 )
{
   float  **D;
   float  **E;
   float  **F;
   int      i, j, length1,length2;
   float    temp;
    
   length1 = strlen(str1);
   length2 = strlen(str2);
  
   D = space((length1+1)*sizeof(float *));
   for(i=0;i<=length1;i++) D[i] = space((length2+1)*sizeof(float));
   E = space((length1+1)*sizeof(float *));
   for(i=0;i<=length1;i++) E[i] = space((length2+1)*sizeof(float));
   F = space((length1+1)*sizeof(float *));
   for(i=0;i<=length1;i++) F[i] = space((length2+1)*sizeof(float));

   D[0][0] = 0.; E[0][0] = 0.; F[0][0] = 0.;
   
   for(i=1;i<=length1;i++) {
      D[i][0] = StrEdit_GotohAlpha + StrEdit_GotohBeta*((float)(i-1));
      E[i][0] = StrEdit_GotohAlpha + StrEdit_GotohBeta*((float)(i-1));
      F[i][0] = 0.;
   }
   for(j=1;j<=length2;j++) {
      D[0][j] = StrEdit_GotohAlpha + StrEdit_GotohBeta*((float)(j-1));
      E[0][j] = 0.;
      F[0][j] = StrEdit_GotohAlpha + StrEdit_GotohBeta*((float)(j-1));
   }
   for(i=1;i<=length1;i++) {
      for(j=1;j<=length2;j++) {
         E[i][j] = MIN(  (D[i][j-1]+StrEdit_GotohAlpha), 
                         (E[i][j-1]+StrEdit_GotohBeta)  );
         F[i][j] = MIN(  (D[i-1][j]+StrEdit_GotohAlpha),
                         (F[i-1][j]+StrEdit_GotohBeta)  );
         D[i][j] = MIN3(  E[i][j], F[i][j], 
                         (D[i-1][j-1]+StrEditCost(i,j,str1,str2)) );
      }
   }
   temp = D[length1][length2];
   for(i=0;i<=length1;i++) {
      free(D[i]); free(E[i]); free(F[i]);
   }
   free(D); free(E); free(F);
   
   return temp;
}

/* -------------------------------------------------------------------------- */

PRIVATE float StrEditCost(int i, int j, char *T1, char *T2)
{
   /* positions i,j from [1..length]; i,j=0 implies Gap */
   int i1,j1;
   if((i==0)&&(j==0)) nrerror("Edit Cost: Aligned gap characters !!!");
   if(i>0) i1 = decode(T1[i-1]); else i1 = 0;
   if(j>0) j1 = decode(T2[j-1]); else j1 = 0;
   if(StrEdit_CostMatrix==NULL) {
      if(i&&j) return (float)(i1!=j1);
      else     return (float) StrEdit_GapCost;
   }
   else return (float) StrEdit_CostMatrix[i1][j1];
}

/* -------------------------------------------------------------------------- */

PRIVATE int decode(char id)
{
   int   n,alen;
   
   if (!StrEdit_ValidAlphabet) return (int)id;
   alen = strlen(StrEdit_ValidAlphabet);
   if(!alen) return (int)id;
   
   for(n=0;n<alen;n++) {
      if(id==StrEdit_ValidAlphabet[n]) return n;
   }
   fprintf(stderr,"Warning: Invalid character in DECODE -> set to ~gap~\n");
   return 0;
}   

/* -------------------------------------------------------------------------- */

PUBLIC void  Set_StrEdit_CostMatrix(char type)
{
   int     i,j;
   if(StrEdit_ValidAlphabet) { 
      free(StrEdit_ValidAlphabet); 
      StrEdit_ValidAlphabet = NULL;
   }
   if(StrEdit_CostMatrix) { 
      free(StrEdit_CostMatrix); 
      StrEdit_CostMatrix = NULL;
   }
   switch(type){
      case 'D' :
        StrEdit_ValidAlphabet = (char*) space((20+2)*sizeof(char));
        strcpy(StrEdit_ValidAlphabet,StrEdit_DayhoffA);
        StrEdit_CostMatrix    = (float**) space((20+1)*sizeof(float*));
        for(i=0;i<=20;i++) 
           StrEdit_CostMatrix[i] = (float*)space((20+1)*sizeof(float));
        for(i=1;i<=20;i++) { 
           for(j=1;j<=20;j++) {
              StrEdit_CostMatrix[i][j] = 
                 MAX(StrEdit_DayhoffM[i][i],StrEdit_DayhoffM[j][j])-
                                                StrEdit_DayhoffM[i][j];
           }
           StrEdit_CostMatrix[i][0] = StrEdit_DayhoffM[i][i];
           StrEdit_CostMatrix[0][i] = StrEdit_DayhoffM[i][i];
        }
        StrEdit_CostMatrix[0][0] = StrEdit_DayhoffM[0][0];
        break;
      case 'A' :
	StrEdit_ValidAlphabet = (char*) space((20+2)*sizeof(char));
        strcpy(StrEdit_ValidAlphabet,StrEdit_GLHA);
	StrEdit_CostMatrix    = (float**) space((20+1)*sizeof(float*));
        for(i=0;i<=20;i++) 
           StrEdit_CostMatrix[i] = (float*)space((20+1)*sizeof(float));
        for(i=1;i<=20;i++) { 
           for(j=1;j<=20;j++) 
              StrEdit_CostMatrix[i][j] = StrEdit_GLHM[i][j];
	   StrEdit_CostMatrix[i][0] = StrEdit_CostMatrix[0][i] =
	      StrEdit_GapCost;
        }
        StrEdit_CostMatrix[0][0] = StrEdit_GLHM[0][0];
        break;
      case 'B' :
        StrEdit_ValidAlphabet = space((4+2)*sizeof(char));
        strcpy(StrEdit_ValidAlphabet,StrEdit_BinCodeA);
        StrEdit_CostMatrix    = (float**) space((4+1)*sizeof(float*));
        for(i=0;i<=4;i++) 
           StrEdit_CostMatrix[i] = (float*)space((4+1)*sizeof(float));
        for(i=0;i<=4;i++) 
           for(j=0;j<=4;j++) 
              StrEdit_CostMatrix[i][j] = StrEdit_BinCodeM[i][j];
        break;
      case 'H' :
        StrEdit_ValidAlphabet = space((4+2)*sizeof(char));
        strcpy(StrEdit_ValidAlphabet,StrEdit_HogewegA);
        StrEdit_CostMatrix    = (float**) space((4+1)*sizeof(float*));
        for(i=0;i<=4;i++) 
           StrEdit_CostMatrix[i] = (float*)space((4+1)*sizeof(float));
        for(i=0;i<=4;i++) 
           for(j=0;j<=4;j++) 
              StrEdit_CostMatrix[i][j] = StrEdit_HogewegM[i][j];
        StrEdit_GotohAlpha    = 3.;
        StrEdit_GotohBeta     = 0.;
        break;
   default: 
        if(!StrEdit_GapCost) StrEdit_GapCost = 1.;  /* This is the simple distance */
   }
}

/* -------------------------------------------------------------------------- */

PUBLIC   void    Set_StrEdit_GapCosts(float per_digit, float per_gap)
{
   if(per_gap==0.) per_gap = per_digit;
   if(per_digit<0) nrerror("Gap Costs invalid.");
   if(per_digit>per_gap) nrerror("Gap Costs invalid.");

   StrEdit_GapCost     = per_digit;
   StrEdit_GotohAlpha  = per_digit;   /* Gotoh gap function g(k) = a + b(k-1) */
   StrEdit_GotohBeta   = per_gap;

}
     
/* -------------------------------------------------------------------------- */
