/* Last Changed Time-stamp: <2001-09-07 10:17:47 ivo> */
/* This program converts the bracket notation for RNA secondary structures
   produced by RNAfold to .ct files used by Michael Zukers Program.
   To compile enter:
                    cc -o b2ct b2ct.c
   And use as
                    b2ct < structure_file > ct_file.ct
   or
                    RNAfold < sequence_file | b2ct > ct_file.ct
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXLENGTH 30000

void write_ct_file(char *fname, char *sequence, char *structure, char *name,
		   float energy);
short *make_pair_table(const char *structure);
void *space(unsigned size);
void nrerror(char *message);

int main(void)
{
   char line[MAXLENGTH+1];
   char *string=NULL, *structure=NULL, *name=NULL;
   float energy;
   int n=0;
   
   while (fgets(line, MAXLENGTH, stdin)!=NULL) {
       if (strcmp(line,"@")==0) break;
      
       switch (line[0]) {
       case '>': name = (char *) space(strlen(line));
	 sscanf(line,"> %s", name);
	 break;
       case '.':
       case '(':
       case ')': structure = (char *) space(strlen(line));
	 if (sscanf(line,"%s (%f)", structure, &energy)!=2) {
	   free(structure); structure=NULL; break;
	 }
	 n++;
	 break;
       default: string = (char *) space(strlen(line)+1);
	 sscanf(line, "%s", string);
       }
       if (structure!=NULL) {
	 if (name==NULL) {
	   name = (char *) space(10);
	   sprintf(name,"%d",n);
	 }
	 write_ct_file("-", string, structure, name, energy);
      	 free(string); 
	 free(structure); 
	 free(name);
	 string = structure = name = NULL;
       }
   }
   return 0;
}

void write_ct_file(char *fname, char *sequence, char *structure, char *name,
		   float energy)
{
   int i, length;
   short *table;
   FILE *ct;

   length = strlen(structure);
   if ((table = make_pair_table(structure))==NULL) {
     return;
   }
   if (length!=strlen(sequence))
     nrerror("sequence and structure have unequal length");
   
   if (strcmp(fname,"-")==0) 
     ct = stdout;
   else {
     ct = fopen(fname, "a");
     if (ct==NULL) nrerror("can't open .ct file");
   }

   fprintf(ct, "%5d ENERGY = %7.1f    %s\n", length, energy, name);
   for (i=1; i<=length; i++) 
     fprintf(ct, "%5d %c   %5d %4d %4d %4d\n",
	     i, sequence[i-1], i-1, (i+1)%(length+1), table[i], i);
   if (strcmp(fname,"-"))
     fclose(ct);
   else fflush(ct);
}

short *make_pair_table(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   int i,j,hx;
   int length;
   short *stack;
   short *table;
   
   length = strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;
   
   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(': 
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "unbalanced brackets in %s\n", structure);
	    free(stack); free(table); return NULL;
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   free(stack);
   if (hx!=0) {
      fprintf(stderr, "unbalanced brackets %s\n", structure);
      free(table);
      return NULL;
   }
   return(table);
}

void *space(unsigned size)
{
    void *pointer;
    
    if ( (pointer = (void *) calloc(1, size)) == NULL) {
       fprintf(stderr,"SPACE: requested size: %d\n", size);
       nrerror("SPACE allocation failure -> no memory");
    }
    return  pointer;
}


void nrerror(char *message)       /* output message upon error */
{
    fprintf(stderr, "\n%s\n", message);
    exit(0);
}

