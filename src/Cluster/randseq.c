
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXALPHA 100

#define PUBLIC
#define PRIVATE static

PRIVATE void usage(void);
PRIVATE char randaa(void);

extern long lrand48(void);
extern double drand48(void);
extern void srand48(long int seedval);

int main(int argc, char *argv[])
{
   int i,j,length=100, num=1, a_len, AA=0;
   char c[MAXALPHA], aa;
   unsigned long r;
   
   /* defaults */

   a_len = 4;   strcpy(c,"GCAU");
   length = 100;
   num    = 1;

   for (i=1; i<argc; i++) {
      if (argv[i][0] == '-') {
         if(strlen(argv[i])>2) usage();
         switch(argv[i][1]){
            case 'l' : sscanf(argv[++i], "%d", &length);
               break;
            case 'n' : sscanf(argv[++i], "%d", &num);
               break;
            case 'a' : sscanf(argv[++i], "%s", c);
               a_len = strlen(c);
               break;
	    case 'A' : AA = 1;
	      break;
            default:
               usage();
         }
      }
   }
   
   srand48(time(NULL));
   
   for (i=1; i<=num; i++) {
      if (AA) for (j=0; j<length; j++) {
	 aa = randaa();
	 printf("%c",aa);
      }  else
	 for (j=0; j<length; j++) {
	    r = (lrand48());
	    printf("%c",c[r % a_len]);
	 }
      printf("\n");
   }
   return 0;
}

PRIVATE char randaa(void)
{
   double MeanFreq[20] = {
      .0760, .0176, .0529, .0628, .0401, .0695, .0224, .0561, .0584, .0922,
      .0236, .0448, .0500, .0403, .0523, .0715, .0581, .0652, .0128, .0321};
   char AminoAcid[21] = "ACDEFGHIKLMNPQRSTVWY";
   double r;
   static double sum=0;
   int i;

   if (sum==0) for (i=0; i<20; i++) sum+=MeanFreq[i];
   i = 0;
   r = drand48()*sum;
   while ((r -= MeanFreq[i])>1e-10) i++;
   return AminoAcid[i];
}

PRIVATE void usage(void)
{
   fprintf(stderr,
	   "\nusage: randseq [-l length] [-n number] [-A] [-a alphabet]\n");
   exit(13);
}
