/*
  Plot RNA structures using different layout algorithms
  Last changed Time-stamp: <2003-09-10 13:55:01 ivo> 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "utils.h"
#include "PS_dot.h"

#define PRIVATE static

PRIVATE char  scale[] = ".........1.........2.........3.........4.........5.........6.........7.........8";

PRIVATE void usage(void);


/*--------------------------------------------------------------------------*/

extern int svg_rna_plot(char *string, char *structure, char *ssfile);

int main(int argc, char *argv[])
{
   
   char *string=NULL, *line;
   char *structure=NULL, *pre=NULL, *post=NULL;
   char  fname[13], ffname[20];
   int   i, r;
   float energy;
   int   istty;
   char  format[5]="ps";
   
   string=NULL;
   for (i=1; i<argc; i++) {
     if (argv[i][0]=='-') 
       switch ( argv[i][1] ) {
       case 't':
	 r = sscanf(argv[++i], "%d", &rna_plot_type);
	 if (r!=1) usage();
	 break;
       case 'o':
	 if (i==argc-1) usage();
	 strncpy(format, argv[i+1], 4);
	 break;
       case '-': /* long option */
	 if (strcmp(argv[i], "--pre")==0) pre=argv[++i];
	 else if (strcmp(argv[i], "--post")==0) post=argv[++i];
	 else usage();
	 break;
       default: usage(); 
       } 
   }
   
   istty = isatty(fileno(stdin));
   
   do {				/* main loop: continue until end of file */
     if (istty) {
       printf("\nInput sequence and structure; @ to quit\n");
       printf("%s\n", scale);
     }
     fname[0]='\0';
     if ((line = get_line(stdin))==NULL) break;
     
     /* skip comment lines and get filenames */
     while ((*line=='*')||(*line=='\0')||(*line=='>')) {
       if (*line=='>')
	 sscanf(line, ">%12s", fname);
       printf("%s\n", line);
       free(line);
       if ((line = get_line(stdin))==NULL) line = "@";
     } 
     
     if (strcmp(line, "@") == 0) break;
     
     string = (char *) space(strlen(line)+1);
     sscanf(line,"%s",string);
     free(line);
     
     if ((line = get_line(stdin))==NULL) break;
     structure = (char *) space(strlen(line)+1);
     sscanf(line,"%s (%f)", structure, &energy);
     free(line);

     if (strlen(string)!=strlen(structure)) 
       nrerror("sequence and structure have unequal length!");

     if (fname[0]!='\0') {
       strcpy(ffname, fname);
       strcat(ffname, "_ss");
     } else
       strcpy(ffname, "rna");

     switch (format[0]) {
     case 'p':
       strcat(ffname, ".ps");
       PS_rna_plot_a(string, structure, ffname, pre, post);
       break;
     case 'g':
       strcat(ffname, ".gml");
       gmlRNA(string, structure, ffname, 'x');
       break;
     case 'x':
       strcat(ffname, ".ss");
       xrna_plot(string, structure, ffname);
       break;
     case 's':
       strcat(ffname, ".svg");
       svg_rna_plot(string, structure, ffname);
       break;
     default:
       usage();
     }
     fflush(stdout);
     free(string);
     free(structure); 
   } while (1);
   return 0;
}

PRIVATE void usage(void)
{
  nrerror("usage: RNAplot [-t 0|1] [-o ps|gml|xrna|svg] [--pre <string>] [--post <string>]");
}
