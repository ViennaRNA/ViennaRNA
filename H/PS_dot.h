#ifndef __VIENNA_RNA_PACKAGE_PS_DOT_H__
#define __VIENNA_RNA_PACKAGE_PS_DOT_H__

/* routines from PS_dot.c */
int PS_rna_plot(char *string, char *structure, char *file);
/* write PostScript drawing of structure to file */
int PS_rna_plot_a(char *string, char *structure, char *file, char *pre, char *post);
/* write PostScript drawing of structure to file with annotation */
int gmlRNA(char *string, char *structure, char *ssfile, char option);
/* structure drawing in gml */
int ssv_rna_plot(char *string, char *structure, char *ssfile);
/*write coord file for SStructView */
int svg_rna_plot(char *string, char *structure, char *ssfile);
/*write RNAplot in SVG */
int xrna_plot(char *string, char *structure, char *ssfile);
/*write .ss file for further editing in XRNA */
int PS_dot_plot(char *string, char *file);
/* produce a PostScript dot plot of the pair probability matix */
int rna_plot_type;   /* 0= simple coordinates, 1= naview */

typedef struct cpair {
  int i,j,mfe;
  float p, hue, sat;
} cpair;
int PS_color_dot_plot(char *string, cpair *pi, char *filename);
int PS_color_dot_plot_turn(char *seq, cpair *pi, char *filename, int winSize);

typedef struct plist {
  int i;
  int j;
  float p;
}plist;
int PS_dot_plot_list(char *seq, char *filename, struct plist *pl,
         struct plist *mf, char *comment);
int PS_dot_plot_turn(char *seq, struct plist *pl, char *filename,
         int winSize);
int PS_color_aln(const char *structure, const char *filename, 
			const char *seqs[], const char *names[]);

#endif
