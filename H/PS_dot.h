/* routines from PS_dot.c */
extern int PS_rna_plot(char *string, char *structure, char *file);
/* write PostScript drawing of structure to file */
extern int PS_rna_plot_a(char *string, char *structure, char *file, char *pre, char *post);
/* write PostScript drawing of structure to file with annotation */
extern int gmlRNA(char *string, char *structure, char *ssfile, char option);
/* structure drawing in gml */
extern int ssv_rna_plot(char *string, char *structure, char *ssfile);
/*write coord file for SStructView */
extern int svg_rna_plot(char *string, char *structure, char *ssfile);
/*write RNAplot in SVG */
extern int xrna_plot(char *string, char *structure, char *ssfile);
/*write .ss file for further editing in XRNA */
extern int PS_dot_plot(char *string, char *file);
/* produce a PostScript dot plot of the pair probability matix */
extern int rna_plot_type;   /* 0= simple coordinates, 1= naview */

typedef struct cpair {
  int i,j,mfe;
  float p, hue, sat;
} cpair;
extern int PS_color_dot_plot(char *string, cpair *pi, char *filename);

typedef struct plist {
  int i;
  int j;
  float p;
}plist;
extern int PS_dot_plot_list(char *seq, char *filename, struct plist *pl,
			    struct plist *mf, char *comment);
extern int PS_dot_plot_turn(char *seq, struct plist *pl, char *filename,
			    int winSize);
