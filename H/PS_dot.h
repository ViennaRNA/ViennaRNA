/* routines from PS_dot.c */
extern int PS_rna_plot(char *string, char *structure, char *file);
/* write PostScript drawing of structure to file */
extern int gmlRNA(char *string, char *structure, char *ssfile, char option);
/* structure drawing in gml */
extern int ssv_rna_plot(char *string, char *structure, char *ssfile);
/*write coord file for SStructView */
extern int xrna_plot(char *string, char *structure, char *ssfile);
/*write .ss file for further editing in XRNA */
extern int PS_dot_plot(char *string, char *file);
/* produce a PostScript dot plot of the pair probability matix */
extern int rna_plot_type;   /* 0= simple coordinates, 1= naview */
