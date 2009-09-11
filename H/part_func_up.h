#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_UP_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_UP_H__

typedef struct pu_contrib { /* contributions to prob_unpaired in */
  double **H; /* hairpin loops */
  double **I; /* interior loops */
  double **M; /* multi loops */
  double **E; /* exterior loop */
  int length; /* length of the input sequence */
  int w;      /* longest unpaired region */
} pu_contrib;

typedef struct interact { /* contributions to prob_unpaired in */
  double *Pi; /* probabilities of interaction */
  double *Gi; /* free energies of interaction */
  double Gikjl; /* full free energy for interaction between [k,i] k<i
		   in longer seq and [j,l] j<l in shorter seq */
  double Gikjl_wo; /* Gikjl without contributions for prob_unpaired */
  int i; /* k<i in longer seq */
  int k; /* k<i in longer seq */
  int j; /*j<l in shorter seq */
  int l; /*j<l in shorter seq */
  int length; /* length of longer sequence */  
} interact;

/* prob of unpaired region of length w */
pu_contrib *pf_unstru(char *sequence, int max_w);
/* prob. of intermolecular interaction between two sequences of maximal length w*  prob unpaired from pf_unpaired */
interact *pf_interact(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int max_w, char *cstruc, int incr3, int incr5);
void free_pu_contrib(pu_contrib *p_con);
void free_interact(interact *pin);
int Up_plot(pu_contrib *p_c, pu_contrib *p_c_sh, interact *pint, char *ofile, int *u_vals, char *select_contrib, char *head);

#endif
