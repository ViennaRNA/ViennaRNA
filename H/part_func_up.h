typedef struct pu_contrib { /* contributions to prob_unpaired in */
  double **H; /* hairpin loops */
  double **I; /* interior loops */
  double **M; /* multi loops */
  double **E; /* exterior loop */
  int length; /* length of the input sequence */
} pu_contrib;

extern pu_contrib *pf_unstru(char *sequence, char *structure, int max_w);
extern void free_pf_unstru(void);
/* prob. of intermolecular interaction between two sequences *  prob
   from pf_unpaired */
extern double **pf_interact(const char *s1, const char *s2, pu_contrib *p_c, int max_w, int incr3, int incr5);
/* free_pf_two: first argument output of pf_unstru() !!!! */ 
extern  void free_pf_two(pu_contrib *p_con, double **pin);
extern  int Up_plot(pu_contrib *p_c, double **pint, int len, char *ofile, int max_w, char *select_contrib);
extern void scale_stru_pf_params(unsigned int length);



