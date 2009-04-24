extern float  alifold(char **strings, char *structure);
extern  void  free_alifold_arrays(void);
extern double cv_fact /* =1 */;
extern double nc_fact /* =1 */;
typedef struct {
   short i;        /* i,j in [0, n-1] */
   short j;
   float p;      /* probability */
   float ent;    /* pseudo entropy for p(i,j) = S_i + S_j - p_ij*ln(p_ij) */
   short bp[8];  /* frequencies of pair_types */
   char comp;    /* 1 iff pair is in mfe structure */
}  pair_info;

/* extern float aliLfold(char **strings, char *structure, int maxdist); */
extern float alipf_fold(char **sequences, char *structure, struct plist **pl);
/* extern float alipfW_fold(char **sequences, char *structure, struct plist **pl,int winSize, float cutoff, int pairsize); */
/* extern struct plist *get_plistW(struct plist *pl, int length, double cutoff, int start, FLT_OR_DBL **Tpr, int winSize); */
extern char *centroid_ali(int length, double *dist,struct plist *pl);
extern float **readribosum(char *name);
extern char *alipbacktrack(double *prob) ;
extern void  free_alipf_arrays(void);
extern float  energy_of_alistruct(char **sequences, const char *structure, int n_seq, float *CVenergy);
extern float circalifold(const char **strings, char *structure);
extern float alipf_circ_fold(char **sequences, char *structure, struct plist **pl);
