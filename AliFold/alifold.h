extern float  alifold(char **strings, char *structure);
extern  void  free_alifold_arrays(void);
extern  void  update_alifold_params(void);
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
extern float alipf_fold(char **sequences, char *structure, pair_info **pi);
