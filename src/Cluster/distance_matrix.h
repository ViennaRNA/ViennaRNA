extern   float **read_distance_matrix(char type[]);
extern   char  **read_sequence_list(int *n_of_seqs,char *mask);
extern   float **Hamming_Distance_Matrix(char **seqs, int n_of_seqs);
extern   float **StrEdit_SimpleDistMatrix(char **seqs, int n_of_seqs);
extern   float **StrEdit_GotohDistMatrix(char **seqs, int n_of_seqs);
extern   char   *get_taxon_label(int whoami);
extern   void    free_distance_matrix(float **x);
extern   void    printf_distance_matrix(float **x);
extern   void    printf_taxa_list(void);
extern   float   StrEdit_SimpleDist(char *str1, char *str2);
extern   float   StrEdit_GotohDist(char *str1, char *str2);
extern   void    Set_StrEdit_CostMatrix(char type);
extern   void    Set_StrEdit_GapCosts(float per_digit, float per_gap);

