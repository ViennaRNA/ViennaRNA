#ifndef VIENNA_RNA_PACKAGE_PK_PLEX_H
#define VIENNA_RNA_PACKAGE_PK_PLEX_H

#include <ViennaRNA/datastructures/basic.h>

typedef int (vrna_callback_pk_plex_score)(const short *pt,
                                          int         i,
                                          int         j,
                                          int         k,
                                          int         l,
                                          void        *data);

typedef struct vrna_pk_plex_option_s  *vrna_pk_plex_opt_t;
typedef struct vrna_pk_plex_result_s  vrna_pk_plex_result_t;

struct vrna_pk_plex_result_s {
  char          *structure;
  double        energy;
  double        dGpk;
  double        dGint;
  double        dG1;
  double        dG2;
  unsigned int  start_5;
  unsigned int  end_5;
  unsigned int  start_3;
  unsigned int  end_3;
};

/**
 *  @brief 
 */
vrna_pkplex_t *
vrna_pk_plex(vrna_fold_compound_t *fc,
             const int            **accessibility,
             vrna_pk_plex_opt_t   options);


int **
vrna_pk_plex_accessibility(const char    *sequence,
                           unsigned int  unpaired,
                           double        cutoff);

vrna_pk_plex_opt_t
vrna_pk_plex_opt_defaults(void);


vrna_pk_plex_opt_t
vrna_pk_plex_opt(unsigned int delta,
                 unsigned int max_interaction_length,
                 int      pk_penalty);


vrna_pk_plex_opt_t
vrna_pk_plex_opt_fun(unsigned int                 delta,
                     unsigned int                 max_interaction_length,
                     vrna_callback_pk_plex_score  *scoring_function,
                     void                         *scoring_data);


#endif
