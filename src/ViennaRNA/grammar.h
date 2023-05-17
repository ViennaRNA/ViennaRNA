#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_H

/**
 *  @file     grammar.h
 *  @ingroup  grammar
 *  @brief    Implementations for the RNA folding grammar
 */

/**
 * @addtogroup grammar
 * @{
 *    @brief  The RNA folding grammar as implemented in RNAlib
 */

#include <ViennaRNA/fold_compound.h>

typedef int (*vrna_grammar_rule_f)(vrna_fold_compound_t  *fc,
                                    int                   i,
                                    int                   j,
                                    void                  *data);


typedef void (*vrna_grammar_rule_f_aux)(vrna_fold_compound_t  *fc,
                                    int                   i,
                                    int                   j,
                                    void                  *data);


typedef FLT_OR_DBL (*vrna_grammar_rule_f_exp)(vrna_fold_compound_t *fc,
                                               int                  i,
                                               int                  j,
                                               void                 *data);


typedef void (*vrna_grammar_rule_f_aux_exp)(vrna_fold_compound_t *fc,
                                               int                  i,
                                               int                  j,
                                               void                 *data);


typedef void (*vrna_grammar_cond_f)(vrna_fold_compound_t *fc,
                                     unsigned char        stage,
                                     void                 *data);


/**
 *  @brief  Free auxiliary data
 *  @param  data  The auxiliary data to be free'd
 */
typedef void (*vrna_grammar_data_free_f)(void *data);


typedef struct vrna_gr_aux_s vrna_gr_aux_t;


struct vrna_gr_aux_s {
  vrna_grammar_cond_f       cb_proc; /**< @brief A callback for pre- and post-processing of auxiliary grammar rules */

  vrna_grammar_rule_f       cb_aux_f;
  vrna_grammar_rule_f       cb_aux_c;
  vrna_grammar_rule_f       cb_aux_m;
  vrna_grammar_rule_f       cb_aux_m1;
  vrna_grammar_rule_f_aux       cb_aux;

  vrna_grammar_rule_f_exp   cb_aux_exp_f;
  vrna_grammar_rule_f_exp   cb_aux_exp_c;
  vrna_grammar_rule_f_exp     cb_aux_exp_m;
  vrna_grammar_rule_f_exp     cb_aux_exp_m1;
  vrna_grammar_rule_f_aux_exp   cb_aux_exp;

  void                        *data;
  vrna_grammar_data_free_f  free_data;
};


int
vrna_gr_set_aux_f(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f cb);


int
vrna_gr_set_aux_exp_f(vrna_fold_compound_t      *fc,
                      vrna_grammar_rule_f_exp cb);


int
vrna_gr_set_aux_c(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f cb);


int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t      *fc,
                      vrna_grammar_rule_f_exp cb);


int
vrna_gr_set_aux_m(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f cb);


int
vrna_gr_set_aux_exp_m(vrna_fold_compound_t      *fc,
                      vrna_grammar_rule_f_exp cb);


int
vrna_gr_set_aux_m1(vrna_fold_compound_t   *fc,
                   vrna_grammar_rule_f  cb);


int
vrna_gr_set_aux_exp_m1(vrna_fold_compound_t       *fc,
                       vrna_grammar_rule_f_exp  cb);


int
vrna_gr_set_aux(vrna_fold_compound_t  *fc,
                vrna_grammar_rule_f_aux cb);


int
vrna_gr_set_aux_exp(vrna_fold_compound_t      *fc,
                    vrna_grammar_rule_f_aux_exp cb);


int
vrna_gr_set_data(vrna_fold_compound_t       *fc,
                 void                       *data,
                 vrna_grammar_data_free_f free_data);


int
vrna_gr_set_cond(vrna_fold_compound_t   *fc,
                 vrna_grammar_cond_f  cb);


int
vrna_gr_reset(vrna_fold_compound_t *fc);


/**
 * @}
 */

/**
 *  @addtogroup domains
 *  @brief This module covers simple and straight-forward extensions to the RNA folding grammar.
 */

#endif
