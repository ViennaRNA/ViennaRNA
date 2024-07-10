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
 */

#include <ViennaRNA/fold_compound.h>

typedef int (*vrna_grammar_rule_f)(vrna_fold_compound_t *fc,
                                   int                  i,
                                   int                  j,
                                   void                 *data);


typedef void (*vrna_grammar_rule_f_aux)(vrna_fold_compound_t  *fc,
                                        int                   i,
                                        int                   j,
                                        void                  *data);


typedef FLT_OR_DBL (*vrna_grammar_rule_f_exp)(vrna_fold_compound_t  *fc,
                                              int                   i,
                                              int                   j,
                                              void                  *data);


typedef void (*vrna_grammar_rule_f_aux_exp)(vrna_fold_compound_t  *fc,
                                            int                   i,
                                            int                   j,
                                            void                  *data);


typedef struct vrna_gr_aux_s vrna_gr_aux_t;


struct vrna_gr_aux_s {
  vrna_grammar_rule_f         cb_aux_f;
  vrna_grammar_rule_f         cb_aux_c;
  vrna_grammar_rule_f         cb_aux_m;
  vrna_grammar_rule_f         cb_aux_m1;
  vrna_grammar_rule_f_aux     cb_aux;

  vrna_grammar_rule_f_exp     cb_aux_exp_f;
  vrna_grammar_rule_f_exp     cb_aux_exp_c;
  vrna_grammar_rule_f_exp     cb_aux_exp_m;
  vrna_grammar_rule_f_exp     cb_aux_exp_m1;
  vrna_grammar_rule_f_aux_exp cb_aux_exp;

  vrna_recursion_status_f     cb_proc; /**< @brief A callback for pre- and post-processing of auxiliary grammar rules */
  void                        *data;
  vrna_auxdata_prepare_f      prepare_data;
  vrna_auxdata_free_f         free_data;
};


int
vrna_gr_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options);


int
vrna_gr_set_aux_f(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f   cb);


int
vrna_gr_set_aux_exp_f(vrna_fold_compound_t    *fc,
                      vrna_grammar_rule_f_exp cb);


int
vrna_gr_set_aux_c(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f   cb);


int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t    *fc,
                      vrna_grammar_rule_f_exp cb);


int
vrna_gr_set_aux_m(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f   cb);


int
vrna_gr_set_aux_exp_m(vrna_fold_compound_t    *fc,
                      vrna_grammar_rule_f_exp cb);


int
vrna_gr_set_aux_m1(vrna_fold_compound_t *fc,
                   vrna_grammar_rule_f  cb);


int
vrna_gr_set_aux_exp_m1(vrna_fold_compound_t     *fc,
                       vrna_grammar_rule_f_exp  cb);


int
vrna_gr_set_aux(vrna_fold_compound_t    *fc,
                vrna_grammar_rule_f_aux cb);


int
vrna_gr_set_aux_exp(vrna_fold_compound_t        *fc,
                    vrna_grammar_rule_f_aux_exp cb);


int
vrna_gr_set_data(vrna_fold_compound_t   *fc,
                 void                   *data,
                 vrna_auxdata_prepare_f prepare_cb,
                 vrna_auxdata_free_f    free_cb);


int
vrna_gr_set_cond(vrna_fold_compound_t     *fc,
                 vrna_recursion_status_f  cb);


int
vrna_gr_reset(vrna_fold_compound_t *fc);


/**
 * @}
 */

#endif
