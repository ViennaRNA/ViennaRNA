
struct vrna_gr_aux_s {
  vrna_grammar_rule_f         cb_aux_f;
  vrna_grammar_rule_f         cb_aux_c;
  vrna_grammar_rule_f         cb_aux_m;
  vrna_grammar_rule_f         cb_aux_m1;
  vrna_grammar_rule_f_aux     cb_aux;

  vrna_grammar_bt_f           cb_bt_f;
  vrna_grammar_bt_f           cb_bt_c;
  vrna_grammar_bt_f           cb_bt_m;
  vrna_grammar_bt_f           cb_bt_m1;
  vrna_grammar_bt_f           cb_bt_aux;


  vrna_grammar_rule_f_exp     cb_aux_exp_f;
  vrna_grammar_rule_f_exp     cb_aux_exp_c;
  vrna_grammar_rule_f_exp     cb_aux_exp_m;
  vrna_grammar_rule_f_exp     cb_aux_exp_m1;
  vrna_grammar_rule_f_aux_exp cb_aux_exp;

  vrna_recursion_status_f     cb_proc; /**< @brief A callback for pre- and post-processing of auxiliary grammar rules */
  void                        *data;
  vrna_auxdata_prepare_f      prepare_data;
  vrna_auxdata_free_f         free_data;

  vrna_array(unsigned int)    bt_flags;
};

