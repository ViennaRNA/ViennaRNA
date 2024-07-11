
typedef struct {
  vrna_grammar_rule_f     cb;
  vrna_grammar_bt_f       cb_bt;
  void                    *data;
  vrna_auxdata_prepare_f  prepare;
  vrna_auxdata_free_f     release;
} grammar_rule_wrap_t;


typedef struct {
  vrna_gr_inside_exp_f    cb;
  vrna_gr_outside_exp_f   cb_out;
  void                    *data;
  vrna_auxdata_prepare_f  prepare;
  vrna_auxdata_free_f     release;
} grammar_rule_exp_wrap_t;


struct vrna_gr_aux_s {
  vrna_array(grammar_rule_wrap_t) f;
  vrna_array(grammar_rule_wrap_t) c;
  vrna_array(grammar_rule_wrap_t) m;
  vrna_array(grammar_rule_wrap_t) m1;
  vrna_array(grammar_rule_wrap_t) aux;

  vrna_array(grammar_rule_exp_wrap_t) exp_f;
  vrna_array(grammar_rule_exp_wrap_t) exp_c;
  vrna_array(grammar_rule_exp_wrap_t) exp_m;
  vrna_array(grammar_rule_exp_wrap_t) exp_m1;
  vrna_array(grammar_rule_exp_wrap_t) exp_aux;

  vrna_array(vrna_recursion_status_f) cbs_status;
  vrna_array(void *)                  datas;
  vrna_array(vrna_auxdata_prepare_f)  prepare_datas;
  vrna_array(vrna_auxdata_free_f)     free_datas;
};
