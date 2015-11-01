/**********************************************/
/* BEGIN interface for structure constraints  */
/**********************************************/

%ignore vrna_hc_t;
%ignore vrna_hc_s;
%ignore vrna_sc_t;
%ignore vrna_sc_s;
%ignore vrna_message_constraint_options;
%ignore vrna_message_constraint_options_all;
%ignore vrna_constraints_add;
%ignore vrna_hc_init;
%ignore vrna_hc_add_up;
%ignore vrna_hc_add_bp;
%ignore vrna_hc_add_bp_nonspecific;
%ignore vrna_hc_free;
%ignore vrna_hc_add_f;
%ignore vrna_hc_add_pre;
%ignore vrna_hc_add_post;
%ignore vrna_sc_init;
%ignore vrna_sc_add_bp;
%ignore vrna_sc_add_up;
%ignore vrna_sc_remove;
%ignore vrna_sc_free;
%ignore vrna_sc_add_SHAPE_deigan;
%ignore vrna_sc_add_SHAPE_deigan_ali;
%ignore vrna_sc_add_SHAPE_zarringhalam;
%ignore vrna_sc_SHAPE_parse_method;
%ignore vrna_sc_SHAPE_to_pr;
%ignore vrna_sc_add_f;
%ignore vrna_sc_add_bt;
%ignore vrna_sc_add_exp_f;
%ignore vrna_sc_add_pre;
%ignore vrna_sc_add_post;
%ignore print_tty_constraint;
%ignore print_tty_constraint_full;
%ignore constrain_ptypes;

%include  "../src/ViennaRNA/constraints.h"
