/**********************************************/
/* BEGIN interface for energy evaluation      */
/**********************************************/

%ignore vrna_eval_structure;
%ignore vrna_eval_covar_structure;
%ignore vrna_eval_structure_simple;
%ignore vrna_eval_structure_verbose;
%ignore vrna_eval_structure_simple_verbose;
%ignore vrna_eval_structure_pt;
%ignore vrna_eval_structure_pt_simple;
%ignore vrna_eval_structure_pt_verbose;
%ignore vrna_eval_structure_pt_simple_verbose;
%ignore vrna_eval_loop_pt;
%ignore vrna_eval_move;
%ignore vrna_eval_move_pt;
%ignore vrna_eval_move_pt_simple;
%ignore energy_of_struct_par;
%ignore energy_of_circ_struct_par;
%ignore energy_of_gquad_struct_par;
%ignore energy_of_struct_pt_par;

%include  "../src/ViennaRNA/eval.h"
