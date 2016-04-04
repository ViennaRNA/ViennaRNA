/**********************************************/
/* BEGIN interface for structure constraints  */
/**********************************************/


%extend vrna_fold_compound_t {

  void constraints_add(const char *constraint, unsigned int options=VRNA_OPTION_MFE){
    vrna_constraints_add($self,constraint, options);
  }
}

%ignore print_tty_constraint;
%ignore print_tty_constraint_full;
%ignore constrain_ptypes;

%constant int DECOMP_PAIR_HP          = (int)VRNA_DECOMP_PAIR_HP;
%constant int DECOMP_PAIR_IL          = (int)VRNA_DECOMP_PAIR_IL;
%constant int DECOMP_PAIR_ML          = (int)VRNA_DECOMP_PAIR_ML;
%constant int DECOMP_ML_ML_ML         = (int)VRNA_DECOMP_ML_ML_ML;
%constant int DECOMP_ML_STEM          = (int)VRNA_DECOMP_ML_STEM;
%constant int DECOMP_ML_ML            = (int)VRNA_DECOMP_ML_ML;
%constant int DECOMP_ML_UP            = (int)VRNA_DECOMP_ML_UP;
%constant int DECOMP_ML_ML_STEM       = (int)VRNA_DECOMP_ML_ML_STEM;
%constant int DECOMP_ML_COAXIAL       = (int)VRNA_DECOMP_ML_COAXIAL;
%constant int DECOMP_EXT_EXT          = (int)VRNA_DECOMP_EXT_EXT;
%constant int DECOMP_EXT_UP           = (int)VRNA_DECOMP_EXT_UP;
%constant int DECOMP_EXT_STEM         = (int)VRNA_DECOMP_EXT_STEM;
%constant int DECOMP_EXT_EXT_EXT      = (int)VRNA_DECOMP_EXT_EXT_EXT;
%constant int DECOMP_EXT_STEM_EXT     = (int)VRNA_DECOMP_EXT_STEM_EXT;
%constant int DECOMP_EXT_STEM_OUTSIDE = (int)VRNA_DECOMP_EXT_STEM_OUTSIDE;
%constant int DECOMP_EXT_EXT_STEM     = (int)VRNA_DECOMP_EXT_EXT_STEM;
%constant int DECOMP_EXT_EXT_STEM1    = (int)VRNA_DECOMP_EXT_EXT_STEM1;

%include  <ViennaRNA/constraints.h>

