/**********************************************/
/* BEGIN interface for structure constraints  */
/**********************************************/


%extend vrna_fold_compound_t {
	
  void constraints_add(const char *constraint, unsigned int options=VRNA_OPTION_MFE)
  {
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


%constant unsigned int CONSTRAINT_DB            = VRNA_CONSTRAINT_DB;
%constant unsigned int CONSTRAINT_DB_ENFORCE_BP = VRNA_CONSTRAINT_DB_ENFORCE_BP;
%constant unsigned int CONSTRAINT_DB_PIPE       = VRNA_CONSTRAINT_DB_PIPE;
%constant unsigned int CONSTRAINT_DB_DOT        = VRNA_CONSTRAINT_DB_DOT;
%constant unsigned int CONSTRAINT_DB_X          = VRNA_CONSTRAINT_DB_X;
%constant unsigned int CONSTRAINT_DB_ANG_BRACK  = VRNA_CONSTRAINT_DB_ANG_BRACK;
%constant unsigned int CONSTRAINT_DB_RND_BRACK  = VRNA_CONSTRAINT_DB_RND_BRACK;
%constant unsigned int CONSTRAINT_DB_INTRAMOL   = VRNA_CONSTRAINT_DB_INTRAMOL;
%constant unsigned int CONSTRAINT_DB_INTERMOL   = VRNA_CONSTRAINT_DB_INTERMOL;
%constant unsigned int CONSTRAINT_DB_GQUAD      = VRNA_CONSTRAINT_DB_GQUAD;
%constant unsigned int CONSTRAINT_DB_DEFAULT    = VRNA_CONSTRAINT_DB_DEFAULT;

%constant int CONSTRAINT_CONTEXT_EXT_LOOP     = (int)VRNA_CONSTRAINT_CONTEXT_EXT_LOOP;
%constant int CONSTRAINT_CONTEXT_HP_LOOP      = (int)VRNA_CONSTRAINT_CONTEXT_HP_LOOP;
%constant int CONSTRAINT_CONTEXT_INT_LOOP     = (int)VRNA_CONSTRAINT_CONTEXT_INT_LOOP;
%constant int CONSTRAINT_CONTEXT_INT_LOOP_ENC = (int)VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
%constant int CONSTRAINT_CONTEXT_MB_LOOP      = (int)VRNA_CONSTRAINT_CONTEXT_MB_LOOP;
%constant int CONSTRAINT_CONTEXT_MB_LOOP_ENC  = (int)VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;
%constant int CONSTRAINT_CONTEXT_ENFORCE      = (int)VRNA_CONSTRAINT_CONTEXT_ENFORCE;
%constant int CONSTRAINT_CONTEXT_NO_REMOVE    = (int)VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;
%constant int CONSTRAINT_CONTEXT_ALL_LOOPS    = (int)VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;




%include  <ViennaRNA/constraints.h>

