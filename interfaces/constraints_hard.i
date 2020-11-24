/**********************************************/
/* BEGIN interface for hard constraints       */
/**********************************************/

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") hc_add_up;
%feature("kwargs") hc_add_up;
%feature("autodoc") hc_add_bp_nonspecific;
%feature("kwargs") hc_add_bp_nonspecific;
%feature("autodoc") hc_add_bp;
%feature("kwargs") hc_add_bp;
%feature("autodoc") hc_add_from_db;
%feature("kwargs") hc_add_from_db;
#endif

  void
  hc_init()
  {
    vrna_hc_init($self);
  }

  /* Make a certain nucleotide unpaired */
  void
  hc_add_up(int          i,
            unsigned int option = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
  {
    vrna_hc_add_up($self,i, (unsigned char)option);
  }

  /* Enforce a nucleotide to be paired (upstream/downstream) */
  void
  hc_add_bp_nonspecific(int           i,
                        int           d,
                        unsigned int  option = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
  {
    vrna_hc_add_bp_nonspecific($self,i,d, (unsigned char)option);
  }
  
  /* Favorize/Enforce  a certain base pair (i,j) */
  void
  hc_add_bp(int           i,
            int           j,
            unsigned int  option = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
  {
    vrna_hc_add_bp($self,i,j,(unsigned char)option);
  } 
  
  int
  hc_add_from_db(const char   *constraint,
                 unsigned int options = VRNA_CONSTRAINT_DB_DEFAULT)
  {
    return vrna_hc_add_from_db($self,constraint, options);
  }
}

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

%constant unsigned int CONSTRAINT_CONTEXT_EXT_LOOP        = VRNA_CONSTRAINT_CONTEXT_EXT_LOOP;
%constant unsigned int CONSTRAINT_CONTEXT_HP_LOOP         = VRNA_CONSTRAINT_CONTEXT_HP_LOOP;
%constant unsigned int CONSTRAINT_CONTEXT_INT_LOOP        = VRNA_CONSTRAINT_CONTEXT_INT_LOOP;
%constant unsigned int CONSTRAINT_CONTEXT_INT_LOOP_ENC    = VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
%constant unsigned int CONSTRAINT_CONTEXT_MB_LOOP         = VRNA_CONSTRAINT_CONTEXT_MB_LOOP;
%constant unsigned int CONSTRAINT_CONTEXT_MB_LOOP_ENC     = VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;
%constant unsigned int CONSTRAINT_CONTEXT_ENFORCE         = VRNA_CONSTRAINT_CONTEXT_ENFORCE;
%constant unsigned int CONSTRAINT_CONTEXT_NO_REMOVE       = VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;
%constant unsigned int CONSTRAINT_CONTEXT_ALL_LOOPS       = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
%constant unsigned int CONSTRAINT_CONTEXT_NONE            = VRNA_CONSTRAINT_CONTEXT_NONE;
%constant unsigned int CONSTRAINT_CONTEXT_CLOSING_LOOPS   = VRNA_CONSTRAINT_CONTEXT_CLOSING_LOOPS;
%constant unsigned int CONSTRAINT_CONTEXT_ENCLOSED_LOOPS  = VRNA_CONSTRAINT_CONTEXT_ENCLOSED_LOOPS;

%include  <ViennaRNA/constraints/hard.h>
