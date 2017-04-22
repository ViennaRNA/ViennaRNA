/**********************************************/
/* BEGIN interface for hard constraints       */
/**********************************************/

%extend vrna_fold_compound_t {

  void hc_init(){
    vrna_hc_init($self);
  }

  /* Make a certain nucleotide unpaired */
  void hc_add_up(int i, int option=(unsigned char)VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS){
    vrna_hc_add_up($self,i, (char)option);
  }

  /* Enforce a nucleotide to be paired (upstream/downstream) */
  void hc_add_bp_nonspecific(int i, int d, int option=(unsigned char)VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS){
    vrna_hc_add_bp_nonspecific($self,i,d, (char)option);
  }
  
  /* Favorize/Enforce  a certain base pair (i,j) */
  void hc_add_bp(int i, int j, int option=(unsigned char)VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS){
    vrna_hc_add_bp($self,i,j,(char)option);
  } 
  
  int hc_add_from_db(const char *constraint, unsigned int options=VRNA_CONSTRAINT_DB_DEFAULT){
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

%constant int CONSTRAINT_CONTEXT_EXT_LOOP     = (unsigned char)VRNA_CONSTRAINT_CONTEXT_EXT_LOOP;
%constant int CONSTRAINT_CONTEXT_HP_LOOP      = (unsigned char)VRNA_CONSTRAINT_CONTEXT_HP_LOOP;
%constant int CONSTRAINT_CONTEXT_INT_LOOP     = (unsigned char)VRNA_CONSTRAINT_CONTEXT_INT_LOOP;
%constant int CONSTRAINT_CONTEXT_INT_LOOP_ENC = (unsigned char)VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
%constant int CONSTRAINT_CONTEXT_MB_LOOP      = (unsigned char)VRNA_CONSTRAINT_CONTEXT_MB_LOOP;
%constant int CONSTRAINT_CONTEXT_MB_LOOP_ENC  = (unsigned char)VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;
%constant int CONSTRAINT_CONTEXT_ENFORCE      = (unsigned char)VRNA_CONSTRAINT_CONTEXT_ENFORCE;
%constant int CONSTRAINT_CONTEXT_NO_REMOVE    = (unsigned char)VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;
%constant int CONSTRAINT_CONTEXT_ALL_LOOPS    = (unsigned char)VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

%include  <ViennaRNA/constraints_hard.h>
