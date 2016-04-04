/**********************************************/
/* BEGIN interface for structure hard constraints */
/**********************************************/

%extend vrna_fold_compound_t {

  void hc_init()
  {
	  vrna_hc_init($self);
  }
  /*Make a certain nucleotide unpaired*/
  void hc_add_up(int i, char option=VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
  {
	  vrna_hc_add_up($self,i,option);
  }
  /*Enforce a nucleotide to be paired (upstream/downstream)*/
  void hc_add_bp_nonspecific(int i, int d, char option=VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
  {
	  vrna_hc_add_bp_nonspecific($self,i,d,option);
  }
  
  /*Favorize/Enforce  a certain base pair (i,j)*/
  void hc_add_bp(int i, int j, char option= VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
  {
	  vrna_hc_add_bp($self,i,j,option);
  } 
  
  
  int hc_add_from_db(const char *constraint, unsigned int options=VRNA_CONSTRAINT_DB_DEFAULT)
  {
	  return vrna_hc_add_from_db($self,constraint,options);
  }
  
}


%include  <ViennaRNA/constraints_hard.h>
