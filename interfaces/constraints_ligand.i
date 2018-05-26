/**********************************************/
/* BEGIN interface for structure ligand constraints */
/**********************************************/

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") sc_add_hi_motif;
%feature("kwargs") sc_add_hi_motif;
#endif

  int sc_add_hi_motif(const char *seq,
                      const char *structure,
                      FLT_OR_DBL energy,
                      unsigned int options=VRNA_OPTION_DEFAULT){

    return vrna_sc_add_hi_motif($self,seq,structure,energy,options);
  }
}

%include  <ViennaRNA/constraints/ligand.h>
