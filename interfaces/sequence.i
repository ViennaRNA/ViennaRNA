/******************************************************/
/* BEGIN interface for Combinatorics Implementations  */
/******************************************************/

%extend vrna_fold_compound_t {

  int
  sequence_add(std::string sequence,
               unsigned int options = VRNA_SEQUENCE_RNA)
  {
    return vrna_sequence_add($self,
                             sequence.c_str(),
                             options);
  }


  int
  sequence_remove(unsigned int i)
  {
    return vrna_sequence_remove($self,
                                i);
  }


  void
  sequence_remove_all(void)
  {
    vrna_sequence_remove_all($self);
  }


  void
  sequence_prepare(void)
  {
    vrna_sequence_prepare($self);
  }
}


/*
 *  Rename all the preprocessor macros defined in sequence.h
 *  (wrapped as constants)
 */
%constant unsigned int SEQUENCE_RNA = VRNA_SEQUENCE_RNA;
%constant unsigned int SEQUENCE_DNA = VRNA_SEQUENCE_DNA;

%include <ViennaRNA/sequence.h>
