/******************************************************/
/* BEGIN interface for Combinatorics Implementations  */
/******************************************************/

%extend vrna_fold_compound_t {

  int
  sequence_add(std::string  sequence,
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

%rename (seq_encode) my_seq_encode;

#ifdef SWIGPYTHON
%feature("autodoc") my_seq_encode;
%feature("kwargs") my_seq_encode;
#endif

%{
#include <vector>

  std::vector<int>
  my_seq_encode(std::string sequence,
                vrna_md_t   *md_p = NULL)
  {
    short             *s;
    int               n;
    std::vector<int>  encoding;
    vrna_md_t         md;

    if (!md_p) {
      vrna_md_set_default(&md);
      md_p = &md;
    }

    n = sequence.length();
    s = vrna_seq_encode(sequence.c_str(),
                        md_p);
  
    encoding.push_back(n);
    for (int i = 1; i <= n; i++)
      encoding.push_back(s[i]);

    free(s);

    return encoding;
  }

%}

std::vector<int>
my_seq_encode(std::string sequence,
              vrna_md_t   *md_p = NULL);

/*
 *  Rename all the preprocessor macros defined in sequence.h
 *  (wrapped as constants)
 */
%constant unsigned int SEQUENCE_RNA = VRNA_SEQUENCE_RNA;
%constant unsigned int SEQUENCE_DNA = VRNA_SEQUENCE_DNA;

%include <ViennaRNA/sequence.h>
