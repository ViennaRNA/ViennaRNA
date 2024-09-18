/**********************************************/
/* BEGIN interface for Backtracking           */
/**********************************************/

/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::backtrack;

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") backtrack;
#endif

  char *backtrack(unsigned int length,
                  float *OUTPUT) {
    char *structure = (char *)vrna_alloc(sizeof(char) * (length + 1));
    *OUTPUT = vrna_backtrack5($self, length, structure);
    return structure;
  }

  char *backtrack(float *OUTPUT) {
    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_backtrack5($self, $self->length, structure);
    return structure;
  }

}

%include  <ViennaRNA/backtrack/global.h>
%include  <ViennaRNA/backtrack/exterior.h>
%include  <ViennaRNA/backtrack/hairpin.h>
%include  <ViennaRNA/backtrack/internal.h>
%include  <ViennaRNA/backtrack/multibranch.h>
%include  <ViennaRNA/backtrack/gquad.h>
