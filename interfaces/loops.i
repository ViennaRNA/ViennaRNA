/**********************************************/
/* BEGIN interface for loop energies          */
/**********************************************/


%rename(E_ext_stem)       vrna_E_ext_stem;
%rename(exp_E_ext_stem)   vrna_exp_E_ext_stem;


%extend vrna_fold_compound_t{

#ifdef SWIGPYTHON
%feature("autodoc") E_ext_loop;
%feature("kwargs") E_ext_loop;
%feature("autodoc") eval_hp_loop;
%feature("kwargs") eval_hp_loop;
%feature("autodoc") eval_int_loop;
%feature("kwargs") eval_int_loop;
#endif

  int
  eval_ext_stem(int i,
                int j)
  {
    return vrna_eval_ext_stem($self, i, j);
  }

  int
  E_hp_loop(int i,
            int j)
  {
    return vrna_E_hp_loop($self, i, j);
  }

  int
  E_ext_hp_loop(int i,
                int j)
  {
    return vrna_E_ext_hp_loop($self, i, j);
  }

  int
  eval_ext_hp_loop(int i,
                   int j)
  {
    return vrna_eval_ext_hp_loop($self, i, j);
  }

  int
  eval_hp_loop(int i,
               int j)
  {
    return vrna_eval_hp_loop($self, i, j);
  }

  double
  exp_E_hp_loop(int i,
                int j)
  {
    return (double)vrna_exp_E_hp_loop($self, i, j);
  }

  int
  E_int_loop(int i,
             int j)
  {
    return vrna_E_int_loop($self, i, j);
  }

  int
  eval_int_loop(int i,
                int j,
                int k,
                int l)
  {
    return vrna_eval_int_loop($self, i, j, k, l);
  }

  %apply int *OUTPUT { int *ip, int *iq };

  int
  E_ext_int_loop(int i,
                 int j,
                 int *ip,
                 int *iq)
  {
    return vrna_E_ext_int_loop($self, i, j, ip, iq);
  }

  %clear  int *ip, int *iq;

  int
  E_stack(int i,
          int j)
  {
    return vrna_E_stack($self, i, j);
  }

  double
  exp_E_int_loop(int i,
                 int j)
  {
    return (double)vrna_exp_E_int_loop($self, i, j);
  }

  double
  exp_E_interior_loop(int i,
                      int j,
                      int k,
                      int l)
  {
    return (double)vrna_exp_E_interior_loop($self, i, j, k, l);
  }

  double
  exp_E_ext_stem(int i,
                 int j)
  {
    unsigned int type;
    int enc5, enc3;
    enc5 = enc3 = -1;
    
    type = vrna_get_ptype_md($self->sequence_encoding2[i],
                             $self->sequence_encoding2[j],
                             &($self->params->model_details));

    if (i > 1)
      enc5 = $self->sequence_encoding[i - 1];
    if (j < $self->length)
      enc3 = $self->sequence_encoding[j + 1];

    return (double)vrna_exp_E_ext_stem(type,
                                       enc5,
                                       enc3,
                                       $self->exp_params);
  }
}

%include  <ViennaRNA/loops/external.h>
%include  <ViennaRNA/loops/hairpin.h>
%include  <ViennaRNA/loops/internal.h>
%include  <ViennaRNA/loops/multibranch.h>
