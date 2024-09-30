/**********************************************/
/* BEGIN interface for structure soft constraints */
/**********************************************/

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") sc_set_bp;
%feature("kwargs") sc_set_bp;
%feature("autodoc") sc_set_up;
%feature("kwargs") sc_set_up;
%feature("autodoc") sc_set_stack;
%feature("kwargs") sc_set_stack;
%feature("autodoc") sc_add_stack;
%feature("kwargs") sc_add_stack;
#endif

  void
  sc_remove()
  {
    vrna_sc_remove($self);
  }

  void
  sc_init()
  {
    vrna_sc_init($self);
  }
  
  int
  sc_add_up(unsigned int  i,
            double        energy,
            unsigned int  options = VRNA_OPTION_DEFAULT)
  {
    return vrna_sc_add_up($self, i, energy, options);
  }

  int
  sc_add_up(std::vector<double> constraints,
            unsigned int        options = VRNA_OPTION_DEFAULT)
  {
    std::vector<double>::iterator it;
    unsigned int  i = 1;
    int           ret = 1;
    it = constraints.begin();
    for(it++; it != constraints.end(); it++, i++){
      ret &= (vrna_sc_add_up($self, i, *it, options)) ? 1 : 0;
    }

    return ret;
  }

  int
  sc_add_bp(unsigned int  i,
            unsigned int  j,
            double        energy,
            unsigned int  options = VRNA_OPTION_DEFAULT)
  {
    return vrna_sc_add_bp($self, i, j, energy, options);
  }

  int
  sc_add_bp(std::vector<std::vector<double> > constraints,
            unsigned int                      options = VRNA_OPTION_DEFAULT)
  {
    int ret = 1;
    for (size_t i = 1; i < constraints.size(); i++)
      for (size_t j = i + 1; j < constraints[i].size(); j++)
        if (constraints[i][j] != 0)
          ret &= (vrna_sc_add_bp($self, i, j, constraints[i][j], options)) ? 1 : 0;

    return ret;
  }

  int
  sc_set_bp(std::vector<std::vector<double> > constraints,
            unsigned int                      options = VRNA_OPTION_DEFAULT)
  {
    int ret = 0;

    /* make sure that the constraints matrix is large enough */
    FLT_OR_DBL **c = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * ($self->length + 1));
    for(unsigned int i = 0; i <= $self->length; i++)
      c[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * ($self->length + 1));

    /* copy input data (missing values have value 0 */
    for(unsigned int i = 0; (i < constraints.size()) && (i <= $self->length); i++)
      for(unsigned int j = i; (j < constraints[i].size()) && (j <= $self->length); j++)
        c[i][j] = (FLT_OR_DBL)constraints[i][j];

    ret = vrna_sc_set_bp($self, (const FLT_OR_DBL **)c, options);

    /* cleanup */
    for(unsigned int i = 0; i <= $self->length; i++)
      free(c[i]);
    free(c);

    return ret;
  }

  int
  sc_set_up(std::vector<double> constraints,
            unsigned int        options = VRNA_OPTION_DEFAULT)
  {
    std::vector<FLT_OR_DBL>  v;
    transform(constraints.begin(), constraints.end(), std::back_inserter(v), convert_vecdbl2vecFLR_OR_DBL);
    return vrna_sc_set_up($self, (const FLT_OR_DBL *)&v[0], options);
  }

  int
  sc_set_stack(std::vector<double> constraints,
               unsigned int        options = VRNA_OPTION_DEFAULT)
  {
    std::vector<FLT_OR_DBL>  v;
    transform(constraints.begin(), constraints.end(), std::back_inserter(v), convert_vecdbl2vecFLR_OR_DBL);
    return vrna_sc_set_stack($self, (const FLT_OR_DBL *)&v[0], options);
  }

  int
  sc_set_stack(std::vector<std::vector<double> >  constraints,
               unsigned int                       options = VRNA_OPTION_DEFAULT)
  {
    int ret = 0;

    if ($self->type == VRNA_FC_TYPE_COMPARATIVE) {
      FLT_OR_DBL **c = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * ($self->n_seq + 1));

      for(unsigned int s = 0; s <= $self->n_seq; s++)
        c[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * ($self->length + 1));

      /* copy input data (missing values have value 0 */
      for(unsigned int s = 0; (s < constraints.size()) && (s <= $self->n_seq); s++)
        for(unsigned int i = 1; (i < constraints[s].size()) && (i <= $self->length); i++)
          c[s][i] = (FLT_OR_DBL)constraints[s][i];

      ret = vrna_sc_set_stack_comparative($self, (const FLT_OR_DBL **)c, options);

      /* cleanup */
      for(unsigned int i = 0; i <= $self->length; i++)
        free(c[i]);
      free(c);
    }

    return ret;
  }

  int
  sc_add_stack(unsigned int i,
               double       energy,
               unsigned int options = VRNA_OPTION_DEFAULT)
  {
    return vrna_sc_add_stack($self, i, energy, options);
  }

  int
  sc_add_stack(std::vector<unsigned int>  i,
               std::vector<double>        energies,
               unsigned int               options = VRNA_OPTION_DEFAULT)
  {
    std::vector<FLT_OR_DBL>  v;
    transform(energies.begin(), energies.end(), std::back_inserter(v), convert_vecdbl2vecFLR_OR_DBL);
    return vrna_sc_add_stack_comparative($self, (unsigned int*)&i[0], (const FLT_OR_DBL *)&v[0], options);
  }
}

%include  <ViennaRNA/constraints/soft.h>
