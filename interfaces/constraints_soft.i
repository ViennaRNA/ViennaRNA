/**********************************************/
/* BEGIN interface for structure soft constraints */
/**********************************************/

%extend vrna_fold_compound_t {

  void sc_remove(){
    vrna_sc_remove($self);
  }

  void sc_init(){
    vrna_sc_init($self);
  }
  
  void sc_add_bp( std::vector<std::vector<double> > constraints,
                  unsigned int options=VRNA_OPTION_DEFAULT){

    /* make sure that the constraints matrix is large enough */
    FLT_OR_DBL **c = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * ($self->length + 1));
    for(unsigned int i = 0; i <= $self->length; i++)
      c[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * ($self->length + 1));

    /* copy input data (missing values have value 0 */
    for(unsigned int i = 0; (i < constraints.size()) && (i <= $self->length); i++)
      for(unsigned int j = i; (j < constraints[i].size()) && (j <= $self->length); j++)
        c[i][j] = (FLT_OR_DBL)constraints[i][j];

    vrna_sc_add_bp($self, (const FLT_OR_DBL **)c, options);

    /* cleanup */
    for(unsigned int i = 0; i <= $self->length; i++)
      free(c[i]);
    free(c);
  }

  void sc_add_up( std::vector<double> constraints,
                  unsigned int options=VRNA_OPTION_DEFAULT){

    std::vector<FLT_OR_DBL>  v;
    transform(constraints.begin(), constraints.end(), std::back_inserter(v), convert_vecdbl2vecFLR_OR_DBL);
    vrna_sc_add_up($self, (const FLT_OR_DBL *)&v[0], options);
  }
}

%include  <ViennaRNA/constraints_soft.h>
