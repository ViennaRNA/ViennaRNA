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
  
  void sc_add_up(int i, double energy, unsigned int options=VRNA_OPTION_DEFAULT){
    vrna_sc_add_up($self, i, energy, options);
  }

  void sc_add_up( std::vector<double> constraints,
                  unsigned int options=VRNA_OPTION_DEFAULT){

    std::vector<double>::iterator it;
    int i = 1;
    it = constraints.begin();
    for(it++; it != constraints.end(); it++, i++){
      vrna_sc_add_up($self, i, *it, options);
    }
  }

  void sc_add_bp(int i, int j, double energy, unsigned int options=VRNA_OPTION_DEFAULT){
    vrna_sc_add_bp($self, i, j, energy, options);
  }

  void sc_add_bp( std::vector<std::vector<double> > constraints,
                  unsigned int options=VRNA_OPTION_DEFAULT){

    std::vector<std::vector<double> >::iterator it;
    std::vector<double>::iterator it2;
    int i, j;

    i = 1;
    it = constraints.begin();
    for(it++; it != constraints.end(); it++, i++){
      it2 = (*it).begin();
      j   = 1;
      for(it2++; it2 != (*it).end(); it2++, j++){
        vrna_sc_add_bp($self, i, j, *it2, options);
      }
    }
  }

  void sc_set_bp( std::vector<std::vector<double> > constraints,
                  unsigned int options=VRNA_OPTION_DEFAULT){

    /* make sure that the constraints matrix is large enough */
    FLT_OR_DBL **c = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * ($self->length + 1));
    for(unsigned int i = 0; i <= $self->length; i++)
      c[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * ($self->length + 1));

    /* copy input data (missing values have value 0 */
    for(unsigned int i = 0; (i < constraints.size()) && (i <= $self->length); i++)
      for(unsigned int j = i; (j < constraints[i].size()) && (j <= $self->length); j++)
        c[i][j] = (FLT_OR_DBL)constraints[i][j];

    vrna_sc_set_bp($self, (const FLT_OR_DBL **)c, options);

    /* cleanup */
    for(unsigned int i = 0; i <= $self->length; i++)
      free(c[i]);
    free(c);
  }

  void sc_set_up( std::vector<double> constraints,
                  unsigned int options=VRNA_OPTION_DEFAULT){

    std::vector<FLT_OR_DBL>  v;
    transform(constraints.begin(), constraints.end(), std::back_inserter(v), convert_vecdbl2vecFLR_OR_DBL);
    vrna_sc_set_up($self, (const FLT_OR_DBL *)&v[0], options);
  }


}

%include  <ViennaRNA/constraints_soft.h>
