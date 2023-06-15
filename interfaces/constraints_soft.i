/**********************************************/
/* BEGIN interface for structure soft constraints */
/**********************************************/

%rename(sc_mod_param) vrna_sc_mod_param_s;

typedef struct {} vrna_sc_mod_param_s;

/* no default constructor / destructor */
%nodefaultctor vrna_sc_mod_param_s;
%nodefaultdtor vrna_sc_mod_param_s;

%extend vrna_sc_mod_param_s {
  vrna_sc_mod_param_s(std::string json,
                      vrna_md_t   *md = NULL)
  {
    /* try interpreting the json input as file name */
    vrna_sc_mod_param_s *obj = vrna_sc_mod_read_from_jsonfile(json.c_str(), md);

    if (obj == NULL)
      /* second guess is that the json string is actual json data */
      obj = vrna_sc_mod_read_from_json(json.c_str(), md);

    return obj;
  }

  ~vrna_sc_mod_param_s()
  {
    vrna_sc_mod_parameters_free($self);
  }
}

%rename (sc_mod_read_from_jsonfile) my_sc_mod_read_from_jsonfile;
%rename (sc_mod_read_from_json) my_sc_mod_read_from_json;
%rename (sc_mod_parameters_free) vrna_sc_mod_parameters_free;

%{
  vrna_sc_mod_param_t
  my_sc_mod_read_from_jsonfile(std::string  filename,
                               vrna_md_t    *md = NULL)
  {
    return vrna_sc_mod_read_from_jsonfile(filename.c_str(), md);
  }

  vrna_sc_mod_param_t
  my_sc_mod_read_from_json(std::string  json,
                           vrna_md_t    *md = NULL)
  {
    return vrna_sc_mod_read_from_json(json.c_str(), md);
  }
%}

vrna_sc_mod_param_t my_sc_mod_read_from_jsonfile(std::string filename, vrna_md_t *md = NULL);
vrna_sc_mod_param_t my_sc_mod_read_from_json(std::string json, vrna_md_t *md = NULL);


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

%feature("kwargs") sc_mod_json;
%feature("kwargs") sc_mod_jsonfile;
%feature("kwargs") sc_mod;
%feature("kwargs") sc_mod_m6A;
%feature("kwargs") sc_mod_pseudouridine;
%feature("kwargs") sc_mod_inosine;
%feature("kwargs") sc_mod_7DA;
%feature("kwargs") sc_mod_purine;
%feature("kwargs") sc_mod_dihydrouridine;
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
  sc_add_up(int          i,
            double       energy,
            unsigned int options = VRNA_OPTION_DEFAULT)
  {
    return vrna_sc_add_up($self, i, energy, options);
  }

  int
  sc_add_up(std::vector<double> constraints,
            unsigned int        options = VRNA_OPTION_DEFAULT)
  {
    std::vector<double>::iterator it;
    int i = 1, ret = 1;
    it = constraints.begin();
    for(it++; it != constraints.end(); it++, i++){
      ret &= (vrna_sc_add_up($self, i, *it, options)) ? 1 : 0;
    }

    return ret;
  }

  int
  sc_add_bp(int          i,
            int          j,
            double       energy,
            unsigned int options = VRNA_OPTION_DEFAULT)
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
  sc_add_stack(int           i,
               double        energy,
               unsigned int  options = VRNA_OPTION_DEFAULT)
  {
    return vrna_sc_add_stack($self, i, energy, options);
  }

  int
  sc_add_stack(int                  i,
               std::vector<double>  energies,
               unsigned int         options = VRNA_OPTION_DEFAULT)
  {
    std::vector<FLT_OR_DBL>  v;
    transform(energies.begin(), energies.end(), std::back_inserter(v), convert_vecdbl2vecFLR_OR_DBL);
    return vrna_sc_add_stack_comparative($self, i, (const FLT_OR_DBL *)&v[0], options);
  }

  int
  sc_mod_json(std::string json,
              std::vector<unsigned int> modification_sites,
              unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_json($self, json.c_str(), &modification_sites[0], options);
  }

  int
  sc_mod_jsonfile(std::string jsonfile,
                  std::vector<unsigned int> modification_sites,
                  unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_json($self, jsonfile.c_str(), &modification_sites[0], options);
  }

  int
  sc_mod(const vrna_sc_mod_param_t params,
         std::vector<unsigned int> modification_sites,
         unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod($self, params, &modification_sites[0], options);
  }

  int
  sc_mod_m6A(std::vector<unsigned int> modification_sites,
             unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_m6A($self, &modification_sites[0], options);
  }

  int
  sc_mod_pseudouridine(std::vector<unsigned int> modification_sites,
                       unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_pseudouridine($self, &modification_sites[0], options);
  }

  int
  sc_mod_inosine(std::vector<unsigned int> modification_sites,
                 unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_inosine($self, &modification_sites[0], options);
  }

  int
  sc_mod_7DA(std::vector<unsigned int> modification_sites,
             unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_7DA($self, &modification_sites[0], options);
  }

  int
  sc_mod_purine(std::vector<unsigned int> modification_sites,
                unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_purine($self, &modification_sites[0], options);
  }

  int
  sc_mod_dihydrouridine(std::vector<unsigned int> modification_sites,
                        unsigned int options = VRNA_SC_MOD_DEFAULT) {
    modification_sites.push_back(0); /* end marker for C-implementation */
    return vrna_sc_mod_dihydrouridine($self, &modification_sites[0], options);
  }

}

%constant unsigned int SC_MOD_CHECK_FALLBACK  = VRNA_SC_MOD_CHECK_FALLBACK;
%constant unsigned int SC_MOD_CHECK_UNMOD     = VRNA_SC_MOD_CHECK_UNMOD;
%constant unsigned int SC_MOD_SILENT          = VRNA_SC_MOD_SILENT;
%constant unsigned int SC_MOD_DEFAULT         = VRNA_SC_MOD_DEFAULT;

%include  <ViennaRNA/constraints/soft.h>
%include  <ViennaRNA/constraints/soft_special.h>
