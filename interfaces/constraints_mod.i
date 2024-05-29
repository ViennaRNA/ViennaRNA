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

%include  <ViennaRNA/constraints/soft_special.h>
