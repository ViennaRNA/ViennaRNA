/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT


#ifdef VRNA_BACKWARD_COMPAT

#define ggo_get_temperature(ggostruct, dest) ({ \
    if (ggostruct.temp_given) \
      dest = ggostruct.temp_arg; \
    temperature = (double)dest; \
  })

#define ggo_get_dangles(ggostruct, dest) ({ \
    if (ggostruct.dangles_given) \
      dest = ggostruct.dangles_arg; \
    dangles = (int)dest; \
  })

#define ggo_get_special_hp(ggostruct, dest) ({ \
    if (ggostruct.noTetra_given) \
      dest = 0; \
    tetra_loop = (int)dest; \
  })

#define ggo_get_noLP(ggostruct, dest) ({ \
    if (ggostruct.noLP_given) \
      dest = 1; \
    noLonelyPairs = (int)dest; \
  })

#define ggo_get_noGU(ggostruct, dest) ({ \
    if (ggostruct.noGU_given) \
      dest = 1; \
    noGU = (int)dest; \
  })

#define ggo_get_noGUclosure(ggostruct, dest) ({ \
    if (ggostruct.noClosingGU_given) \
      dest = 1; \
    no_closingGU = (int)dest; \
  })

#define ggo_get_gquad(ggostruct, dest) ({ \
    if (ggostruct.gquad_given) \
      dest = 1; \
    gquad = (int)dest; \
  })

#define ggo_get_energyModel(ggostruct, dest) ({ \
    if (ggostruct.energyModel_given) \
      dest = ggostruct.energyModel_arg; \
    energy_set = (int)dest; \
  })

#define ggo_get_maxBPspan(ggostruct, dest)  ({ \
    if (ggostruct.maxBPspan_given) \
      dest = ggostruct.maxBPspan_arg; \
    max_bp_span = (int)dest; \
  })


#else

#define ggo_get_temperature(ggostruct, dest) ({ \
    if (ggostruct.temp_given) \
      dest = ggostruct.temp_arg; \
  })

#define ggo_get_dangles(ggostruct, dest) ({ \
    if (ggostruct.dangles_given) \
      dest = ggostruct.dangles_arg; \
  })

#define ggo_get_special_hp(ggostruct, dest) ({ \
    if (ggostruct.noTetra_given) \
      dest = 0; \
  })

#define ggo_get_noLP(ggostruct, dest) ({ \
    if (ggostruct.noLP_given) \
      dest = 1; \
  })

#define ggo_get_noGUclosure(ggostruct, dest) ({ \
    if (ggostruct.noClosingGU_given) \
      dest = 1; \
  })

#define ggo_get_gquad(ggostruct, dest) ({ \
    if (ggostruct.gquad_given) \
      dest = 1; \
  })

#define ggo_get_energyModel(ggostruct, dest) ({ \
    if (ggostruct.energyModel_given) \
      dest = ggostruct.energyModel_arg; \
  })

#define ggo_get_maxBPspan(ggostruct, dest)  ({ \
    if (ggostruct.maxBPspan_given) \
      dest = ggostruct.maxBPspan_arg; \
  })

#endif


#define ggo_get_pfScale(ggostruct, dest) ({ \
    if (ggostruct.pfScale_given) \
      dest = ggostruct.pfScale_arg; \
  })

#define ggo_get_betaScale(ggostruct, dest) ({ \
    if (ggostruct.betaScale_given) \
      dest = ggostruct.betaScale_arg; \
  })

/* assume RNA sequence to be circular */
#define ggo_get_circ(ggostruct, dest) ({ \
    if (ggostruct.circ_given) \
      dest = 1; \
  })

/* Allow other pairs in addition to the usual AU,GC,and GU pairs */
#define ggo_get_nsp(ggostruct, dest) ({ \
    if (ggostruct.nsp_given) \
      dest = strdup(ggostruct.nsp_arg); \
  })

#define ggo_get_read_paramFile(ggostruct) ({ \
    /* take another energy parameter set */ \
    if (ggostruct.paramFile_given) \
      read_parameter_file(ggostruct.paramFile_arg); \
  })


/*
 *  The following macro automatically sets a basic set of
 *  model details:
 *  - temperature
 *  - dangles
 *  - special_hp
 *  - gquad
 *  - energy_set
 */
#define ggo_get_md_eval(ggostruct, md) ({ \
    /* temperature */ \
    ggo_get_temperature(ggostruct, md.temperature); \
    /* dangles */ \
    ggo_get_dangles(ggostruct, md.dangles); \
    /* special_hp */ \
    ggo_get_special_hp(ggostruct, md.special_hp); \
    /* gquadruplex support */ \
    ggo_get_gquad(ggostruct, md.gquad); \
    /* set energy model */ \
    ggo_get_energyModel(ggostruct, md.energy_set); \
    /* Allow other pairs in addition to the usual AU,GC,and GU pairs */ \
    { char *ns_bases = NULL; \
      ggo_get_nsp(ggostruct, ns_bases); \
      if (ns_bases != NULL) \
        vrna_md_set_nonstandards(&md, ns_bases); \
    } \
    ggo_get_read_paramFile(ggostruct); \
  })


/*
 *  The following macro automatically sets a basic set of
 *  model details required for RNA structure prediction:
 *  - noLP
 *  - noGU
 *  - noGUclosure
 */
#define ggo_get_md_fold(ggostruct, md) ({ \
    /* do not allow weak pairs */ \
    ggo_get_noLP(ggostruct, md.noLP); \
    /* do not allow wobble pairs (GU) */ \
    ggo_get_noGU(ggostruct, md.noGU); \
    /* do not allow weak closing pairs (AU,GU) */ \
    ggo_get_noGUclosure(ggostruct, md.noGUclosure); \
    /* set maximum base pair span */ \
    ggo_get_maxBPspan(ggostruct, md.max_bp_span); \
  })


/*
 *  The following macro automatically sets a basic set of
 *  model details required for RNA structure prediction:
 *  - noLP
 *  - noGU
 *  - noGUclosure
 */
#define ggo_get_md_part(ggostruct, md) ({ \
    /* set pf scaling factor */ \
    ggo_get_pfScale(ggostruct, md.sfact); \
    ggo_get_betaScale(ggostruct, md.betaScale); \
  })


/*
 * HELPER MACRO FOR SHAPE REACTIVITY DATA INCORPORATION
 */

#define ggo_get_SHAPE(ggostruct, \
                      SHAPE_switch, \
                      SHAPE_files, \
                      SHAPE_method, \
                      SHAPE_conversion)  ({ \
    /* SHAPE reactivity data */ \
    if (ggostruct.shape_given) { \
      SHAPE_switch = 1; \
      SHAPE_files = strdup(ggostruct.shape_arg); \
      SHAPE_method = strdup(ggostruct.shapeMethod_arg); \
      SHAPE_conversion = strdup(ggostruct.shapeConversion_arg); \
    } else { \
      SHAPE_switch = 0; \
      SHAPE_files = NULL; \
      SHAPE_method = NULL; \
      SHAPE_conversion = NULL; \
    } \
  })


/*
 * HELPER FUNCTIONS FOR AUTO-ID FEATURE
 */

#define ggo_get_ID_manipulation(ggostruct, \
                                ID_auto_switch, \
                                ID_prefix, ID_prefix_default, \
                                ID_delimiter, ID_delimiter_default, \
                                ID_digits, ID_digits_default, \
                                ID_start, ID_start_default) ({ \
    if (ggostruct.auto_id_given) { ID_auto_switch = 1; } \
    if (ggostruct.id_prefix_given) { \
      ID_prefix = strdup(ggostruct.id_prefix_arg); \
      ID_auto_switch = 1; \
    } else { ID_prefix = strdup(ID_prefix_default); } \
    if (ggostruct.id_delim_given) { \
      ID_delimiter = strdup(ggostruct.id_delim_arg); \
    } else { ID_delimiter = strdup(ID_delimiter_default); } \
    if (ggostruct.id_digits_given) { \
      if ((ggostruct.id_digits_arg > 0) && (ggostruct.id_digits_arg < 19)) { \
        ID_digits = ggostruct.id_digits_arg; } \
      else { \
        vrna_message_warning("ID number digits out of allowed range! Using defaults..."); \
        ID_digits = ID_digits_default; \
      } \
    } else { ID_digits = ID_digits_default; } \
    if (ggostruct.id_start_given) { \
      ID_auto_switch = 1; \
      if ((ggostruct.id_start_arg >= 0) && (ggostruct.id_start_arg <= LONG_MAX)) { \
        ID_start = ggostruct.id_start_arg; \
      } else { \
        vrna_message_warning("ID number start out of allowed range! Using defaults..."); \
        ID_start = ID_start_default; \
      } \
    } else { ID_start = ID_start_default; } \
  })


#define ggo_get_constraints_settings(ggostruct, \
                                     constraint_switch, \
                                     constraint_file, \
                                     constraint_enforce, \
                                     constraint_batch)  ({ \
    /* structure constraint */ \
    if (ggostruct.constraint_given) { \
      constraint_switch = 1; \
      if (ggostruct.constraint_arg[0] != '\0') \
        constraint_file = strdup(ggostruct.constraint_arg); \
      else \
        constraint_file = NULL; \
    } else { constraint_switch = 0; constraint_file = NULL; } \
    /* enforce base pairs given in constraint string rather than weak enforce */ \
    if (ggostruct.enforceConstraint_given) \
      constraint_enforce = 1; \
    else \
      constraint_enforce = 0; \
    /* do batch jobs despite constraints read from input file */ \
    if (ggostruct.batch_given) \
      constraint_batch = 1; \
    else \
      constraint_batch = 0; \
  })
