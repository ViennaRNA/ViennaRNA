#ifndef VRNA_PROBING_DATA_HELPERS
#define VRNA_PROBING_DATA_HELPERS

#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_compound.h"

typedef struct {
  unsigned int        count;
  vrna_array(char *)  files;
  vrna_array(char *)  strategies;
  vrna_array(char *)  preprocessing;
  vrna_array(char *)  prior_unpaired;
  vrna_array(char *)  prior_paired;
} probing_data_t;


#define ggo_get_probing_data(argc, \
                             argv, \
                             ggostruct, \
                             probing_data_p)  { \
    probing_data_p = NULL; \
    if (ggostruct.sp_data_given) \
      probing_data_p = extract_probing_options(argc, argv, ggostruct.sp_data_given); \
    /* backward compatibility/convenience wrapper */ \
    if (ggostruct.shape_given) { \
      if (probing_data_p == NULL) { \
        probing_data_p = (probing_data_t *)vrna_alloc(sizeof(probing_data_t)); \
        probing_data_p->count = 0; \
        vrna_array_init(probing_data_p->files); \
        vrna_array_init(probing_data_p->strategies); \
        vrna_array_init(probing_data_p->preprocessing); \
      } \
      probing_data_p->count++; \
      vrna_array_append(probing_data_p->files, strdup(ggostruct.shape_arg)); \
      vrna_array_append(probing_data_p->strategies, strdup(ggostruct.shapeMethod_arg)); \
      vrna_array_append(probing_data_p->preprocessing, strdup(ggostruct.shapeConversion_arg)); \
    } \
  }


/*
 * HELPER MACRO FOR SHAPE REACTIVITY DATA INCORPORATION
 */

#define ggo_get_SHAPE(ggostruct, \
                      SHAPE_switch, \
                      SHAPE_files, \
                      SHAPE_method, \
                      SHAPE_conversion)  { \
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
  }


void
probing_data_free(probing_data_t *dat);


int
apply_probing_data(vrna_fold_compound_t *fc,
                   probing_data_t       *d);


void
vrna_constraints_add_SHAPE(vrna_fold_compound_t *fc,
                           const char           *shape_file,
                           const char           *shape_method,
                           const char           *shape_conversion,
                           int                  verbose,
                           unsigned int         constraint_type);


void
vrna_constraints_add_SHAPE_ali(vrna_fold_compound_t *fc,
                               const char           *shape_method,
                               const char           **shape_files,
                               const int            *shape_file_association,
                               int                  verbose,
                               unsigned int         constraint_type);

/**
 *  @brief  Parse a character string and extract the encoded SHAPE reactivity conversion
 *          method and possibly the parameters for conversion into pseudo free energies
 *
 *  @ingroup soft_cosntraints
 *
 *  @param  method_string   The string that contains the encoded SHAPE reactivity conversion method
 *  @param  method          A pointer to the memory location where the method character will be stored
 *  @param  param_1         A pointer to the memory location where the first parameter of the corresponding method will be stored
 *  @param  param_2         A pointer to the memory location where the second parameter of the corresponding method will be stored
 *  @return                 1 on successful extraction of the method, 0 on errors
 */
int
vrna_sc_SHAPE_parse_method(const char *method_string,
                           char       *method,
                           float      *param_1,
                           float      *param_2);


double **
vrna_probing_data_load_n_distribute(unsigned int  n_seq,
                                    unsigned int  *ns,
                                    const char    **sequences,
                                    const char    **file_names,
                                    const int     *file_name_association,
                                    unsigned int  options);


probing_data_t *
extract_probing_options(int     argc,
                        char    *argv[],
                        size_t  expected_num_data);

#endif
