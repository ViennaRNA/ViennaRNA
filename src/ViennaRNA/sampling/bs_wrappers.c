/*
 *                Boltzmann Sampling wrappers
 *
 *                Ronny Lorenz
 *                ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sampling/basic.h"


/*
 #################################
 # PREPROCESSOR DEFININTIONS     #
 #################################
 */

struct structure_list {
  unsigned int  num;
  char          **list;
};


/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE void
store_sample_list(const char  *structure,
                  void        *data);


PRIVATE void
store_sample(const char *structure,
             void       *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_pbacktrack5(vrna_fold_compound_t *fc,
                 unsigned int         length)
{
  char          *structure = NULL;
  unsigned int  i;

  i = vrna_pbacktrack5_cb(fc,
                          1,
                          length,
                          &store_sample,
                          (void *)&structure,
                          VRNA_PBACKTRACK_DEFAULT);

  if (i)
    return structure;

  free(structure);

  return NULL;
}


PUBLIC char **
vrna_pbacktrack5_num(vrna_fold_compound_t *fc,
                     unsigned int         num_samples,
                     unsigned int         length,
                     unsigned int         options)
{
  unsigned int          i;
  struct structure_list data;

  data.num      = 0;
  data.list     = (char **)vrna_alloc(sizeof(char *) * num_samples);
  data.list[0]  = NULL;

  i = vrna_pbacktrack5_cb(fc,
                          num_samples,
                          length,
                          &store_sample_list,
                          (void *)&data,
                          options);

  if (i > 0) {
    /* re-allocate memory */
    data.list           = (char **)vrna_realloc(data.list, sizeof(char *) * (data.num + 1));
    data.list[data.num] = NULL;
  } else {
    free(data.list);
    return NULL;
  }

  return data.list;
}


PUBLIC unsigned int
vrna_pbacktrack5_cb(vrna_fold_compound_t              *fc,
                    unsigned int                      num_samples,
                    unsigned int                      length,
                    vrna_bs_result_f  bs_cb,
                    void                              *data,
                    unsigned int                      options)
{
  unsigned int          i;
  vrna_pbacktrack_mem_t nr_mem = NULL;

  i = vrna_pbacktrack5_resume_cb(fc,
                                 num_samples,
                                 length,
                                 bs_cb,
                                 data,
                                 &nr_mem,
                                 options);

  vrna_pbacktrack_mem_free(nr_mem);

  return i;
}


PUBLIC char **
vrna_pbacktrack5_resume(vrna_fold_compound_t  *vc,
                        unsigned int          num_samples,
                        unsigned int          length,
                        vrna_pbacktrack_mem_t *nr_mem,
                        unsigned int          options)
{
  unsigned int          i;
  struct structure_list data;

  if (vc) {
    data.num      = 0;
    data.list     = (char **)vrna_alloc(sizeof(char *) * num_samples);
    data.list[0]  = NULL;

    i = vrna_pbacktrack5_resume_cb(vc,
                                   num_samples,
                                   length,
                                   &store_sample_list,
                                   (void *)&data,
                                   nr_mem,
                                   options);

    if (i > 0) {
      /* re-allocate memory */
      data.list           = (char **)vrna_realloc(data.list, sizeof(char *) * (data.num + 1));
      data.list[data.num] = NULL;
    } else {
      free(data.list);
      return NULL;
    }

    return data.list;
  }

  return NULL;
}


PUBLIC unsigned int
vrna_pbacktrack5_resume_cb(vrna_fold_compound_t             *fc,
                           unsigned int                     num_samples,
                           unsigned int                     end,
                           vrna_bs_result_f bs_cb,
                           void                             *data,
                           vrna_pbacktrack_mem_t            *nr_mem,
                           unsigned int                     options)
{
  return vrna_pbacktrack_sub_resume_cb(fc,
                                       num_samples,
                                       1,
                                       end,
                                       bs_cb,
                                       data,
                                       nr_mem,
                                       options);
}


/*
 * stochastic backtracking in pf_fold arrays
 * returns random structure S with Boltzman probabilty
 * p(S) = exp(-E(S)/kT)/Z
 */
PUBLIC char *
vrna_pbacktrack(vrna_fold_compound_t *fc)
{
  if (fc)
    return vrna_pbacktrack5(fc, fc->length);

  return NULL;
}


PUBLIC char **
vrna_pbacktrack_num(vrna_fold_compound_t  *fc,
                    unsigned int          num_samples,
                    unsigned int          options)
{
  if (fc) {
    return vrna_pbacktrack5_num(fc,
                                num_samples,
                                fc->length,
                                options);
  }

  return NULL;
}


PUBLIC unsigned int
vrna_pbacktrack_cb(vrna_fold_compound_t             *fc,
                   unsigned int                     num_samples,
                   vrna_bs_result_f bs_cb,
                   void                             *data,
                   unsigned int                     options)
{
  if (fc) {
    return vrna_pbacktrack5_cb(fc,
                               num_samples,
                               fc->length,
                               bs_cb,
                               data,
                               options);
  }

  return 0;
}


PUBLIC char **
vrna_pbacktrack_resume(vrna_fold_compound_t   *fc,
                       unsigned int           num_samples,
                       vrna_pbacktrack_mem_t  *nr_mem,
                       unsigned int           options)
{
  if (fc) {
    return vrna_pbacktrack5_resume(fc,
                                   num_samples,
                                   fc->length,
                                   nr_mem,
                                   options);
  }

  return NULL;
}


PUBLIC unsigned int
vrna_pbacktrack_resume_cb(vrna_fold_compound_t              *fc,
                          unsigned int                      num_samples,
                          vrna_bs_result_f  bs_cb,
                          void                              *data,
                          vrna_pbacktrack_mem_t             *nr_mem,
                          unsigned int                      options)
{
  if (fc) {
    return vrna_pbacktrack5_resume_cb(fc,
                                      num_samples,
                                      fc->length,
                                      bs_cb,
                                      data,
                                      nr_mem,
                                      options);
  }

  return 0;
}


PUBLIC char *
vrna_pbacktrack_sub(vrna_fold_compound_t  *fc,
                    unsigned int          start,
                    unsigned int          end)
{
  char          *structure = NULL;
  unsigned int  i;

  i = vrna_pbacktrack_sub_cb(fc,
                             1,
                             start,
                             end,
                             &store_sample,
                             (void *)&structure,
                             VRNA_PBACKTRACK_DEFAULT);

  if (i)
    return structure;

  free(structure);

  return NULL;
}


PUBLIC char **
vrna_pbacktrack_sub_num(vrna_fold_compound_t  *fc,
                        unsigned int          num_samples,
                        unsigned int          start,
                        unsigned int          end,
                        unsigned int          options)
{
  unsigned int          i;
  struct structure_list data;

  data.num      = 0;
  data.list     = (char **)vrna_alloc(sizeof(char *) * num_samples);
  data.list[0]  = NULL;

  i = vrna_pbacktrack_sub_cb(fc,
                             num_samples,
                             start,
                             end,
                             &store_sample_list,
                             (void *)&data,
                             options);

  if (i > 0) {
    /* re-allocate memory */
    data.list           = (char **)vrna_realloc(data.list, sizeof(char *) * (data.num + 1));
    data.list[data.num] = NULL;
  } else {
    free(data.list);
    return NULL;
  }

  return data.list;
}


PUBLIC unsigned int
vrna_pbacktrack_sub_cb(vrna_fold_compound_t             *fc,
                       unsigned int                     num_samples,
                       unsigned int                     start,
                       unsigned int                     end,
                       vrna_bs_result_f bs_cb,
                       void                             *data,
                       unsigned int                     options)
{
  unsigned int          i;
  vrna_pbacktrack_mem_t nr_mem = NULL;

  i = vrna_pbacktrack_sub_resume_cb(fc,
                                    num_samples,
                                    start,
                                    end,
                                    bs_cb,
                                    data,
                                    &nr_mem,
                                    options);

  vrna_pbacktrack_mem_free(nr_mem);

  return i;
}


PUBLIC char **
vrna_pbacktrack_sub_resume(vrna_fold_compound_t   *vc,
                           unsigned int           num_samples,
                           unsigned int           start,
                           unsigned int           end,
                           vrna_pbacktrack_mem_t  *nr_mem,
                           unsigned int           options)
{
  unsigned int          i;
  struct structure_list data;

  if (vc) {
    data.num      = 0;
    data.list     = (char **)vrna_alloc(sizeof(char *) * num_samples);
    data.list[0]  = NULL;

    i = vrna_pbacktrack_sub_resume_cb(vc,
                                      num_samples,
                                      start,
                                      end,
                                      &store_sample_list,
                                      (void *)&data,
                                      nr_mem,
                                      options);

    if (i > 0) {
      /* re-allocate memory */
      data.list           = (char **)vrna_realloc(data.list, sizeof(char *) * (data.num + 1));
      data.list[data.num] = NULL;
    } else {
      free(data.list);
      return NULL;
    }

    return data.list;
  }

  return NULL;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
store_sample_list(const char  *structure,
                  void        *data)
{
  struct structure_list *d = (struct structure_list *)data;

  if (structure)
    d->list[d->num++] = strdup(structure);
  else
    d->list[d->num++] = NULL;
}


PRIVATE void
store_sample(const char *structure,
             void       *data)
{
  char **s = (char **)data;

  if (structure)
    *s = strdup(structure);
  else
    *s = NULL;
}
