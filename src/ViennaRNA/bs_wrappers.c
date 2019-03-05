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
#include "ViennaRNA/boltzmann_sampling.h"


/*
 #################################
 # PREPROCESSOR DEFININTIONS     #
 #################################
 */

struct nr_structure_list {
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
save_nr_samples(const char  *structure,
                void        *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

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
vrna_pbacktrack_nr(vrna_fold_compound_t *fc,
                   unsigned int         num_samples)
{
  char                    **structures;
  struct vrna_nr_memory_s *nr_mem;

  nr_mem      = NULL;
  structures  = vrna_pbacktrack_nr_resume(fc, num_samples, &nr_mem);

  vrna_pbacktrack_nr_free(nr_mem);

  return structures;
}


PUBLIC char **
vrna_pbacktrack_nr_resume(vrna_fold_compound_t  *vc,
                          unsigned int          num_samples,
                          vrna_nr_memory_t      *nr_mem)
{
  unsigned int              i;
  struct nr_structure_list  data;

  data.num      = 0;
  data.list     = (char **)vrna_alloc(sizeof(char *) * num_samples);
  data.list[0]  = NULL;

  i = vrna_pbacktrack_nr_resume_cb(vc,
                                   num_samples,
                                   &save_nr_samples,
                                   (void *)&data,
                                   nr_mem);

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
vrna_pbacktrack_nr_cb(vrna_fold_compound_t              *vc,
                      unsigned int                      num_samples,
                      vrna_boltzmann_sampling_callback  *bs_cb,
                      void                              *data)
{
  unsigned int            i;
  struct vrna_nr_memory_s *nr_mem;

  nr_mem  = NULL;
  i       = vrna_pbacktrack_nr_resume_cb(vc, num_samples, bs_cb, data, &nr_mem);

  vrna_pbacktrack_nr_free(nr_mem);

  return i;
}


PRIVATE void
save_nr_samples(const char  *structure,
                void        *data)
{
  struct nr_structure_list *d = (struct nr_structure_list *)data;

  if (structure)
    d->list[d->num++] = strdup(structure);
  else
    d->list[d->num++] = NULL;
}
