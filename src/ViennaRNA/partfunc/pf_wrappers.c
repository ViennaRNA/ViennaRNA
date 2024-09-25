#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/problist.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/partfunc/local.h"

PUBLIC float
vrna_pf_fold(const char *seq,
             char       *structure,
             vrna_ep_t  **pl)
{
  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if (!pl) /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;

  vc  = vrna_fold_compound(seq, &md, 0);
  mfe = (double)vrna_mfe(vc, NULL);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if (pl)
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  vrna_fold_compound_free(vc);

  return free_energy;
}


PUBLIC float
vrna_pf_circfold(const char *seq,
                 char       *structure,
                 vrna_ep_t  **pl)
{
  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if (!pl) /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;

  vc  = vrna_fold_compound(seq, &md, 0);
  mfe = (double)vrna_mfe(vc, NULL);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if (pl)
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  vrna_fold_compound_free(vc);

  return free_energy;
}


PUBLIC float
vrna_pf_alifold(const char  **strings,
                char        *structure,
                vrna_ep_t   **pl)
{
  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if (!pl) /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;

  vc  = vrna_fold_compound_comparative(strings, &md, VRNA_OPTION_DEFAULT);
  mfe = (double)vrna_pf(vc, structure);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if (pl)
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  vrna_fold_compound_free(vc);

  return free_energy;
}


PUBLIC float
vrna_pf_circalifold(const char  **sequences,
                    char        *structure,
                    vrna_ep_t   **pl)
{
  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if (!pl) /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;

  vc  = vrna_fold_compound_comparative(sequences, &md, VRNA_OPTION_DEFAULT);
  mfe = (double)vrna_mfe(vc, structure);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if (pl)
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  vrna_fold_compound_free(vc);

  return free_energy;
}


PUBLIC int
vrna_pfl_fold_cb(const char                 *sequence,
                 int                        window_size,
                 int                        max_bp_span,
                 vrna_probs_window_f cb,
                 void                       *data)
{
  unsigned int          options;
  int                   r;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);       /* get default parameters */

  md.compute_bpp  = 1;            /* turn on base pair probability computations */
  md.window_size  = window_size;  /* set size of sliding window */
  md.max_bp_span  = max_bp_span;  /* set maximum base pair span */

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);

  options = VRNA_PROBS_WINDOW_BPP; /* always compute base pair probabilities */

  r = vrna_probs_window(vc, 0, options, cb, data);

  vrna_fold_compound_free(vc);

  return r;
}


PUBLIC int
vrna_pfl_fold_up_cb(const char                  *sequence,
                    int                         ulength,
                    int                         window_size,
                    int                         max_bp_span,
                    vrna_probs_window_f  cb,
                    void                        *data)
{
  unsigned int          options;
  int                   r;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);       /* get default parameters */

  md.compute_bpp  = 1;            /* turn on base pair probability computations */
  md.window_size  = window_size;  /* set size of sliding window */
  md.max_bp_span  = max_bp_span;  /* set maximum base pair span */

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);

  options = VRNA_PROBS_WINDOW_UP; /* compute unpaired probabilties */

  r = vrna_probs_window(vc, ulength, options, cb, data);

  vrna_fold_compound_free(vc);

  return r;
}
