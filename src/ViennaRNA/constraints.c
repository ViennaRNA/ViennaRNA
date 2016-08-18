/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/constraints.h"


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

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC  void
vrna_constraints_add( vrna_fold_compound_t *vc,
                      const char *constraint,
                      unsigned int options){

  int         i, d;

  if(vc){
    if(!vc->hc)
      vrna_hc_init(vc);

    if(options & VRNA_CONSTRAINT_DB){ /* apply hard constraints from dot-bracket notation */
      vrna_hc_add_from_db(vc, constraint, options);
    } else { /* constraints from file is the default */
      plist *p, *c = vrna_file_constraints_read(constraint, vc->length, 0);

      /* now do something with the constraints we've just read */
      if(c){
        FLT_OR_DBL  **sc_bp       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (vc->length + 1));
        FLT_OR_DBL  *sc_up        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
        int     sc_up_present = 0;
        int     sc_bp_present = 0;
        int     num_hc_up     = 0;
        int     mem_hc_up     = 10;
        vrna_hc_up_t  *hc_up  = NULL;

        hc_up = vrna_alloc(sizeof(vrna_hc_up_t) * mem_hc_up);

        for(i = 0; i <= vc->length; i++)
          sc_bp[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));

        for(p = c; p->i; p++){
          if(p->type & 4096){ /* soft constraint */

            if(num_hc_up > 0){ /* apply previously collected hc_up */
              hc_up[num_hc_up].position = 0; /* mark end of list */
              vrna_hc_add_up_batch(vc, hc_up);
              num_hc_up   = 0;
            }

            if(p->j == 0){  /* pseudo energy for unpairedness */
              sc_up_present = 1;
              sc_up[p->i] += (FLT_OR_DBL)p->p;
            } else { /* pseudo energy for base pair */
              sc_bp_present = 1;
              sc_bp[p->i][p->j] += (FLT_OR_DBL)p->p;
            }
          } else {  /* hard constraint */
            if(p->j == 0){ /* collect successive hard constraints for single unpaired nucleotides */
              hc_up[num_hc_up].position = p->i;
              hc_up[num_hc_up].options  = (char)p->type;
              num_hc_up++;
              if(num_hc_up == mem_hc_up){
                mem_hc_up = (int)(1.2 * mem_hc_up);
                hc_up  = (vrna_hc_up_t *)vrna_realloc(hc_up, sizeof(vrna_hc_up_t) * mem_hc_up);
              }
            } else {

              if(num_hc_up > 0){ /* apply previously collected hc_up */
                hc_up[num_hc_up].position = 0; /* mark end of list */
                vrna_hc_add_up_batch(vc, hc_up);
                num_hc_up   = 0;
              }

              if(p->i == p->j){ 
                d = 0;
                if(1024 & p->type)
                  d = -1;
                else if(2048 & p->type)
                  d = 1;
                if(p->type & 8192)
                  vrna_hc_add_bp_nonspecific(vc, p->i, d, (char)(p->type) | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
                else
                  vrna_hc_add_bp_nonspecific(vc, p->i, d, (char)(p->type));
              } else {
                if(p->type & 8192){
                  vrna_hc_add_bp(vc, p->i, p->j, (char)(p->type) | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
                } else {
                  vrna_hc_add_bp(vc, p->i, p->j, (char)(p->type));
                }
              }
            }
          }
        }

        /* ############################### */
        /* add hard constraints            */
        /* ############################### */
        if(num_hc_up > 0){
          hc_up[num_hc_up].position = 0; /* mark end of list */
          vrna_hc_add_up_batch(vc, hc_up);
        }
        free(hc_up);

        /* ############################### */
        /* init empty soft constraints     */
        /* ############################### */
        if(sc_up_present || sc_bp_present){
          vrna_sc_init(vc);
          if(sc_bp_present)
            vrna_sc_add_bp(vc, (const FLT_OR_DBL **)sc_bp, options);
          if(sc_up_present)
            vrna_sc_add_up(vc, (const FLT_OR_DBL *)sc_up, options);
        }
        /* clean up */
        for(i = 0; i <= vc->length; i++)
          free(sc_bp[i]);
        free(sc_bp);
        free(sc_up);
      }

      free(c);
    }
  }
}
