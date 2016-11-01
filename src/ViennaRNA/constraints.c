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
#include "ViennaRNA/commands.h"
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

  if(vc){
    if(!vc->hc)
      vrna_hc_init(vc);

    if(options & VRNA_CONSTRAINT_DB){ /* apply hard constraints from dot-bracket notation */
      vrna_hc_add_from_db(vc, constraint, options);
    } else { /* constraints from file is the default */
      vrna_file_commands_apply(vc, constraint, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);
    }
  }
}
