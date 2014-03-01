/* Functions for Loop energy computations */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/loop_energies.h"



/* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */
