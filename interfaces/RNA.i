%module RNA
//%pragma(perl5)  modulecode="@EXPORT=qw(fold);"
%pragma(perl5)  include="RNA.pod"

%{
#include  "../src/ViennaRNA/data_structures.h"
#include  "../src/ViennaRNA/dp_matrices.h"
#include  "../src/ViennaRNA/model.h"
#include  "../src/ViennaRNA/utils.h"
#include  "../src/ViennaRNA/structure_utils.h"
#include  "../src/ViennaRNA/string_utils.h"
#include  "../src/ViennaRNA/fold_vars.h"
#include  "../src/ViennaRNA/constraints.h"
#undef fold
#include  "../src/ViennaRNA/mfe.h"
#include  "../src/ViennaRNA/fold.h"
#include  "../src/ViennaRNA/eval.h"
#include  "../src/ViennaRNA/cofold.h"
#include  "../src/ViennaRNA/part_func.h"
#include  "../src/ViennaRNA/part_func_co.h"
#include  "../src/ViennaRNA/naview.h"
#include  "../src/ViennaRNA/plot_layouts.h"
#include  "../src/ViennaRNA/plot_structure.h"
#include  "../src/ViennaRNA/plot_aln.h"
#include  "../src/ViennaRNA/PS_dot.h"
#include  "../src/ViennaRNA/inverse.h"
#include  "../src/ViennaRNA/RNAstruct.h"
#include  "../src/ViennaRNA/treedist.h"
#include  "../src/ViennaRNA/stringdist.h"
#include  "../src/ViennaRNA/profiledist.h"
#include  "../src/ViennaRNA/dist_vars.h"
#include  "../src/ViennaRNA/pair_mat.h"
#include  "../src/ViennaRNA/subopt.h"
#include  "../src/ViennaRNA/energy_const.h"
#include  "../src/ViennaRNA/params.h"
#include  "../src/ViennaRNA/duplex.h"
#include  "../src/ViennaRNA/alifold.h"
#include  "../src/ViennaRNA/aln_util.h"
#include  "../src/ViennaRNA/findpath.h"
#include  "../src/ViennaRNA/Lfold.h"
#include  "../src/ViennaRNA/read_epars.h"
#include  "../src/ViennaRNA/move_set.h"
#include  "../src/ViennaRNA/ligand.h"

%}
//
%include carrays.i
%array_functions(int, intP);
%array_class(int, intArray);
%array_functions(float, floatP);
%array_class(float, floatArray);
%array_functions(double, doubleP);
%array_class(double, doubleArray);
%array_functions(unsigned short, ushortP);
%array_functions(short, shortP);
%include cdata.i

#ifdef LARGE_PF
#undef FLT_OR_DBL
#define FLT_OR_DBL  double
#else
#undef FLT_OR_DBL
#define FLT_OR_DBL  float
#endif

%constant double VERSION = 0.3;
%include typemaps.i
%include tmaps.i  // additional typemaps

//%title "Interface to the Vienna RNA library"

/* do not wrap any function prefixed by 'vrna_' */
%rename("$ignore",  %$isfunction, regextarget=1) "^vrna_";

/* do not wrap any data structure, typedef, or enum prefixed by 'vrna_' */
%rename("$ignore",  %$isclass, regextarget=1) "^vrna_";
%rename("$ignore",  %$istypedef, regextarget=1) "^vrna_";
%rename("$ignore",  %$isenum, regextarget=1) "^vrna_";

/*############################################*/
/* Include all relevant interface definitions */
/*############################################*/
%include params.i
%include model_details.i
%include fold_compound.i
%include utils.i
%include plotting.i
%include constraints.i
%include eval.i
%include mfe.i
%include part_func.i
%include subopt.i
%include inverse.i
%include compare.i

/**********************************************/
/* BEGIN interface for data structures        */
/**********************************************/

%ignore folden;
%ignore node;
%ignore snoopT;
%ignore dupVar;

%ignore plist;
%ignore cpair;


%include "../src/ViennaRNA/data_structures.h"

%include "../src/ViennaRNA/dp_matrices.h"

//%subsection "Global Variables to Modify Folding"
//extern double *pr;  /*  base pairing prob. matrix */

%include  "../src/ViennaRNA/fold_vars.h"
%extend bondT {
	bondT *get(int i) {
	   return self+i;
	}
}


typedef struct {
  int i;
  int j;
  char *structure;
  float energy;
} duplexT;

// problem: wrapper will not free the the structure hidden in duplexT
// even more problem for duplex_subopt()

%include "../src/ViennaRNA/duplex.h"
%include "../src/ViennaRNA/Lfold.h"

/**********************************************/
/* BEGIN interface for findpath heursitic     */
/**********************************************/

/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(path) vrna_path_t;

/* no default constructor / destructor */
%nodefaultdtor vrna_path_t;

typedef struct {
  double en;  /**<  @brief  Free energy of current structure */
  char *s;    /**<  @brief  Secondary structure in dot-bracket notation */
} vrna_path_t;

%extend vrna_path_t {
  vrna_path_t *get(int i) {
    return $self+i;
  }

  int size() {
    path_t *st;
    for (st=$self; st->s; st++);
    return (int)(st-$self);
  }

  ~vrna_path_t() {
    vrna_path_t *st;
    for (st=$self; st->s; st++)
      free(st->s);
    free($self);
  }
}

%newobject get_path;

%include "../src/ViennaRNA/findpath.h"

%include  "../src/ViennaRNA/move_set.h"

%ignore move_gradient;
%ignore move_first;
%ignore move_adaptive;
%ignore browse_neighs_pt;
%ignore browse_neighs;
%ignore print_stren;
%ignore print_str;
%ignore copy_arr;
%ignore allocopy;

