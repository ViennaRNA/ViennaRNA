%module RNA
//%pragma(perl5)  modulecode="@EXPORT=qw(fold);"
%pragma(perl5)  include="RNA.pod"

%{

extern "C" {
#include  <ViennaRNA/model.h>
#include  <ViennaRNA/data_structures.h>
#include  <ViennaRNA/dp_matrices.h>
#include  <ViennaRNA/alphabet.h>
#include  <ViennaRNA/grammar.h>
#include  <ViennaRNA/unstructured_domains.h>
#include  <ViennaRNA/structured_domains.h>
#include  <ViennaRNA/commands.h>

#include  <ViennaRNA/utils.h>
#include  <ViennaRNA/structure_utils.h>
#include  <ViennaRNA/string_utils.h>
#include  <ViennaRNA/fold_vars.h>

#include  <ViennaRNA/constraints.h>
#include  <ViennaRNA/constraints_hard.h>
#include  <ViennaRNA/constraints_soft.h>
#include  <ViennaRNA/constraints_SHAPE.h>
#include  <ViennaRNA/constraints_ligand.h>

#undef fold
#include  <ViennaRNA/mfe.h>
#include  <ViennaRNA/fold.h>
#include  <ViennaRNA/eval.h>
#include  <ViennaRNA/cofold.h>
#include  <ViennaRNA/part_func.h>
#include  <ViennaRNA/equilibrium_probs.h>
#include  <ViennaRNA/boltzmann_sampling.h>
#include  <ViennaRNA/part_func_co.h>
#include  <ViennaRNA/naview.h>
#include  <ViennaRNA/plot_layouts.h>
#include  <ViennaRNA/plot_structure.h>
#include  <ViennaRNA/plot_aln.h>
#include  <ViennaRNA/PS_dot.h>
#include  <ViennaRNA/inverse.h>
#include  <ViennaRNA/RNAstruct.h>
#include  <ViennaRNA/treedist.h>
#include  <ViennaRNA/stringdist.h>
#include  <ViennaRNA/profiledist.h>
#include  <ViennaRNA/dist_vars.h>
#include  <ViennaRNA/pair_mat.h>
#include  <ViennaRNA/subopt.h>
#include  <ViennaRNA/energy_const.h>
#include  <ViennaRNA/params.h>
#include  <ViennaRNA/duplex.h>
#include  <ViennaRNA/alifold.h>
#include  <ViennaRNA/aln_util.h>
#include  <ViennaRNA/findpath.h>
#include  <ViennaRNA/Lfold.h>
#include  <ViennaRNA/read_epars.h>
#include  <ViennaRNA/move_set.h>
#include  <ViennaRNA/file_formats.h>
#include  <ViennaRNA/file_formats_msa.h>

#include  <ViennaRNA/exterior_loops.h>
#include  <ViennaRNA/hairpin_loops.h>
#include  <ViennaRNA/interior_loops.h>
#include  <ViennaRNA/multibranch_loops.h>

#include  <ViennaRNA/combinatorics.h>
}

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

/* %constant double VERSION = 0.3; */
%include version.i
%include typemaps.i

// Typemaps that are independent of scripting language

// This cleans up the char ** array after the function call
%typemap(freearg) char ** {
         free($1);
}

// Additional target language specific typemaps
%include tmaps.i

/* handle exceptions */
%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

/* prepare conversions to native types, such as lists */
%include "std_pair.i";
%include "std_vector.i";
%include "std_string.i";

namespace std {
  %template(DoublePair) std::pair<double,double>;
  %template(IntVector) std::vector<int>;
  %template(UIntVector) std::vector<unsigned int>;
  %template(DoubleVector) std::vector<double>;
  %template(StringVector) std::vector<string>;
  %template(ConstCharVector) std::vector<const char*>;
  %template(SOLUTIONVector) std::vector<SOLUTION>;
  %template(CoordinateVector) std::vector<COORDINATE>;
  %template(DoubleDoubleVector) std::vector< std::vector<double> > ;
  %template(IntIntVector) std::vector<std::vector<int> > ;
};

%{
#include <string>
#include <cstring>

  const char *convert_vecstring2veccharcp(const std::string & s){
    return s.c_str();
  }

  char *convert_vecstring2veccharp(const std::string & s){
    char *pc = new char[s.size()+1];
    std::strcpy(pc, s.c_str());
    return pc;
  }
  
  short convert_vecint2vecshort(const int & i){
    return (short) i;
  }

  FLT_OR_DBL convert_vecdbl2vecFLR_OR_DBL(const double & d){
    return (FLT_OR_DBL) d;
  }

%}

//%title "Interface to the Vienna RNA library"

/* do not wrap any function prefixed by 'vrna_' */
%rename("$ignore",  %$isfunction, regextarget=1) "^vrna_";

/* do not wrap any data structure, typedef, enum, or constant prefixed by 'vrna_' || 'VRNA_' */
%rename("$ignore",  %$isclass, regextarget=1) "^vrna_";
%rename("$ignore",  %$istypedef, regextarget=1) "^vrna_";
%rename("$ignore",  %$isenum, regextarget=1) "^vrna_";
%rename("$ignore",  %$isconstant, regextarget=1) "^VRNA_";
%rename("$ignore",  %$isconstant, regextarget=1) "^vrna_";

/*############################################*/
/* Include all relevant interface definitions */
/*############################################*/
%include params.i
%include model_details.i
%include utils.i
%include plotting.i
%include constraints.i
%include constraints_hard.i
%include constraints_soft.i
%include constraints_SHAPE.i
%include constraints_ligand.i
%include eval.i
%include mfe.i
%include part_func.i
%include subopt.i
%include inverse.i
%include compare.i
%include file_formats.i
%include grammar.i
%include commands.i
%include combinatorics.i
%include fold_compound.i


/**********************************************/
/* BEGIN interface for data structures        */
/**********************************************/

%ignore folden;
%ignore node;
%ignore snoopT;
%ignore dupVar;

%ignore plist;
%ignore cpair;


%include <ViennaRNA/data_structures.h>

%include <ViennaRNA/dp_matrices.h>

//%subsection "Global Variables to Modify Folding"
//extern double *pr;  /*  base pairing prob. matrix */

%include  <ViennaRNA/fold_vars.h>
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

%include <ViennaRNA/duplex.h>
%include <ViennaRNA/Lfold.h>

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

%include <ViennaRNA/findpath.h>

%include  <ViennaRNA/move_set.h>

%ignore move_gradient;
%ignore move_first;
%ignore move_adaptive;
%ignore browse_neighs_pt;
%ignore browse_neighs;
%ignore print_stren;
%ignore print_str;
%ignore copy_arr;
%ignore allocopy;

