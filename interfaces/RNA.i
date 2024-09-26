#ifdef SWIGPYTHON
%module(moduleimport="from . import _RNA") RNA
#else
%module RNA
#endif

//%pragma(perl5)  modulecode="@EXPORT=qw(fold);"
%pragma(perl5)  include="RNA.pod"

// ignore SWIG Warning 312 from nested unions
#pragma SWIG nowarn=312


%{

extern "C" {
#include  <ViennaRNA/model.h>
#include  <ViennaRNA/datastructures/basic.h>
#include  <ViennaRNA/datastructures/array.h>
#include  <ViennaRNA/datastructures/sparse_mx.h>
#include  <ViennaRNA/datastructures/dp_matrices.h>
#include  <ViennaRNA/fold_compound.h>
#include  <ViennaRNA/unstructured_domains.h>
#include  <ViennaRNA/structured_domains.h>

#include  <ViennaRNA/grammar/basic.h>
#include  <ViennaRNA/grammar/mfe.h>
#include  <ViennaRNA/grammar/partfunc.h>

#include  <ViennaRNA/utils/basic.h>
#include  <ViennaRNA/utils/strings.h>
#include  <ViennaRNA/utils/log.h>
#include  <ViennaRNA/fold_vars.h>

#include <ViennaRNA/sequences/alphabet.h>
#include <ViennaRNA/sequences/sequence.h>
#include <ViennaRNA/sequences/alignments.h>
#include <ViennaRNA/sequences/utils.h>

#include <ViennaRNA/structures/benchmark.h>
#include <ViennaRNA/structures/dotbracket.h>
#include <ViennaRNA/structures/helix.h>
#include <ViennaRNA/structures/metrics.h>
#include <ViennaRNA/structures/pairtable.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/structures/shapes.h>
#include <ViennaRNA/structures/tree.h>
#include <ViennaRNA/structures/utils.h>
#include <ViennaRNA/structures/centroid.h>
#include <ViennaRNA/structures/mea.h>

#include  <ViennaRNA/params/constants.h>
#include  <ViennaRNA/params/basic.h>
#include  <ViennaRNA/params/io.h>
#include  <ViennaRNA/params/salt.h>
#include  <ViennaRNA/params/default.h>

#include  <ViennaRNA/constraints/basic.h>
#include  <ViennaRNA/constraints/hard.h>
#include  <ViennaRNA/constraints/soft.h>
#include  <ViennaRNA/constraints/soft_special.h>
#include  <ViennaRNA/constraints/ligand.h>

#include  <ViennaRNA/probing/basic.h>
#include  <ViennaRNA/probing/SHAPE.h>


#ifdef VRNA_WITH_NAVIEW_LAYOUT
#include  <ViennaRNA/plotting/naview/naview.h>
#endif
#include  <ViennaRNA/plotting/layouts.h>
#include  <ViennaRNA/plotting/structures.h>
#include  <ViennaRNA/plotting/alignments.h>
#include  <ViennaRNA/plotting/probabilities.h>

#include  <ViennaRNA/io/commands.h>
#include  <ViennaRNA/io/file_formats.h>
#include  <ViennaRNA/io/file_formats_msa.h>
#include  <ViennaRNA/io/utils.h>

#include  <ViennaRNA/eval/basic.h>
#include  <ViennaRNA/eval/structures.h>
#include  <ViennaRNA/eval/exterior.h>
#include  <ViennaRNA/eval/hairpin.h>
#include  <ViennaRNA/eval/internal.h>
#include  <ViennaRNA/eval/multibranch.h>
#include  <ViennaRNA/eval/gquad.h>

#include  <ViennaRNA/mfe/global.h>
#include  <ViennaRNA/mfe/local.h>
#include  <ViennaRNA/mfe/exterior.h>
#include  <ViennaRNA/mfe/internal.h>
#include  <ViennaRNA/mfe/multibranch.h>
#include  <ViennaRNA/mfe/gquad.h>

#include  <ViennaRNA/backtrack/global.h>
#include  <ViennaRNA/backtrack/exterior.h>
#include  <ViennaRNA/backtrack/internal.h>
#include  <ViennaRNA/backtrack/multibranch.h>
#include  <ViennaRNA/backtrack/gquad.h>

#include  <ViennaRNA/partfunc/global.h>
#include  <ViennaRNA/partfunc/local.h>
#include  <ViennaRNA/partfunc/exterior.h>
#include  <ViennaRNA/partfunc/internal.h>
#include  <ViennaRNA/partfunc/multibranch.h>
#include  <ViennaRNA/partfunc/gquad.h>

#include  <ViennaRNA/subopt/basic.h>
#include  <ViennaRNA/subopt/wuchty.h>
#include  <ViennaRNA/subopt/zuker.h>

#include  <ViennaRNA/sampling/basic.h>

#include  <ViennaRNA/probabilities/basepairs.h>
#include  <ViennaRNA/probabilities/structures.h>

#include  <ViennaRNA/fold.h>
#include  <ViennaRNA/cofold.h>
#include  <ViennaRNA/alifold.h>

#include  <ViennaRNA/part_func_co.h>
#include  <ViennaRNA/concentrations.h>
#include  <ViennaRNA/LPfold.h>
#include  <ViennaRNA/heat_capacity.h>

#ifdef VRNA_WITH_SVM
#include  <ViennaRNA/zscore/basic.h>
#endif

#include  <ViennaRNA/inverse/basic.h>

#include  <ViennaRNA/RNAstruct.h>
#include  <ViennaRNA/treedist.h>
#include  <ViennaRNA/stringdist.h>
#include  <ViennaRNA/profiledist.h>
#include  <ViennaRNA/dist_vars.h>
#include  <ViennaRNA/pair_mat.h>
#include  <ViennaRNA/duplex.h>

#include  <ViennaRNA/combinatorics/basic.h>

#include  <ViennaRNA/move_set.h>
#include  <ViennaRNA/landscape/paths.h>
#include  <ViennaRNA/landscape/findpath.h>
#include  <ViennaRNA/landscape/move.h>
#include  <ViennaRNA/landscape/neighbor.h>
#include  <ViennaRNA/landscape/walk.h>

#include  <ViennaRNA/mm.h>

#include  <ViennaRNA/static/energy_parameter_sets.h>
}

%}
//
%include cstring.i

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

%include version.i
%include typemaps.i

/* handle exceptions */
%include "exception.i"
%include "std_except.i"

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
  %template(UIntUIntVector) std::vector<std::vector<unsigned int> >;
  %template(DoubleVector) std::vector<double>;
  %template(StringVector) std::vector<string>;
  %template(ConstCharVector) std::vector<const char*>;
  %template(SOLUTIONVector) std::vector<SOLUTION>;
  %template(CoordinateVector) std::vector<COORDINATE>;
  %template(DoubleDoubleVector) std::vector< std::vector<double> > ;
  %template(IntIntVector) std::vector<std::vector<int> > ;
  %template(ElemProbVector) vector<vrna_ep_t>;
  %template(HelixVector) vector<vrna_hx_t>;
  %template(PathVector) std::vector<vrna_path_t>;
  %template(MoveVector) std::vector<vrna_move_t>;
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


%init %{
  /* always initialize the random number generator when loading the module */
  vrna_init_rand();
%}


// Typemaps that are independent of scripting language

// This cleans up the char ** array after the function call
%typemap(freearg) char ** {
         free($1);
}

%typemap(newfree) char * "free($1);";


// Additional target language specific typemaps
%include "tmaps.i"

/*############################################*/
/* Include all relevant interface definitions */
/*############################################*/
#ifdef SWIGPYTHON
#ifdef DOXY2SWIG
%include documentation.i
#endif
#endif

%include var_arrays.i
%include type_checks.i
%include data_structures.i
%include params.i
%include model_details.i
%include utils.i
%include plotting.i
%include constraints.i
%include constraints_hard.i
%include constraints_soft.i
%include constraints_ligand.i
%include constraints_mod.i
%include probing.i
%include probing_SHAPE.i
%include eval.i
%include loops.i
%include basic_algorithms.i
%include mfe.i
%include mfe_window.i
%include backtrack.i
%include part_func.i
%include boltzmann_sampling.i
%include pf_window.i
%include subopt.i
%include inverse.i
%include compare.i
%include file_formats.i
%include sequence.i
%include grammar.i
%include commands.i
%include combinatorics.i
%include duplex.i
%include move.i
%include neighbor.i
%include walk.i
%include paths.i
%include heat_capacity.i
%include fold_compound.i
%include dp_matrices.i
%include parameter_sets.i

/**********************************************/
/* BEGIN interface for data structures        */
/**********************************************/

//%subsection "Global Variables to Modify Folding"
//extern double *pr;  /*  base pairing prob. matrix */
%immutable;
extern bondT      *base_pair;
extern FLT_OR_DBL *pr;
extern int        *iindx;
%mutable;

%include  <ViennaRNA/fold_vars.h>

struct bondT {
  unsigned int  i;
  unsigned int  j;
};

%extend bondT {

  bondT *get(int i) {

    return self+i;
  }
}

