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
#include  <ViennaRNA/fold_compound.h>
#include  <ViennaRNA/dp_matrices.h>
#include  <ViennaRNA/alphabet.h>
#include  <ViennaRNA/sequence.h>
#include  <ViennaRNA/grammar.h>
#include  <ViennaRNA/unstructured_domains.h>
#include  <ViennaRNA/structured_domains.h>
#include  <ViennaRNA/commands.h>

#include  <ViennaRNA/utils/basic.h>
#include  <ViennaRNA/utils/structures.h>
#include  <ViennaRNA/utils/strings.h>
#include  <ViennaRNA/utils/alignments.h>
#include  <ViennaRNA/fold_vars.h>

#include  <ViennaRNA/params/constants.h>
#include  <ViennaRNA/params/basic.h>
#include  <ViennaRNA/params/io.h>

#include  <ViennaRNA/constraints/basic.h>
#include  <ViennaRNA/constraints/hard.h>
#include  <ViennaRNA/constraints/soft.h>
#include  <ViennaRNA/constraints/SHAPE.h>
#include  <ViennaRNA/constraints/ligand.h>

#include  <ViennaRNA/plotting/naview.h>
#include  <ViennaRNA/plotting/layouts.h>
#include  <ViennaRNA/plotting/structures.h>
#include  <ViennaRNA/plotting/alignments.h>
#include  <ViennaRNA/plotting/probabilities.h>

#include  <ViennaRNA/io/file_formats.h>
#include  <ViennaRNA/io/file_formats_msa.h>
#include  <ViennaRNA/io/utils.h>

#include  <ViennaRNA/loops/external.h>
#include  <ViennaRNA/loops/hairpin.h>
#include  <ViennaRNA/loops/internal.h>
#include  <ViennaRNA/loops/multibranch.h>

#include  <ViennaRNA/mfe.h>
#include  <ViennaRNA/mfe_window.h>
#include  <ViennaRNA/fold.h>
#include  <ViennaRNA/eval.h>
#include  <ViennaRNA/cofold.h>
#include  <ViennaRNA/alifold.h>

#include  <ViennaRNA/part_func.h>
#include  <ViennaRNA/part_func_window.h>
#include  <ViennaRNA/part_func_co.h>
#include  <ViennaRNA/equilibrium_probs.h>
#include  <ViennaRNA/boltzmann_sampling.h>
#include  <ViennaRNA/concentrations.h>
#include  <ViennaRNA/LPfold.h>
#include  <ViennaRNA/centroid.h>
#include  <ViennaRNA/MEA.h>

#include  <ViennaRNA/inverse.h>
#include  <ViennaRNA/RNAstruct.h>
#include  <ViennaRNA/treedist.h>
#include  <ViennaRNA/stringdist.h>
#include  <ViennaRNA/profiledist.h>
#include  <ViennaRNA/dist_vars.h>
#include  <ViennaRNA/pair_mat.h>
#include  <ViennaRNA/subopt.h>
#include  <ViennaRNA/duplex.h>

#include  <ViennaRNA/combinatorics.h>

#include  <ViennaRNA/move_set.h>
#include  <ViennaRNA/landscape/paths.h>
#include  <ViennaRNA/landscape/findpath.h>
#include  <ViennaRNA/landscape/move.h>
#include  <ViennaRNA/landscape/neighbor.h>
#include  <ViennaRNA/landscape/walk.h>

#include  <ViennaRNA/mm.h>
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

%include version.i
%include typemaps.i

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
  %template(ElemProbVector) std::vector<vrna_ep_t>;
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


// Typemaps that are independent of scripting language

// This cleans up the char ** array after the function call
%typemap(freearg) char ** {
         free($1);
}


// Additional target language specific typemaps
%include "tmaps.i"

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
%include loops.i
%include basic_algorithms.i
%include mfe.i
%include mfe_window.i
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
%include data_structures.i
%include fold_compound.i


/**********************************************/
/* BEGIN interface for data structures        */
/**********************************************/

%include <ViennaRNA/dp_matrices.h>

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

