/**********************************************/
/* BEGIN interface for RNAxplorer functions   */
/**********************************************/
%module RNAxplorer

%feature("autodoc", "3");

%{
extern "C" {
#include "../src/barrier_lower_bound.h"
#include "../src/distorted_sampling.h"
#include "../src/distorted_samplingMD.h"
#include "../src/dist_class_sc.h"
#include "../src/repellant_sampling.h"
#include "../src/meshpoint.h"
#include "../src/PathFinder.h"
#include "../src/RNAwalk.h"
}
%}

%include "typemaps.i"
%include "std_vector.i";
%include "std_string.i";

%template(StringVector) std::vector<std::string>;
%template(DoubleVector) std::vector<double>;

%include barrier_lower_bound.i
%include distorted_sampling.i
%include distorted_samplingMD.i
%include meshpoint.i
%include RNAwalk.i
%include paths.i


/* add interface for repulsive sampling */

%rename (add_repulsion) rnax_add_repulsion;

%{
  int add_repulsion(vrna_fold_compound_t *fc,
                    const char *structure,
                    double     strength)
  {
    return rnax_add_repulsion(fc, structure, strength);
  }

  int change_repulsion(vrna_fold_compound_t *fc,
                       int                   id,
                       double     strength)
  {
    return rnax_change_repulsion(fc, id, strength);
  }
%}

int add_repulsion(vrna_fold_compound_t *fc, const char *structure, double     strength);
int change_repulsion(vrna_fold_compound_t *fc, int id, double     strength);

%include "../src/repellant_sampling.h"

/* 
rnax_add_repulsion(vrna_fold_compound_t *fc,
                   const char *structure,
                   double     strength)
*/
