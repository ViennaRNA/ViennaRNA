%module RNA
%{
#include  "../H/utils.h"
#include  "../H/fold_vars.h"
#undef fold
#include  "../H/fold.h"
#include  "../H/part_func.h"
#include  "../H/PS_dot.h"
#include  "../H/inverse.h"
#include  "../H/RNAstruct.h"
#include  "../H/treedist.h"
#include  "../H/stringdist.h"
#include  "../H/profiledist.h"
#include  "../H/dist_vars.h"
%}
//
%title "Interface to the Vienna RNA library"
%section "Folding Routines"
%subsection "Minimum free Energy Folding"
%include  "../H/fold.h"
%subsection "Partition function Folding"
%include  "../H/part_func.h"
%subsection "Inverse Folding"
%include  "../H/inverse.h"
%subsection "Global Variables to Modify Folding"
extern float *pr;  /*  base pairing prob. matrix */
%include  "../H/fold_vars.h"
%{
float get_pr(int i, int j) {
   int ii;
  if (i>j) {ii=i; i=j; j=ii;}
  return pr[iindx[i]-j];
}
%}
float get_pr(int i, int j);
/* Get probability of pair i.j from the pr array */

%section "Parsing and Comparing Structures"
%include  "../H/RNAstruct.h"
%include  "../H/treedist.h"
%include  "../H/stringdist.h"
%include  "../H/profiledist.h"
// from dist_vars.h
extern int   edit_backtrack;  /* set to 1 if you want backtracking */ 
extern char *aligned_line[2]; /* containes alignment after backtracking */
extern int  cost_matrix;      /* 0 usual costs (default), 1 Shapiro's costs */

%section "Utilities"
%include  "../H/utils.h"
%include  "../H/PS_dot.h"

// from read_epars.c
extern void  read_parameter_file(const char fname[]);
/* read energy parameters from file */
extern void  write_parameter_file(const char fname[]);
/* write energy parameters to file */
//%include array.i
%include pointer.i
%inline %{
void *deref_any(void **ptr, int index) {
   /* dereference arbitray pointer */
   return (void *) ptr[index];
}
%}
