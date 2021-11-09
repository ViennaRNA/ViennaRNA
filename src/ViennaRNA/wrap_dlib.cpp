// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*
 *
 *  This is an example illustrating the use the general purpose non-linear
 *  optimization routines from the dlib C++ Library.
 *
 *  The library provides implementations of many popular algorithms such as L-BFGS
 *  and BOBYQA.  These algorithms allow you to find the minimum or maximum of a
 *  function of many input variables.  This example walks though a few of the ways
 *  you might put these routines to use.
 *
 */


#include <iostream>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "ViennaRNA/utils/basic.h"
#ifdef __cplusplus
}
#endif

#include "wrap_dlib.h"

using namespace std;
using namespace dlib;

// ----------------------------------------------------------------------------------------

// In dlib, most of the general purpose solvers optimize functions that take a
// column vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.
typedef matrix<double, 0, 1> column_vector;

/*
 *  function to minimize to obtain equlibrium concentrations
 *  of multistrand systems. We use the transformation
 *
 *  L_a = lambda_a + ln Z_a
 *
 *  such that h(L) reads
 *
 *  h(L) = -\sum_a (c_a L_a - exp(L_a)) + sum_k K_k exp(sum_b L_b A_{b,k}
 *
 *  with total concentration c_a of strand a, equilibrium constant
 *  K_k of strand k, and membership matrix A[b][k] denoting the number
 *  of strands b in complex k
 *
 *  Note, here we minimize h(L) due to implementation issues whereas
 *  in our publication we've written h'(L) = -h(L) to effectively
 *  maximize the function instead.
 */
PRIVATE double
h(const column_vector&  L,
  const double          *eq_constants,
  const double          *concentration_strands_tot,
  const unsigned int    **A,
  size_t                strands,
  size_t                complexes)
{
  double h, hh, *K, maxK;

  K = (double *)vrna_alloc(sizeof(double) * complexes);
  h = 0.;
  maxK = (double)(-INF);

  for (size_t a = 0; a < strands; a++) {
//    printf("L[%u] = %g\n", a, L(a));
    maxK = (maxK < L(a)) ? L(a) : maxK;
  }

  for (size_t k = 0; k < complexes; k++) {
    K[k] = log(eq_constants[k]);

    for (size_t a = 0; a < strands; a++) {
      K[k] += L(a) *
              (double)A[a][k];
    }

    maxK = (maxK < K[k]) ? K[k] : maxK;
  }

  for (size_t a = 0; a < strands; a++)
    h -= concentration_strands_tot[a] *
         L(a);

//  printf("h = %g\n", h);
  hh = 0;

  for (size_t a = 0; a < strands; a++)
    hh += exp(L(a) - maxK);

  for (size_t k = 0; k < complexes; k++)
    hh += exp(K[k] - maxK);

  h += exp(maxK + log(hh));
//  printf("h = %g\n", h);

  free(K);

  return h;
}


/*
 *  Get gradient of h(L)
 */
PRIVATE const column_vector
h_derivative(const column_vector& L,
             const double         *eq_constants,
             const double         *concentration_strands_tot,
             const unsigned int   **A,
             size_t               strands,
             size_t               complexes)
{
  double        *K, *maxK;
  column_vector g(strands);

  K     = (double *)vrna_alloc(sizeof(double) * complexes);
  maxK  = (double *)vrna_alloc(sizeof(double) * strands);

  for (size_t a = 0; a < strands; a++)
    maxK[a] = L(a);

  for (size_t k = 0; k < complexes; k++) {
    K[k] = log(eq_constants[k]);
    for (size_t a = 0; a < strands; a++)
      K[k] += L(a) *
              (double)A[a][k];

    for (size_t a = 0; a < strands; a++)
      if (A[a][k] > 0)
        maxK[a] = (maxK[a] < K[k] + log((double)A[a][k])) ?
                  K[k] + log((double)A[a][k]) :
                  maxK[a];
  }

  for (size_t a = 0; a < strands; a++) {
    g(a) = -concentration_strands_tot[a];

    double hh = exp(L(a) - maxK[a]);
    for (size_t k = 0; k < complexes; k++)
      if (A[a][k] > 0)
        hh += exp(log((double)A[a][k]) +
                  K[k] -
                  maxK[a]);

    g(a) += exp(maxK[a] + log(hh));
//    printf("g(%u) = %g\n", a, g(a));
  }

  free(K);
  free(maxK);

  return g;
}


/*
 *  Get Hessian of h(L)
 */
PRIVATE matrix<double>
h_hessian(const column_vector&  L,
          const double          *eq_constants,
          const unsigned int    **A,
          size_t                strands,
          size_t                complexes)
{
  double                  *K, **xs;

  PRIVATE matrix<double>  H(strands, strands);

  K = (double *)vrna_alloc(sizeof(double) * complexes);
  xs = (double **)vrna_alloc(sizeof(double *) * strands);

  for (size_t a = 0; a < strands; a++) {
    xs[a] = (double *)vrna_alloc(sizeof(double) * strands);
    for (size_t b = 0; b < strands; b++)
      xs[a][b] = (a == b) ? L(a) : (double)(-INF);
  }

  for (size_t k = 0; k < complexes; k++) {
    K[k] = log(eq_constants[k]);
    for (size_t a = 0; a < strands; a++)
      K[k] += L(a) *
              (double)A[a][k];

    for (size_t a = 0; a < strands; a++)
      for (size_t b = 0; b < strands; b++) {
        if ((A[a][k] > 0) && (A[b][k] > 0))
          xs[a][b] = (xs[a][b] < K[k] + log((double)A[a][k]) + log((double)A[b][k])) ?
                      K[k] + log((double)A[a][k]) + log((double)A[b][k]) :
                      xs[a][b];
      }
  }

  for (size_t a = 0; a < strands; a++) {
    for (size_t b = 0; b < strands; b++) {
      double hh = (a == b) ? exp(L(a) - xs[a][b]) : 0.;

      for (size_t k = 0; k < complexes; k++)
        if ((A[a][k] > 0) && (A[b][k] > 0))
          hh += exp(log((double)A[a][k]) +
                    log((double)A[b][k]) +
                    K[k] -
                    xs[a][b]);

      H(a,b) = exp(xs[a][b] + log(hh));
//      H(b,a) = H(a,b);
//      printf("H(%u,%u) = %g\n", a, b, H(a,b));
    }
  }

  free(K);
  for (size_t a = 0; a < strands; a++)
    free(xs[a]);
  free(xs);

  return H;
}


class h_model
{
/*!
 *  This object is a "function model" which can be used with the
 *  find_min_trust_region() routine.
 * !*/

public:
typedef ::column_vector column_vector;
typedef matrix<double> general_matrix;

const double *eq_constants;
const double *concentration_strands_tot;
const unsigned int **A;
size_t strands;
size_t complexes;

double
operator()(const column_vector& x) const
{
  return h(x, eq_constants, concentration_strands_tot, A, strands, complexes);
}


void
init(const double       *eq_constants,
     const double       *concentration_strands_tot,
     const unsigned int **A,
     size_t             strands,
     size_t             complexes)
{
  this->eq_constants              = eq_constants;
  this->concentration_strands_tot = concentration_strands_tot;
  this->A                         = A;
  this->strands                   = strands;
  this->complexes                 = complexes;
}


void
get_derivative_and_hessian(const column_vector& x,
                           column_vector&       der,
                           general_matrix&      hess) const
{
  der   = h_derivative(x, eq_constants, concentration_strands_tot, A, strands, complexes);
  hess  = h_hessian(x, eq_constants, A, strands, complexes);
}
};


/*
 *  Get concentrations of single strands from
 *  a given vector L (that minimizes h(L))
 */
PRIVATE double *
conc_single_strands(const column_vector&  L,
                    size_t                strands)
{
  double *c = (double *)vrna_alloc(sizeof(double) * strands);

  for (size_t a = 0; a < strands; a++)
    c[a] = exp(L(a));

  return c;
}


/*
 *  Get concentrations of complexes from
 *  a given vector L (that minimizes h(L))
 */
PRIVATE double *
conc_complexes(const column_vector& L,
               const double         *eq_const,
               const unsigned int   **A,
               size_t               strands,
               size_t               complexes)
{
  double *c;

  c = (double *)vrna_alloc(sizeof(double) * complexes);

  for (size_t k = 0; k < complexes; k++) {
    c[k] = log(eq_const[k]);

    for (size_t a = 0; a < strands; a++)
      c[k] += (double)A[a][k] * L(a);

    c[k] = exp(c[k]);
  }

  return c;
}


double *
vrna_equilibrium_conc(const double        *eq_constants,
                      double              *concentration_strands,
                      const unsigned int  **A,
                      size_t              num_strands,
                      size_t              num_complexes)
{
  double        *r = NULL;

  column_vector starting_point;

  h_model       h;

  h.init(eq_constants,
         concentration_strands,
         A,
         num_strands,
         num_complexes);

  starting_point.set_size(num_strands);

  for (size_t a = 0; a < num_strands; a++)
    starting_point(a) = 0.;

  find_min_trust_region(objective_delta_stop_strategy(1e-18),
                        h,
                        starting_point,
                        1   // initial trust region radius
                        );

  double *conc_monomers = conc_single_strands(starting_point, num_strands);

  for (size_t a = 0; a < num_strands; a++)
    concentration_strands[a] = conc_monomers[a];

  r = conc_complexes(starting_point, eq_constants, A, num_strands, num_complexes);

  free(conc_monomers);

  return r;
}
