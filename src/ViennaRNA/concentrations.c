/*
 *                concentrations for RNA-RNA interaction
 *
 *                Ivo L Hofacker
 *                Stephan Bernhart
 *                Ronny Lorenz
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/concentrations.h"

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
PRIVATE double *
Newton_Conc(double  ZAB,
            double  ZAA,
            double  ZBB,
            double  concA,
            double  concB);


PRIVATE double *
conc_single_strands(double  *L,
                    size_t  strands);


PRIVATE double *
conc_complexes(double       *L,
               double       *eq_const,
               unsigned int **A,
               size_t       strands,
               size_t       complexes) VRNA_UNUSED;


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_dimer_conc_t *
vrna_pf_dimer_concentrations(double                 FcAB,
                             double                 FcAA,
                             double                 FcBB,
                             double                 FEA,
                             double                 FEB,
                             const double           *startconc,
                             const vrna_exp_param_t *exp_params)
{
  /*
   * takes an array of start concentrations, computes equilibrium concentrations of dimers,
   * monomers, returns array of concentrations in strucutre vrna_dimer_conc_t
   */
  double            *ConcVec;
  int               i;
  vrna_dimer_conc_t *Concentration;
  double            KAA, KAB, KBB, kT;

  kT            = exp_params->kT / 1000.;
  Concentration = (vrna_dimer_conc_t *)vrna_alloc(20 * sizeof(vrna_dimer_conc_t));
  /*
   * Compute equilibrium constants
   * again note the input free energies are not from the null model (without DuplexInit)
   */

  KAA = exp((2.0 * FEA - FcAA) / kT);
  KBB = exp((2.0 * FEB - FcBB) / kT);
  KAB = exp((FEA + FEB - FcAB) / kT);
  /* printf("Kaa..%g %g %g\n", KAA, KBB, KAB); */
  for (i = 0; ((startconc[i] != 0) || (startconc[i + 1] != 0)); i += 2) {
    ConcVec = Newton_Conc(KAB,
                          KAA,
                          KBB,
                          startconc[i],
                          startconc[i + 1]);
    Concentration[i / 2].Ac_start = startconc[i];
    Concentration[i / 2].Bc_start = startconc[i + 1];
    Concentration[i / 2].ABc      = ConcVec[0];
    Concentration[i / 2].AAc      = ConcVec[1];
    Concentration[i / 2].BBc      = ConcVec[2];
    Concentration[i / 2].Ac       = ConcVec[3];
    Concentration[i / 2].Bc       = ConcVec[4];

    if (!(((i + 2) / 2) % 20))
      Concentration =
        (vrna_dimer_conc_t *)vrna_realloc(Concentration,
                                          ((i + 2) / 2 + 20) * sizeof(vrna_dimer_conc_t));

    free(ConcVec);
  }

  return Concentration;
}


PRIVATE double *
Newton_Conc(double  KAB,
            double  KAA,
            double  KBB,
            double  concA,
            double  concB)
{
  double  TOL, EPS, xn, yn, det, cA, cB, *ConcVec;
  int     i;

  i = 0;
  /* Newton iteration for computing concentrations */
  cA      = concA;
  cB      = concB;
  TOL     = 1e-6;                                     /* Tolerance for convergence */
  ConcVec = (double *)vrna_alloc(5 * sizeof(double)); /* holds concentrations */
  do {
    /* det = (4.0 * KAA * cA + KAB *cB + 1.0) * (4.0 * KBB * cB + KAB *cA + 1.0) - (KAB *cB) * (KAB *cA); */
    det = 1 + 16. * KAA * KBB * cA * cB + KAB * (cA + cB) + 4. * KAA * cA + 4. * KBB * cB + 4. *
          KAB * (KBB * cB * cB + KAA * cA * cA);
    /*
     * xn  = ( (2.0 * KBB * cB*cB + KAB *cA *cB + cB - concB) * (KAB *cA) -
     *       (2.0 * KAA * cA*cA + KAB *cA *cB + cA - concA) * (4.0 * KBB * cB + KAB *cA + 1.0) ) / det;
     */
    xn = ((2.0 * KBB * cB * cB + cB - concB) * (KAB * cA) - KAB * cA * cB * (4. * KBB * cB + 1.) -
          (2.0 * KAA * cA * cA + cA - concA) * (4.0 * KBB * cB + KAB * cA + 1.0)) / det;
    /*
     * yn  = ( (2.0 * KAA * cA*cA + KAB *cA *cB + cA - concA) * (KAB *cB) -
     *       (2.0 * KBB * cB*cB + KAB *cA *cB + cB - concB) * (4.0 * KAA * cA + KAB *cB + 1.0) ) / det;
     */
    yn = ((2.0 * KAA * cA * cA + cA - concA) * (KAB * cB) - KAB * cA * cB * (4. * KAA * cA + 1.) -
          (2.0 * KBB * cB * cB + cB - concB) * (4.0 * KAA * cA + KAB * cB + 1.0)) / det;
    EPS = fabs(xn / cA) + fabs(yn / cB);
    cA  += xn;
    cB  += yn;
    i++;
    if (i > 10000) {
      vrna_log_warning("Newton did not converge after %d steps!!", i);
      break;
    }
  } while (EPS > TOL);

  ConcVec[0]  = cA * cB * KAB;  /* AB concentration */
  ConcVec[1]  = cA * cA * KAA;  /* AA concentration */
  ConcVec[2]  = cB * cB * KBB;  /* BB concentration */
  ConcVec[3]  = cA;             /* A concentration */
  ConcVec[4]  = cB;             /* B concentration */

  return ConcVec;
}


/*
 *  Compute the equilibrium constants for each complex
 *  from its respective ensemble free energy and the
 *  ensemble free energies of the single strands
 */
PUBLIC double *
vrna_equilibrium_constants(const double       *dG_complexes,
                           const double       *dG_strands,
                           const unsigned int **A,
                           double             kT,
                           size_t             strands,
                           size_t             complexes)
{
  double *K, tmp;

  K = (double *)vrna_alloc(sizeof(double) * complexes);

  for (size_t k = 0; k < complexes; k++) {
    tmp = 0;
    for (size_t a = 0; a < strands; a++)
      tmp += (double)A[a][k] * dG_strands[a];

    K[k] = exp((tmp - dG_complexes[k]) / kT);
  }

  return K;
}


/*
 *  Get concentrations of single strands from
 *  a given vector L (that minimizes h(L))
 */
PRIVATE double *
conc_single_strands(double  *L,
                    size_t  strands)
{
  double *c = (double *)vrna_alloc(sizeof(double) * strands);

  for (size_t a = 0; a < strands; a++)
    c[a] = exp(L[a]);

  return c;
}


/*
 *  Get concentrations of complexes from
 *  a given vector L (that minimizes h(L))
 */
PRIVATE double *
conc_complexes(double       *L,
               double       *eq_const,
               unsigned int **A,
               size_t       strands,
               size_t       complexes)
{
  double *c, *c_free;

  c       = (double *)vrna_alloc(sizeof(double) * complexes);
  c_free  = conc_single_strands(L, strands);

  for (size_t k = 0; k < complexes; k++) {
    c[k] = eq_const[k];

    for (size_t a = 0; a < strands; a++)
      c[k] *= pow(c_free[a], (double)A[a][k]);
  }

  free(c_free);

  return c;
}
