#ifndef VIENNA_RNA_PACKAGE_CONCENTRATIONS_H
#define VIENNA_RNA_PACKAGE_CONCENTRATIONS_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file     concentrations.h
 *  @ingroup  pf_fold cofold pf_cofold
 *  @brief    Concentration computations for RNA-RNA interactions
 */

/**
 *  @{
 *  @ingroup  pf_cofold
 */

/** @brief Typename for the data structure that stores the dimer concentrations, #vrna_dimer_conc_s, as required by vrna_pf_dimer_concentration() */
typedef struct vrna_dimer_conc_s vrna_dimer_conc_t;


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Backward compatibility typedef for #vrna_dimer_conc_s
 */
typedef struct vrna_dimer_conc_s ConcEnt;

#endif

#include <ViennaRNA/params/basic.h>

/**
 *  @brief  Data structure for concentration dependency computations
 */
struct vrna_dimer_conc_s {
  double  Ac_start;   /**< @brief start concentration A */
  double  Bc_start;   /**< @brief start concentration B */
  double  ABc;        /**< @brief End concentration AB */
  double  AAc;
  double  BBc;
  double  Ac;
  double  Bc;
};


/**
 *  @brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 *
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #vrna_dimer_pf_t struct.
 *
 *  @param FcAB       Free energy of AB dimer (FcAB entry)
 *  @param FcAA       Free energy of AA dimer (FcAB entry)
 *  @param FcBB       Free energy of BB dimer (FcAB entry)
 *  @param FEA        Free energy of monomer A
 *  @param FEB        Free energy of monomer B
 *  @param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  @param exp_params The precomputed Boltzmann factors
 *  @return vrna_dimer_conc_t array containing the equilibrium energies and start concentrations
 */
vrna_dimer_conc_t *vrna_pf_dimer_concentrations(double                  FcAB,
                                                double                  FcAA,
                                                double                  FcBB,
                                                double                  FEA,
                                                double                  FEB,
                                                const double            *startconc,
                                                const vrna_exp_param_t  *exp_params);

double *
vrna_equilibrium_constants(const double        *dG_complexes,
                      const double        *dG_strands,
                      const unsigned int  **A,
                      double        kT,
                      size_t        strands,
                      size_t        complexes);

/**
 *  @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 #################################################
 # DEPRECATED FUNCTIONS                          #
 #################################################
 */

/**
 *  @brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 *
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #vrna_dimer_pf_t struct.
 *
 *  @deprecated{ Use vrna_pf_dimer_concentrations() instead!}
 *
 *  @param FEAB       Free energy of AB dimer (FcAB entry)
 *  @param FEAA       Free energy of AA dimer (FcAB entry)
 *  @param FEBB       Free energy of BB dimer (FcAB entry)
 *  @param FEA        Free energy of monomer A
 *  @param FEB        Free energy of monomer B
 *  @param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  @return vrna_dimer_conc_t array containing the equilibrium energies and start concentrations
 */
DEPRECATED(vrna_dimer_conc_t *get_concentrations(double FEAB,
                                                 double FEAA,
                                                 double FEBB,
                                                 double FEA,
                                                 double FEB,
                                                 double *startconc),
          "Use vrna_pf_dimer_concentrations() instead");

#endif

#endif
