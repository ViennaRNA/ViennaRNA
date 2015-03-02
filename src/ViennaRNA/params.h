#ifndef __VIENNA_RNA_PACKAGE_PARAMS_H__
#define __VIENNA_RNA_PACKAGE_PARAMS_H__

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  \addtogroup energy_parameters
 *  \brief All relevant functions to retrieve and copy precalculated energy parameter sets as well as
 *  reading/writing the energy parameter set from/to file(s).
 *
 *  This module covers all relevant functions for precalculation of the energy parameters
 *  necessary for the folding routines provided by RNAlib. Furthermore, the energy parameter set
 *  in the RNAlib can be easily exchanged by a user-defined one. It is also possible to write the
 *  current energy parameter set into a text file.
 *  @{
 *
 *  \file params.h
 */

#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/model.h>

#define   VRNA_GQUAD_MAX_STACK_SIZE     7
#define   VRNA_GQUAD_MIN_STACK_SIZE     2
#define   VRNA_GQUAD_MAX_LINKER_LENGTH  15
#define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
#define   VRNA_GQUAD_MIN_BOX_SIZE       ((4*VRNA_GQUAD_MIN_STACK_SIZE)+(3*VRNA_GQUAD_MIN_LINKER_LENGTH))
#define   VRNA_GQUAD_MAX_BOX_SIZE       ((4*VRNA_GQUAD_MAX_STACK_SIZE)+(3*VRNA_GQUAD_MAX_LINKER_LENGTH))

/**
 *  \brief The datastructure that contains temperature scaled energy parameters.
 */
typedef struct vrna_param_t{
  int     id;
  int     stack[NBPAIRS+1][NBPAIRS+1];
  int     hairpin[31];
  int     bulge[MAXLOOP+1];
  int     internal_loop[MAXLOOP+1];
  int     mismatchExt[NBPAIRS+1][5][5];
  int     mismatchI[NBPAIRS+1][5][5];
  int     mismatch1nI[NBPAIRS+1][5][5];
  int     mismatch23I[NBPAIRS+1][5][5];
  int     mismatchH[NBPAIRS+1][5][5];
  int     mismatchM[NBPAIRS+1][5][5];
  int     dangle5[NBPAIRS+1][5];
  int     dangle3[NBPAIRS+1][5];
  int     int11[NBPAIRS+1][NBPAIRS+1][5][5];
  int     int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  int     int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  int     ninio[5];
  double  lxc;
  int     MLbase;
  int     MLintern[NBPAIRS+1];
  int     MLclosing;
  int     TerminalAU;
  int     DuplexInit;
  int     Tetraloop_E[200];
  char    Tetraloops[1401];
  int     Triloop_E[40];
  char    Triloops[241];
  int     Hexaloop_E[40];
  char    Hexaloops[1801];
  int     TripleC;
  int     MultipleCA;
  int     MultipleCB;
  int     gquad [VRNA_GQUAD_MAX_STACK_SIZE + 1]
                [3*VRNA_GQUAD_MAX_LINKER_LENGTH + 1];

  double  temperature;            /**<  \brief  Temperature used for loop contribution scaling */

  vrna_md_t model_details;   /**<  \brief  Model details to be used in the recursions */

} vrna_param_t;

/**
 *  \brief  The datastructure that contains temperature scaled Boltzmann weights of the energy parameters.
 */
typedef struct vrna_exp_param_t{
  int     id;
  double  expstack[NBPAIRS+1][NBPAIRS+1];
  double  exphairpin[31];
  double  expbulge[MAXLOOP+1];
  double  expinternal[MAXLOOP+1];
  double  expmismatchExt[NBPAIRS+1][5][5];
  double  expmismatchI[NBPAIRS+1][5][5];
  double  expmismatch23I[NBPAIRS+1][5][5];
  double  expmismatch1nI[NBPAIRS+1][5][5];
  double  expmismatchH[NBPAIRS+1][5][5];
  double  expmismatchM[NBPAIRS+1][5][5];
  double  expdangle5[NBPAIRS+1][5];
  double  expdangle3[NBPAIRS+1][5];
  double  expint11[NBPAIRS+1][NBPAIRS+1][5][5];
  double  expint21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  double  expint22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  double  expninio[5][MAXLOOP+1];
  double  lxc;
  double  expMLbase;
  double  expMLintern[NBPAIRS+1];
  double  expMLclosing;
  double  expTermAU;
  double  expDuplexInit;
  double  exptetra[40];
  double  exptri[40];
  double  exphex[40];
  char    Tetraloops[1401];
  double  expTriloop[40];
  char    Triloops[241];
  char    Hexaloops[1801];
  double  expTripleC;
  double  expMultipleCA;
  double  expMultipleCB;
  double  expgquad[VRNA_GQUAD_MAX_STACK_SIZE + 1]
                  [3*VRNA_GQUAD_MAX_LINKER_LENGTH + 1];

  double  kT;
  double  pf_scale;     /**<  \brief    Scaling factor to avoid over-/underflows */

  double  temperature;  /**<  \brief    Temperature used for loop contribution scaling */
  double  alpha;        /**<  \brief    Scaling factor for the thermodynamic temperature
                              \details  This allows for temperature scaling in Boltzmann
                                        factors independently from the energy contributions.
                                        The resulting Boltzmann factors are then computed by
                                        \f$ e^{-E/(\alpha \cdot K \cdot T)} \f$
                        */

  vrna_md_t model_details; /**<  \brief  Model details to be used in the recursions */

} vrna_exp_param_t;

/**
 * \brief Get precomputed energy contributions for all the known loop types
 *
 *  \note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See get_scaled_parameters() for a completely threadsafe
 *        implementation.
 *
 * \return     A set of precomputed energy contributions
 */
vrna_param_t *scale_parameters(void);


vrna_param_t *vrna_get_energy_contributions(vrna_md_t md);


/**
 * \brief Get precomputed energy contributions for all the known loop types
 *
 *  Call this function to retrieve precomputed energy contributions, i.e. scaled
 *  according to the temperature passed. Furthermore, this function assumes a
 *  data structure that contains the model details as well, such that subsequent
 *  folding recursions are able to retrieve the correct model settings
 *
 *  \see #vrna_md_t, set_model_details()
 *
 *  \param temperature  The temperature in degrees Celcius
 *  \param md           The model details
 *  \return             precomputed energy contributions and model settings
 */
vrna_param_t *get_scaled_parameters(double temperature,
                              vrna_md_t md);

vrna_param_t *get_parameter_copy(vrna_param_t *par);

/**
 *  get a datastructure of type \ref vrna_exp_param_t which contains
 *  the Boltzmann weights of several energy parameters scaled
 *  according to the current temperature
 *  \return The datastructure containing Boltzmann weights for use in partition function calculations
 */
vrna_exp_param_t *get_scaled_pf_parameters(void);

vrna_exp_param_t *vrna_get_boltzmann_factors(vrna_md_t md);

/**
 *  \brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions with independent thermodynamic
 *  temperature
 *
 *  This function returns a data structure that contains
 *  all necessary precalculated Boltzmann factors for each
 *  loop type contribution.<br>
 *  In contrast to get_scaled_pf_parameters(), this function
 *  enables setting of independent temperatures for both, the
 *  individual energy contributions as well as the thermodynamic
 *  temperature used in
 *  \f$ exp(-\Delta G / kT) \f$
 *
 *  \see get_scaled_pf_parameters(), get_boltzmann_factor_copy()
 *
 *  \param  temperature   The temperature in degrees Celcius used for (re-)scaling the energy contributions
 *  \param  betaScale     A scaling value that is used as a multiplication factor for the absolute
 *                        temperature of the system
 *  \param  md            The model details to be used
 *  \param  pf_scale      The scaling factor for the Boltzmann factors
 *  \return               A set of precomputed Boltzmann factors
 */
vrna_exp_param_t *get_boltzmann_factors( double temperature,
                                  double betaScale,
                                  vrna_md_t md,
                                  double pf_scale);

/**
 *  \brief Get a copy of already precomputed Boltzmann factors
 *
 *  \see get_boltzmann_factors(), get_scaled_pf_parameters()
 *
 *  \param  parameters  The input data structure that shall be copied
 *  \return             A copy of the provided Boltzmann factor dataset
 */
vrna_exp_param_t *get_boltzmann_factor_copy(vrna_exp_param_t *parameters);

/**
 *  \brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions (alifold variant)
 *
 */
vrna_exp_param_t *get_scaled_alipf_parameters(unsigned int n_seq);

/**
 *  \brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions (alifold variant) with
 *  independent thermodynamic temperature
 *
 */
vrna_exp_param_t *get_boltzmann_factors_ali( unsigned int n_seq,
                                      double temperature,
                                      double betaScale,
                                      vrna_md_t md,
                                      double pf_scale);

vrna_exp_param_t *vrna_get_boltzmann_factors_ali(unsigned int n_seq,
                                          vrna_md_t md);

/**
 *  @}
 */

#ifdef  VRNA_BACKWARD_COMPAT

#define paramT      vrna_param_t        /* restore compatibility of struct rename */
#define pf_paramT   vrna_exp_param_t    /* restore compatibility of struct rename */

DEPRECATED(vrna_param_t     *copy_parameters(void));
DEPRECATED(vrna_param_t     *set_parameters(vrna_param_t *dest));
DEPRECATED(vrna_exp_param_t *scale_pf_parameters(void));
DEPRECATED(vrna_exp_param_t *copy_pf_param(void));
DEPRECATED(vrna_exp_param_t *set_pf_param(vrna_param_t *dest));

#endif



#endif
