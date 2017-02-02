#ifndef VIENNA_RNA_PACKAGE_PARAMS_H
#define VIENNA_RNA_PACKAGE_PARAMS_H

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @file     params.h
 *  @ingroup  energy_parameters
 *  @brief    Functions to deal with sets of energy parameters
 */

/**
 *  @addtogroup energy_parameters
 *  @brief All relevant functions to retrieve and copy pre-calculated energy parameter sets as well as
 *  reading/writing the energy parameter set from/to file(s).
 *
 *  This module covers all relevant functions for pre-calculation of the energy parameters
 *  necessary for the folding routines provided by RNAlib. Furthermore, the energy parameter set
 *  in the RNAlib can be easily exchanged by a user-defined one. It is also possible to write the
 *  current energy parameter set into a text file.
 *  @{
 *  @ingroup  energy_parameters
 */

/** @brief Typename for the free energy parameter data structure #vrna_params */
typedef struct  vrna_param_s       vrna_param_t;
/** @brief Typename for the Boltzmann factor data structure #vrna_exp_params */
typedef struct  vrna_exp_param_s   vrna_exp_param_t;

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
 *  @brief The datastructure that contains temperature scaled energy parameters.
 */
struct vrna_param_s {
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

  double  temperature;            /**<  @brief  Temperature used for loop contribution scaling */

  vrna_md_t model_details;   /**<  @brief  Model details to be used in the recursions */
};

/**
 *  @brief  The data structure that contains temperature scaled Boltzmann weights of the energy parameters.
 */
struct vrna_exp_param_s {
  int     id;   /**<  @brief  An identifier for the data structure
                      @deprecated This attribute will be removed in version 3
                */
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
  double  pf_scale;     /**<  @brief    Scaling factor to avoid over-/underflows */

  double  temperature;  /**<  @brief    Temperature used for loop contribution scaling */
  double  alpha;        /**<  @brief    Scaling factor for the thermodynamic temperature
                              @details  This allows for temperature scaling in Boltzmann
                                        factors independently from the energy contributions.
                                        The resulting Boltzmann factors are then computed by
                                        @f$ e^{-E/(\alpha \cdot K \cdot T)} @f$
                        */

  vrna_md_t model_details; /**<  @brief  Model details to be used in the recursions */

};


/**
 *  @brief  Get a data structure containing prescaled free energy parameters
 *
 *  If a NULL pointer is passed for the model details parameter, the default
 *  model parameters are stored within the requested #vrna_param_t structure.
 *
 *  @see #vrna_md_t, vrna_md_set_default(), vrna_exp_params()
 *
 *  @param  md  A pointer to the model details to store inside the structure (Maybe NULL)
 *  @return     A pointer to the memory location where the requested parameters are stored
 */
vrna_param_t *
vrna_params(vrna_md_t *md);

/**
 *  @brief Get a copy of the provided free energy parameters
 *
 *  If NULL is passed as parameter, a default set of energy parameters is created
 *  and returned.
 *
 *  @see vrna_params(), #vrna_param_t
 *
 *  @param  par   The free energy parameters that are to be copied (Maybe NULL)
 *  @return       A copy or a default set of the (provided) parameters
 */
vrna_param_t *
vrna_params_copy(vrna_param_t *par);

/**
 *  @brief  Get a data structure containing prescaled free energy parameters
 *          already transformed to Boltzmann factors
 *
 *  This function returns a data structure that contains all necessary precomputed
 *  energy contributions for each type of loop.
 *
 *  In contrast to vrna_params(), the free energies within this data structure
 *  are stored as their Boltzmann factors, i.e.
 *
 *  @f$ exp(-E / kT) @f$
 *
 *  where @f$ E @f$ is the free energy.
 *
 *  If a NULL pointer is passed for the model details parameter, the default
 *  model parameters are stored within the requested #vrna_exp_param_t structure.
 *
 *  @see #vrna_md_t, vrna_md_set_default(), vrna_params(), vrna_rescale_pf_params()
 *
 *  @param  md  A pointer to the model details to store inside the structure (Maybe NULL)
 *  @return     A pointer to the memory location where the requested parameters are stored
 */
vrna_exp_param_t *
vrna_exp_params(vrna_md_t *md);

/**
 *  @brief  Get a data structure containing prescaled free energy parameters
 *          already transformed to Boltzmann factors (alifold version)
 *
 *  If a NULL pointer is passed for the model details parameter, the default
 *  model parameters are stored within the requested #vrna_exp_param_t structure.
 *
 *  @see #vrna_md_t, vrna_md_set_default(), vrna_exp_params(), vrna_params()
 *
 *  @param  n_seq   The number of sequences in the alignment
 *  @param  md      A pointer to the model details to store inside the structure (Maybe NULL)
 *  @return         A pointer to the memory location where the requested parameters are stored
 */
vrna_exp_param_t *
vrna_exp_params_comparative(unsigned int n_seq,
                            vrna_md_t *md);

/**
 *  @brief Get a copy of the provided free energy parameters (provided as Boltzmann factors)
 *
 *  If NULL is passed as parameter, a default set of energy parameters is created
 *  and returned.
 *
 *  @see vrna_exp_params(), #vrna_exp_param_t
 *
 *  @param  par   The free energy parameters that are to be copied (Maybe NULL)
 *  @return       A copy or a default set of the (provided) parameters
 */
vrna_exp_param_t *
vrna_exp_params_copy(vrna_exp_param_t *par);

/**
 *  @brief  Update/Reset energy parameters data structure within a #vrna_fold_compound_t
 *
 *  Passing NULL as second argument leads to a reset of the energy parameters within
 *  vc to their default values. Otherwise, the energy parameters provided will be copied
 *  over into vc.
 *
 *  @see vrna_params_reset(), #vrna_param_t, #vrna_md_t, vrna_params()
 *
 *  @param  vc    The #vrna_fold_compound_t that is about to receive updated energy parameters
 *  @param  par   The energy parameters used to substitute those within vc (Maybe NULL)
 */
void
vrna_params_subst( vrna_fold_compound_t *vc,
                    vrna_param_t *par);

/**
 *  @brief Update the energy parameters for subsequent partition function computations
 *
 *  This function can be used to properly assign new energy parameters for partition
 *  function computations to a #vrna_fold_compound_t. For this purpose, the data of the
 *  provided pointer `params`  will be copied into `vc` and a recomputation of the partition
 *  function scaling factor is issued, if the `pf_scale` attribute of `params` is less than `1.0`.
 *
 *  Passing NULL as second argument leads to a reset of the energy parameters within
 *  vc to their default values
 *
 *  @see vrna_exp_params_reset(), vrna_exp_params_rescale(), #vrna_exp_param_t, #vrna_md_t,
 *  vrna_exp_params()
 *
 *  @param  vc      The fold compound data structure
 *  @param  params  A pointer to the new energy parameters
 */
void
vrna_exp_params_subst(vrna_fold_compound_t *vc,
                      vrna_exp_param_t *params);

/**
 *  @brief Rescale Boltzmann factors for partition function computations
 *
 *  This function may be used to (automatically) rescale the Boltzmann factors used
 *  in partition function computations. Since partition functions over subsequences
 *  can easily become extremely large, the RNAlib internally rescales them to avoid
 *  numerical over- and/or underflow. Therefore, a proper scaling factor @f$s@f$ needs to
 *  be chosen that in turn is then used to normalize the corresponding
 *  partition functions @f$\hat{q}[i,j] = q[i,j] / s^{(j-i+1)}@f$.
 *
 *  This function provides two ways to automatically adjust the scaling
 *  factor.
 *  1. Automatic guess
 *  2. Automatic adjustment according to MFE
 *
 *  Passing `NULL` as second parameter activates the _automatic guess mode_. Here,
 *  the scaling factor is recomputed according to a mean free energy of `184.3*length` cal
 *  for random sequences.
 *  @note This recomputation only takes place if the `pf_scale` attribute of the
 *  `exp_params` data structure contained in `vc` has a value below `1.0`.
 *
 *  On the other hand, if the MFE for a sequence is known, it can be used to recompute
 *  a more robust scaling factor, since it represents the lowest free energy of the entire
 *  ensemble of structures, i.e. the highest Boltzmann factor. To activate this second
 *  mode of _automatic adjustment according to MFE_, a pointer to the MFE value needs to
 *  be passed as second argument. This value is then taken to compute the scaling factor
 *  as @f$ s = exp((sfact * MFE) / kT / length )@f$, where sfact is an additional
 *  scaling weight located in the vrna_md_t data structure of `exp_params` in `vc`.
 *
 *  The computed scaling factor @f$s@f$ will be stored as `pf_scale` attribute of the
 *  `exp_params` data structure in `vc`.
 *
 *  @see vrna_exp_params_subst(), vrna_md_t, vrna_exp_param_t, #vrna_fold_compound_t
 *
 *  @param  vc  The fold compound data structure
 *  @param  mfe A pointer to the MFE (in kcal/mol) or NULL
 */
void
vrna_exp_params_rescale(vrna_fold_compound_t *vc,
                        double *mfe);

/**
 *  @brief  Reset free energy parameters within a #vrna_fold_compound_t
 *          according to provided, or default model details
 *
 *  This function allows one to rescale free energy parameters for subsequent structure
 *  prediction or evaluation according to a set of model details, e.g. temperature
 *  values. To do so, the caller provides either a pointer to a set of model details
 *  to be used for rescaling, or NULL if global default setting should be used.
 *
 *  @see vrna_exp_params_reset(), vrna_params_subs()
 *  @param  vc    The fold compound data structure
 *  @param  md_p  A pointer to the new model details (or NULL for reset to defaults)
 */
void vrna_params_reset( vrna_fold_compound_t *vc,
                        vrna_md_t *md_p);

/**
 *  @brief  Reset Boltzmann factors for partition function computations
 *          within a #vrna_fold_compound_t according to provided, or
 *          default model details
 *
 *  This function allows one to rescale Boltzmann factors for subsequent partition
 *  function computations according to a set of model details, e.g. temperature
 *  values. To do so, the caller provides either a pointer to a set of model details
 *  to be used for rescaling, or NULL if global default setting should be used.
 *
 *  @see vrna_params_reset(), vrna_exp_params_subst(), vrna_exp_params_rescale()
 *  @param  vc    The fold compound data structure
 *  @param  md_p  A pointer to the new model details (or NULL for reset to defaults)
 */
void vrna_exp_params_reset( vrna_fold_compound_t *vc,
                            vrna_md_t *md_p);

#ifdef  VRNA_BACKWARD_COMPAT

/**
 *  @brief Old typename of #vrna_param_s
 *  @deprecated Use #vrna_param_t instead!
*/
typedef struct vrna_param_s     paramT;

/**
 *  @brief Old typename of #vrna_exp_param_s
 *  @deprecated Use #vrna_exp_param_t instead!
*/
typedef struct vrna_exp_param_s pf_paramT;

DEPRECATED(vrna_param_t *get_parameter_copy(vrna_param_t *par));

/**
 *  get a data structure of type @ref vrna_exp_param_t which contains
 *  the Boltzmann weights of several energy parameters scaled
 *  according to the current temperature
 *
 *  @deprecated Use vrna_exp_params() instead!
 *
 *  @return The data structure containing Boltzmann weights for use in partition function calculations
 */
DEPRECATED(vrna_exp_param_t *get_scaled_pf_parameters(void));

/**
 *  @brief Get precomputed Boltzmann factors of the loop type
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
 *  @f$ exp(-\Delta G / kT) @f$
 *
 *  @deprecated Use vrna_exp_params() instead!
 *
 *  @see get_scaled_pf_parameters(), get_boltzmann_factor_copy()
 *
 *  @param  temperature   The temperature in degrees Celcius used for (re-)scaling the energy contributions
 *  @param  betaScale     A scaling value that is used as a multiplication factor for the absolute
 *                        temperature of the system
 *  @param  md            The model details to be used
 *  @param  pf_scale      The scaling factor for the Boltzmann factors
 *  @return               A set of precomputed Boltzmann factors
 */
DEPRECATED(vrna_exp_param_t *get_boltzmann_factors(double temperature, double betaScale, vrna_md_t md, double pf_scale));

/**
 *  @brief Get a copy of already precomputed Boltzmann factors
 *
 *  @deprecated Use vrna_exp_params_copy() instead!
 *
 *  @see get_boltzmann_factors(), get_scaled_pf_parameters()
 *
 *  @param  parameters  The input data structure that shall be copied
 *  @return             A copy of the provided Boltzmann factor data set
 */
DEPRECATED(vrna_exp_param_t *get_boltzmann_factor_copy(vrna_exp_param_t *parameters));

/**
 *  @brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions (alifold variant)
 *
 *  @deprecated Use vrna_exp_params_comparative() instead!
 *
 */
DEPRECATED(vrna_exp_param_t *get_scaled_alipf_parameters(unsigned int n_seq));

/**
 *  @brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions (alifold variant) with
 *  independent thermodynamic temperature
 *
 *  @deprecated Use vrna_exp_params_comparative() instead!
 *
 */
DEPRECATED(vrna_exp_param_t *get_boltzmann_factors_ali(unsigned int n_seq, double temperature, double betaScale, vrna_md_t md, double pf_scale));

/**
 * @brief Get precomputed energy contributions for all the known loop types
 *
 *  @note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See get_scaled_parameters() for a completely threadsafe
 *        implementation.
 *
 *  @deprecated Use vrna_params() instead!
 *
 * @return     A set of precomputed energy contributions
 */
DEPRECATED(vrna_param_t *scale_parameters(void));

/**
 * @brief Get precomputed energy contributions for all the known loop types
 *
 *  Call this function to retrieve precomputed energy contributions, i.e. scaled
 *  according to the temperature passed. Furthermore, this function assumes a
 *  data structure that contains the model details as well, such that subsequent
 *  folding recursions are able to retrieve the correct model settings
 *
 *  @deprecated Use vrna_params() instead!
 *
 *  @see #vrna_md_t, set_model_details()
 *
 *  @param temperature  The temperature in degrees Celcius
 *  @param md           The model details
 *  @return             precomputed energy contributions and model settings
 */
DEPRECATED(vrna_param_t *get_scaled_parameters(double temperature, vrna_md_t md));

DEPRECATED(vrna_param_t     *copy_parameters(void));
DEPRECATED(vrna_param_t     *set_parameters(vrna_param_t *dest));
DEPRECATED(vrna_exp_param_t *scale_pf_parameters(void));
DEPRECATED(vrna_exp_param_t *copy_pf_param(void));
DEPRECATED(vrna_exp_param_t *set_pf_param(vrna_param_t *dest));

#endif

/**
 *  @}
 */



#endif
