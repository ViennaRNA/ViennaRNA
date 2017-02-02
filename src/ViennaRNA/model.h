#ifndef VIENNA_RNA_PACKAGE_MODEL_H
#define VIENNA_RNA_PACKAGE_MODEL_H

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
 *  @file     model.h
 *  @ingroup  model_details
 *  @brief    The model details data structure and its corresponding modifiers
 */

/**
 *  @{
 *  @ingroup   model_details
 */

#ifndef NBASES
#define NBASES 8
#endif

/** @brief Typename for the model details data structure #vrna_md_s */
typedef struct vrna_md_s  vrna_md_t;

/**
 *  @brief
 *  @htmlonly Default temperature for structure prediction and free energy evaluation in &#176C @endhtmlonly
 *  @latexonly Default temperature for structure prediction and free energy evaluation in $^\circ C$ @endlatexonly
 *  @see  #vrna_md_t.temperature, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0

/**
 *  @brief  Default scaling factor for partition function computations
 *  @see  #vrna_exp_param_t.pf_scale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_PF_SCALE       -1

/**
 *  @brief  Default scaling factor for absolute thermodynamic temperature in Boltzmann factors
 *  @see    #vrna_exp_param_t.alpha, #vrna_md_t.betaScale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BETA_SCALE     1.

/** @brief  Default dangling end model
 *  @see  #vrna_md_t.dangles, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_DANGLES        2

/**
 *  @brief  Default model behavior for lookup of special tri-, tetra-, and hexa-loops
 *  @see    #vrna_md_t.special_hp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_SPECIAL_HP     1

/**
 *  @brief  Default model behavior for so-called 'lonely pairs'
 *  @see    #vrna_md_t.noLP, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_LP          0

/**
 *  @brief  Default model behavior for G-U base pairs
 *  @see    #vrna_md_t.noGU, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU          0

/**
 *  @brief  Default model behavior for G-U base pairs closing a loop
 *  @see    #vrna_md_t.noGUclosure, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU_CLOSURE  0

/**
 *  @brief  Default model behavior to treat a molecule as a circular RNA (DNA)
 *  @see    #vrna_md_t.circ, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_CIRC           0

/**
 *  @brief  Default model behavior regarding the treatment of G-Quadruplexes
 *  @see    #vrna_md_t.gquad, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_GQUAD          0

#define VRNA_MODEL_DEFAULT_CANONICAL_BP   0

/**
 *  @brief  Default behavior of the model regarding unique multi-branch loop decomposition
 *  @see    #vrna_md_t.uniq_ML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_UNIQ_ML        0

/**
 *  @brief  Default model behavior on which energy set to use
 *  @see    #vrna_md_t.energy_set, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ENERGY_SET     0

/**
 *  @brief  Default model behavior with regards to backtracking of structures
 *  @see    #vrna_md_t.backtrack, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK      1

/**
 *  @brief  Default model behavior on what type of backtracking to perform
 *  @see    #vrna_md_t.backtrack_type, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK_TYPE 'F'

/**
 *  @brief  Default model behavior with regards to computing base pair probabilities
 *  @see    #vrna_md_t.compute_bpp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_COMPUTE_BPP    1

/**
 *  @brief  Default model behavior for the allowed maximum base pair span
 *  @see    #vrna_md_t.max_bp_span, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_MAX_BP_SPAN    -1

/**
 *  @brief  Default model behavior for the sliding window approach
 *  @see    #vrna_md_t.window_size, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_WINDOW_SIZE    -1

/**
 *  @brief  Default model behavior on how to evaluate the energy contribution of multi-branch loops
 *  @see    #vrna_md_t.logML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_LOG_ML         0

/**
 *  @brief  Default model behavior for consensus structure energy evaluation
 *  @see    #vrna_md_t.oldAliEn, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_OLD_EN     0

/**
 *  @brief  Default model behavior for consensus structure co-variance contribution assessment
 *  @see    #vrna_md_t.ribo, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_RIBO       0

/**
 *  @brief  Default model behavior for weighting the co-variance score in consensus structure prediction
 *  @see    #vrna_md_t.cv_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_CV_FACT    1.

/** @brief  Default model behavior for weighting the nucleotide conservation? in consensus structure prediction
 *  @see    #vrna_md_t.nc_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_NC_FACT    1.


#ifdef  VRNA_BACKWARD_COMPAT

#ifndef MAXALPHA
/**
 *  @brief Maximal length of alphabet
 */
#define MAXALPHA              20
#endif

#endif

/**
 *  @brief The data structure that contains the complete model details used throughout the calculations
 *
 *  For convenience reasons, we provide the type name #vrna_md_t to address this data structure
 *  without the use of the struct keyword
 *
 *  @see  vrna_md_set_default(), set_model_details(), vrna_md_update(), #vrna_md_t
 *
 */
struct vrna_md_s {
  double  temperature;                  /**<  @brief  The temperature used to scale the thermodynamic parameters */
  double  betaScale;                    /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
  int     dangles;                      /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)

                                              If set to 0 no stabilizing energies are assigned to bases adjacent to
                                              helices in free ends and multiloops (so called dangling ends). Normally
                                              (dangles = 1) dangling end energies are assigned only to unpaired
                                              bases and a base cannot participate simultaneously in two dangling ends. In
                                              the partition function algorithm vrna_pf() these checks are neglected.
                                              To provide comparability between free energy minimization and partition function
                                              algorithms, the default setting is 2.
                                              This treatment of dangling ends gives more favorable energies to helices
                                              directly adjacent to one another, which can be beneficial since such
                                              helices often do engage in stabilizing interactions through co-axial
                                              stacking.\n
                                              If set to 3 co-axial stacking is explicitly included for
                                              adjacent helices in multiloops. The option affects only mfe folding
                                              and energy evaluation (vrna_mfe() and vrna_eval_structure()), as
                                              well as suboptimal folding (vrna_subopt()) via re-evaluation of energies.
                                              Co-axial stacking with one intervening mismatch is not considered so far.
                                              @note   Some function do not implement all dangle model but only a subset of
                                                      (0,1,2,3). In particular, partition function algorithms can only handle
                                                      0 and 2. Read the documentation of the particular recurrences or
                                                      energy evaluation function for information about the provided dangle
                                                      model.
                                        */
  int     special_hp;                   /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int     noLP;                         /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int     noGU;                         /**<  @brief  Do not allow GU pairs */
  int     noGUclosure;                  /**<  @brief  Do not allow loops to be closed by GU pair */
  int     logML;                        /**<  @brief  Use logarithmic scaling for multiloops */
  int     circ;                         /**<  @brief  Assume RNA to be circular instead of linear */
  int     gquad;                        /**<  @brief  Include G-quadruplexes in structure prediction */
  int     canonicalBPonly;              /**<  @brief  remove non-canonical bp's from constraint structures  */
  int     uniq_ML;                      /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
  int     energy_set;                   /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
  int     backtrack;                    /**<  @brief  Specifies whether or not secondary structures should be backtraced */
  char    backtrack_type;               /**<  @brief  Specifies in which matrix to backtrack */
  int     compute_bpp;                  /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
  char    nonstandards[64];             /**<  @brief  contains allowed non standard bases */
  int     max_bp_span;                  /**<  @brief  maximum allowed base pair span */

  int     min_loop_size;                /**<  @brief  Minimum size of hairpin loops
                                              @note The default value for this field is #TURN, however, it may
                                              be 0 in cofolding context.
                                        */
  int     window_size;                  /**<  @brief  Size of the sliding window for locally optimal structure prediction */
  int     oldAliEn;                     /**<  @brief  Use old alifold energy model */
  int     ribo;                         /**<  @brief  Use ribosum scoring table in alifold energy model */
  double  cv_fact;                      /**<  @brief  Co-variance scaling factor for consensus structure prediction */
  double  nc_fact;                      /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
  double  sfact;                        /**<  @brief  Scaling factor for partition function scaling */
  int     rtype[8];                     /**<  @brief  Reverse base pair type array */
  short   alias[MAXALPHA+1];            /**<  @brief  alias of an integer nucleotide representation */
  int     pair[MAXALPHA+1][MAXALPHA+1]; /**<  @brief  Integer representation of a base pair */
};


/**
 * @brief Apply default model details to a provided #vrna_md_t data structure
 *
 *  Use this function to initialize a #vrna_md_t data structure with
 *  its default values
 *
 *  @param md A pointer to the data structure that is about to be initialized
 */
void
vrna_md_set_default(vrna_md_t *md);

/**
 *  @brief Update the model details data structure
 *
 *  This function should be called after changing the vrna_md_t.energy_set attribute
 *  since it re-initializes base pairing related arrays within the #vrna_md_t data
 *  structure. In particular, #vrna_md_t.pair, #vrna_md_t.alias, and #vrna_md_t.rtype
 *  are set to the values that correspond to the specified #vrna_md_t.energy_set
 *  option
 *
 *  @see  #vrna_md_t, #vrna_md_t.energy_set, #vrna_md_t.pair, #vrna_md_t.rtype,
 *        #vrna_md_t.alias, vrna_md_set_default()
 */
void
vrna_md_update(vrna_md_t *md);

/**
 *  @brief Copy/Clone a #vrna_md_t model
 *
 *  Use this function to clone a given model either inplace (target container @p md_to
 *  given) or create a copy by cloning the source model and returning it (@p md_to == NULL).
 *
 *  @param md_to    The model to be overwritten (if non-NULL and @p md_to != @p md_from)
 *  @param md_from  The model to copy (if non-NULL)
 *  @return         A pointer to the copy model (or NULL if @p md_from == NULL)
 */
vrna_md_t *
vrna_md_copy( vrna_md_t       *md_to,
              const vrna_md_t *md_from);

/**
 *  @brief  Get a corresponding commandline parameter string of the options in a #vrna_md_t
 *
 *  @note This function is not threadsafe!
 */
char *
vrna_md_option_string(vrna_md_t  *md);

void
vrna_md_set_nonstandards(vrna_md_t *md, const char *ns_bases);

/**
 *  @brief  Reset the global default model details to a specific set of parameters, or their initial values
 *
 *  This function resets the global default model details to their initial values,
 *  i.e. as specified by the ViennaRNA Package release, upon passing NULL as argument.
 *  Alternatively it resets them according to a set of provided parameters.
 *
 *  @note The global default parameters affect all function calls of RNAlib where
 *        model details are not explicitly provided. Hence, any change of them
 *        is not considered threadsafe
 *  @warning  This function first resets the global default settings to factory
 *            defaults, and only then applies user provided settings (if any).
 *            User settings that do not meet specifications are skipped.
 *  @see  vrna_md_set_default(), #vrna_md_t
 *
 *  @param md_p A set of model details to use as global default (if NULL is passed, factory defaults are restored)
 */
void
vrna_md_defaults_reset(vrna_md_t *md_p);

/**
 *  @brief  Set default temperature for energy evaluation of loops
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_TEMPERATURE
 *  @param T  Temperature in centigrade
 */
void
vrna_md_defaults_temperature(double T);

/**
 *  @brief  Get default temperature for energy evaluation of loops
 *  @see vrna_md_defaults_temperature(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_TEMPERATURE
 *  @return  The global default settings for temperature in centigrade
 */
double
vrna_md_defaults_temperature_get(void);

/**
 *  @brief  Set default scaling factor of thermodynamic temperature in Boltzmann factors
 *
 *  Bolzmann factors are then computed as @f$ exp(-E / (b \cdot kT))@f$.
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_BETA_SCALE
 *  @param b  The scaling factor, default is 1.0
 */
void
vrna_md_defaults_betaScale(double b);

/**
 *  @brief  Get default scaling factor of thermodynamic temperature in Boltzmann factors
 *
 *  @see vrna_md_defaults_betaScale(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_BETA_SCALE
 *  @return  The global default thermodynamic temperature scaling factor
 */
double
vrna_md_defaults_betaScale_get(void);

/**
 *  @brief  Set default dangle model for structure prediction
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_DANGLES
 *  @param d  The dangle model
 */
void
vrna_md_defaults_dangles(int d);

/**
 *  @brief  Get default dangle model for structure prediction
 *  @see vrna_md_defaults_dangles(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_DANGLES
 *  @return The global default settings for the dangle model
 */
int
vrna_md_defaults_dangles_get(void);

/**
 *  @brief  Set default behavior for lookup of tabulated free energies for special hairpin loops, such as Tri-, Tetra-, or Hexa-loops.
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_SPECIAL_HP
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_special_hp(int flag);

/**
 *  @brief  Get default behavior for lookup of tabulated free energies for special hairpin loops, such as Tri-, Tetra-, or Hexa-loops.
 *  @see vrna_md_defaults_special_hp(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_SPECIAL_HP
 *  @return  The global default settings for the treatment of special hairpin loops
 */
int
vrna_md_defaults_special_hp_get(void);

/**
 *  @brief  Set default behavior for prediction of canonical secondary structures
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_NO_LP
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_noLP(int flag);

/**
 *  @brief  Get default behavior for prediction of canonical secondary structures
 *  @see vrna_md_defaults_noLP(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_NO_LP
 *  @return  The global default settings for predicting canonical secondary structures
 */
int
vrna_md_defaults_noLP_get(void);

/**
 *  @brief  Set default behavior for treatment of G-U wobble pairs
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_NO_GU
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_noGU(int flag);

/**
 *  @brief  Get default behavior for treatment of G-U wobble pairs
 *  @see vrna_md_defaults_noGU(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_NO_GU
 *  @return The global default settings for treatment of G-U wobble pairs
 */
int
vrna_md_defaults_noGU_get(void);

/**
 *  @brief  Set default behavior for G-U pairs as closing pair for loops
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_NO_GU_CLOSURE
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_noGUclosure(int flag);

/**
 *  @brief  Get default behavior for G-U pairs as closing pair for loops
 *  @see vrna_md_defaults_noGUclosure(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_NO_GU_CLOSURE
 *  @return The global default settings for treatment of G-U pairs closing a loop
 */
int
vrna_md_defaults_noGUclosure_get(void);

/**
 *  @brief  Set default behavior recomputing free energies of multi-branch loops using a logarithmic model
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_LOG_ML
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_logML(int flag);

/**
 *  @brief  Get default behavior recomputing free energies of multi-branch loops using a logarithmic model
 *  @see vrna_md_defaults_logML(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_LOG_ML
 *  @return The global default settings for logarithmic model in multi-branch loop free energy evaluation
 */
int
vrna_md_defaults_logML_get(void);

/**
 *  @brief  Set default behavior whether input sequences are circularized
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_CIRC
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_circ(int flag);

/**
 *  @brief  Get default behavior whether input sequences are circularized
 *  @see vrna_md_defaults_circ(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_CIRC
 *  @return The global default settings for treating input sequences as circular
 */
int
vrna_md_defaults_circ_get(void);

/**
 *  @brief  Set default behavior for treatment of G-Quadruplexes
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_GQUAD
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_gquad(int flag);

/**
 *  @brief  Get default behavior for treatment of G-Quadruplexes
 *  @see vrna_md_defaults_gquad(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_GQUAD
 *  @return The global default settings for treatment of G-Quadruplexes
 */
int
vrna_md_defaults_gquad_get(void);

/**
 *  @brief  Set default behavior for creating additional matrix for unique multi-branch loop prediction
 *  @note   Activating this option usually results in higher memory consumption!
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_UNIQ_ML
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_uniq_ML(int flag);

/**
 *  @brief  Get default behavior for creating additional matrix for unique multi-branch loop prediction
 *  @see vrna_md_defaults_uniq_ML(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_UNIQ_ML
 *  @return The global default settings for creating additional matrices for unique multi-branch loop prediction
 */
int
vrna_md_defaults_uniq_ML_get(void);

/**
 *  @brief  Set default energy set
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ENERGY_SET
 *  @param  e   Energy set (0, 1, 2, 3)
 */
void
vrna_md_defaults_energy_set(int e);

/**
 *  @brief  Get default energy set
 *  @see vrna_md_defaults_energy_set(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ENERGY_SET
 *  @return The global default settings for the energy set
 */
int
vrna_md_defaults_energy_set_get(void);

/**
 *  @brief  Set default behavior for whether to backtrack secondary structures
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_BACKTRACK
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_backtrack(int flag);

/**
 *  @brief  Get default behavior for whether to backtrack secondary structures
 *  @see vrna_md_defaults_backtrack(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_BACKTRACK
 *  @return The global default settings for backtracking structures
 */
int
vrna_md_defaults_backtrack_get(void);

/**
 *  @brief  Set default backtrack type, i.e. which DP matrix is used
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_BACKTRACK_TYPE
 *  @param  t   The type ('F', 'C', or 'M')
 */
void
vrna_md_defaults_backtrack_type(char t);

/**
 *  @brief  Get default backtrack type, i.e. which DP matrix is used
 *  @see vrna_md_defaults_backtrack_type(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_BACKTRACK_TYPE
 *  @return The global default settings that specify which DP matrix is used for backtracking
 */
char
vrna_md_defaults_backtrack_type_get(void);

/**
 *  @brief  Set the default behavior for whether to compute base pair probabilities after partition function computation
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_COMPUTE_BPP
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_compute_bpp(int flag);

/**
 *  @brief  Get the default behavior for whether to compute base pair probabilities after partition function computation
 *  @see vrna_md_defaults_compute_bpp(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_COMPUTE_BPP
 *  @return The global default settings that specify whether base pair probabilities are computed together with partition function
 */
int
vrna_md_defaults_compute_bpp_get(void);

/**
 *  @brief  Set default maximal base pair span
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_MAX_BP_SPAN
 *  @param  span  Maximal base pair span
 */
void
vrna_md_defaults_max_bp_span(int span);

/**
 *  @brief  Get default maximal base pair span
 *  @see vrna_md_defaults_max_bp_span(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_MAX_BP_SPAN
 *  @return The global default settings for maximum base pair span
 */
int
vrna_md_defaults_max_bp_span_get(void);

/**
 *  @brief  Set default minimal loop size
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #TURN
 *  @param  size  Minimal size, i.e. number of unpaired nucleotides for a hairpin loop
 */
void
vrna_md_defaults_min_loop_size(int size);

/**
 *  @brief  Get default minimal loop size
 *  @see vrna_md_defaults_min_loop_size(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #TURN
 *  @return The global default settings for minimal size of hairpin loops
 */
int
vrna_md_defaults_min_loop_size_get(void);

/**
 *  @brief  Set default window size for sliding window structure prediction approaches
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_WINDOW_SIZE
 *  @param  size  The size of the sliding window
 */
void
vrna_md_defaults_window_size(int size);

/**
 *  @brief  Get default window size for sliding window structure prediction approaches
 *  @see vrna_md_defaults_window_size(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_WINDOW_SIZE
 *  @return The global default settings for the size of the sliding window
 */
int
vrna_md_defaults_window_size_get(void);

/**
 *  @brief  Set default behavior for whether to use old energy model for comparative structure prediction
 *  @note   This option is outdated. Activating the old energy model usually results in worse consensus
 *          structure predictions.
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_OLD_EN
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_oldAliEn(int flag);

/**
 *  @brief  Get default behavior for whether to use old energy model for comparative structure prediction
 *  @see vrna_md_defaults_oldAliEn(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_OLD_EN
 *  @return The global default settings for using old energy model for comparative structure prediction
 */
int
vrna_md_defaults_oldAliEn_get(void);

/**
 *  @brief  Set default behavior for whether to use Ribosum Scoring in comparative structure prediction
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_RIBO
 *  @param  flag  On/Off switch (0 = OFF, else = ON)
 */
void
vrna_md_defaults_ribo(int flag);

/**
 *  @brief  Get default behavior for whether to use Ribosum Scoring in comparative structure prediction
 *  @see vrna_md_defaults_ribo(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_RIBO
 *  @return The global default settings for using Ribosum scoring in comparative structure prediction
 */
int
vrna_md_defaults_ribo_get(void);

/**
 *  @brief  Set the default co-variance scaling factor used in comparative structure prediction
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_CV_FACT
 *  @param  factor  The co-variance factor
 */
void
vrna_md_defaults_cv_fact(double factor);

/**
 *  @brief  Get the default co-variance scaling factor used in comparative structure prediction
 *  @see vrna_md_defaults_cv_fact(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_CV_FACT
 *  @return The global default settings for the co-variance factor
 */
double
vrna_md_defaults_cv_fact_get(void);

/**
 *  @brief
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_NC_FACT
 *  @param factor
 */
void
vrna_md_defaults_nc_fact(double factor);

/**
 *  @brief
 *  @see vrna_md_defaults_nc_fact(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t, #VRNA_MODEL_DEFAULT_ALI_NC_FACT
 *  @return
 */
double
vrna_md_defaults_nc_fact_get(void);

/**
 *  @brief  Set the default scaling factor used to avoid under-/overflows in partition function computation
 *  @see vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t
 *  @param  factor  The scaling factor  (default: 1.07)
 */
void
vrna_md_defaults_sfact(double factor);

/**
 *  @brief  Get the default scaling factor used to avoid under-/overflows in partition function computation
 *  @see vrna_md_defaults_sfact(), vrna_md_defaults_reset(), vrna_md_set_default(), #vrna_md_t
 *  @return The global default settings of the scaling factor
 */
double
vrna_md_defaults_sfact_get(void);

#ifdef  VRNA_BACKWARD_COMPAT

#define model_detailsT        vrna_md_t               /* restore compatibility of struct rename */

/* BEGIN deprecated global variables: */

/**
 *  @brief Rescale energy parameters to a temperature in degC.
 * 
 *  Default is 37C. You have to call the update_..._params() functions after
 *  changing this parameter.
 *  @deprecated   Use vrna_md_defaults_temperature(), and vrna_md_defaults_temperature_get()
 *                to change, and read the global default temperature settings
 *  @see vrna_md_defaults_temperature(), vrna_md_defaults_temperature_get(), vrna_md_defaults_reset()
 */
extern double temperature;

/**
 *  @brief A scaling factor used by pf_fold() to avoid overflows.
 * 
 *  Should be set to approximately @f$exp{((-F/kT)/length)}@f$, where @f$F@f$ is an estimate
 *  for the ensemble free energy, for example the minimum free energy. You must
 *  call update_pf_params() after changing this parameter.\n
 *  If pf_scale is -1 (the default) , an estimate will be provided
 *  automatically when computing partition functions, e.g. pf_fold()
 *  The automatic estimate is usually insufficient for sequences more
 *  than a few hundred bases long.
 */
extern double pf_scale;

/**
 *  @brief Switch the energy model for dangling end contributions (0, 1, 2, 3)
 * 
 *  If set to 0 no stabilizing energies are assigned to bases adjacent to
 *  helices in free ends and multiloops (so called dangling ends). Normally
 *  (dangles = 1) dangling end energies are assigned only to unpaired
 *  bases and a base cannot participate simultaneously in two dangling ends. In
 *  the partition function algorithm pf_fold() these checks are neglected.
 *  If #dangles is set to 2, all folding routines will follow this convention.
 *  This treatment of dangling ends gives more favorable energies to helices
 *  directly adjacent to one another, which can be beneficial since such
 *  helices often do engage in stabilizing interactions through co-axial
 *  stacking.\n
 *  If dangles = 3 co-axial stacking is explicitly included for
 *  adjacent helices in multiloops. The option affects only mfe folding
 *  and energy evaluation (fold() and energy_of_structure()), as
 *  well as suboptimal folding (subopt()) via re-evaluation of energies.
 *  Co-axial stacking with one intervening mismatch is not considered so far.
 * 
 *  Default is 2 in most algorithms, partition function algorithms can only handle 0 and 2
 */
extern int  dangles;

/**
 *  @brief Include special stabilizing energies for some tri-, tetra- and hexa-loops;
 * 
 *  default is 1.
 */
extern int  tetra_loop;

/**
 *  @brief Global switch to avoid/allow helices of length 1
 * 
 *  Disallow all pairs which can only occur as lonely pairs (i.e. as helix
 *  of length 1). This avoids lonely base pairs in the predicted structures in
 *  most cases.
 */
extern int    noLonelyPairs;

/**
 *  @brief Global switch to forbid/allow GU base pairs at all
 */
extern int  noGU;

/**
 *  @brief GU allowed only inside stacks if set to 1
 */
extern int  no_closingGU;

/**
 *  @brief backward compatibility variable.. this does not effect anything
 */
extern int  circ;

/**
 *  @brief Allow G-quadruplex formation
 */
extern int gquad;

/**
 *  Do not use this variable, it will eventually be removed in one of the next versions
 */
extern int canonicalBPonly;

/**
 *  @brief do ML decomposition uniquely (for subopt)
 */
extern  int uniq_ML;

/**
 *  @brief 0 = BP; 1=any with GC; 2=any with AU-parameter
 * 
 *  If set to 1 or 2: fold sequences from an artificial alphabet ABCD..., where A
 *  pairs B, C pairs D, etc. using either GC (1) or AU parameters (2);
 *  default is 0, you probably don't want to change it.
 */
extern int  energy_set;

/**
 *  @brief do backtracking, i.e. compute secondary structures or base pair probabilities
 * 
 *  If 0, do not calculate pair probabilities in pf_fold(); this is about
 *  twice as fast. Default is 1.
 */
extern int    do_backtrack;

/**
 *  @brief A backtrack array marker for inverse_fold()
 * 
 *  If set to 'C': force (1,N) to be paired,
 *  'M' fold as if the sequence were inside a multiloop. Otherwise ('F') the
 *  usual mfe structure is computed.
 */
extern char backtrack_type;

/**
 *  @brief contains allowed non standard base pairs
 * 
 *  Lists additional base pairs that will be allowed to form in addition to
 *  GC, CG, AU, UA, GU and UG. Nonstandard base pairs are given a stacking
 *  energy of 0.
 */
extern char *nonstandards;

/**
 *  @brief Maximum allowed base pair span
 *
 *  A value of -1 indicates no restriction for distant base pairs.
 */
extern int max_bp_span;

/**
 *  @brief use old alifold energies (with gaps)
 */
extern int oldAliEn;

/**
 *  @brief use ribosum matrices
 */
extern int ribo;            

extern double cv_fact;

extern double nc_fact;

/** @brief if nonzero use logarithmic ML energy in energy_of_struct  */
extern  int logML;

/* END deprecated global variables: */

/**
 * @brief Set default model details
 *
 *  Use this function if you wish to initialize a #vrna_md_t data structure with
 *  its default values, i.e. the global model settings as provided by the deprecated
 *  global variables.
 *
 *  @deprecated This function will vanish as soon as backward compatibility of
 *              RNAlib is dropped (expected in version 3).
 *              Use vrna_md_set_default() instead!
 *
 *  @param md A pointer to the data structure that is about to be initialized
 */
void
set_model_details(vrna_md_t *md);

char *
option_string(void);

#endif
/**
 * @}
 */

#endif
