#ifndef __VIENNA_RNA_PACKAGE_MODEL_H__
#define __VIENNA_RNA_PACKAGE_MODEL_H__

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif


/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  \addtogroup   basic_data_structures
 *
 *  @{
 *
 *  \file model.h
 *  \brief The model details data structure and its corresponding modifiers
 */


#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0
#define VRNA_MODEL_DEFAULT_PF_SCALE       -1
#define VRNA_MODEL_DEFAULT_BETA_SCALE     1.
#define VRNA_MODEL_DEFAULT_DANGLES        2
#define VRNA_MODEL_DEFAULT_SPECIAL_HP     1
#define VRNA_MODEL_DEFAULT_NO_LP          0
#define VRNA_MODEL_DEFAULT_NO_GU          0
#define VRNA_MODEL_DEFAULT_NO_GU_CLOSURE  0
#define VRNA_MODEL_DEFAULT_CIRC           0
#define VRNA_MODEL_DEFAULT_GQUAD          0
#define VRNA_MODEL_DEFAULT_CANONICAL_BP   0
#define VRNA_MODEL_DEFAULT_UNIQ_ML        0
#define VRNA_MODEL_DEFAULT_ENERGY_SET     0
#define VRNA_MODEL_DEFAULT_BACKTRACK      1
#define VRNA_MODEL_DEFAULT_BACKTRACK_TYPE 'F'
#define VRNA_MODEL_DEFAULT_COMPUTE_BPP    1
#define VRNA_MODEL_DEFAULT_MAX_BP_SPAN    -1
#define VRNA_MODEL_DEFAULT_LOG_ML         0
#define VRNA_MODEL_DEFAULT_ALI_OLD_EN     0
#define VRNA_MODEL_DEFAULT_ALI_RIBO       0
#define VRNA_MODEL_DEFAULT_ALI_CV_FACT    1.
#define VRNA_MODEL_DEFAULT_ALI_NC_FACT    1.


#ifdef  VRNA_BACKWARD_COMPAT

#ifndef MAXALPHA
/**
 *  \brief Maximal length of alphabet
 */
#define MAXALPHA              20
#endif

#define model_detailsT        vrna_md_t               /* restore compatibility of struct rename */
#define set_model_details(a)  vrna_md_set_globals(a)  /* restore compatibility of function rename */

/* BEGIN deprecated global variables: */

/**
 *  \brief Rescale energy parameters to a temperature in degC.
 * 
 *  Default is 37C. You have to call the update_..._params() functions after
 *  changing this parameter.
 */
extern double temperature;

/**
 *  \brief A scaling factor used by pf_fold() to avoid overflows.
 * 
 *  Should be set to approximately \f$exp{((-F/kT)/length)}\f$, where \f$F\f$ is an estimate
 *  for the ensemble free energy, for example the minimum free energy. You must
 *  call update_pf_params() after changing this parameter.\n
 *  If pf_scale is -1 (the default) , an estimate will be provided
 *  automatically when computing partition functions, e.g. pf_fold()
 *  The automatic estimate is usually insufficient for sequences more
 *  than a few hundred bases long.
 */
extern double pf_scale;

/**
 *  \brief Switch the energy model for dangling end contributions (0, 1, 2, 3)
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
 *  adjacent helices in mutli-loops. The option affects only mfe folding
 *  and energy evaluation (fold() and energy_of_structure()), as
 *  well as suboptimal folding (subopt()) via re-evaluation of energies.
 *  Co-axial stacking with one intervening mismatch is not considered so far.
 * 
 *  Default is 2 in most algorithms, partition function algorithms can only handle 0 and 2
 */
extern int  dangles;

/**
 *  \brief Include special stabilizing energies for some tri-, tetra- and hexa-loops;
 * 
 *  default is 1.
 */
extern int  tetra_loop;

/**
 *  \brief Global switch to avoid/allow helices of length 1
 * 
 *  Disallow all pairs which can only occur as lonely pairs (i.e. as helix
 *  of length 1). This avoids lonely base pairs in the predicted structures in
 *  most cases.
 */
extern int    noLonelyPairs;

/**
 *  \brief Global switch to forbid/allow GU base pairs at all
 */
extern int  noGU;

/**
 *  \brief GU allowed only inside stacks if set to 1
 */
extern int  no_closingGU;

/**
 *  \brief backward compatibility variable.. this does not effect anything
 */
extern int  circ;

/**
 *  \brief Allow G-quadruplex formation
 */
extern int gquad;

/**
 *  Do not use this variable, it will eventually be removed in one of the next versions
 */
extern int canonicalBPonly;

/**
 *  \brief do ML decomposition uniquely (for subopt)
 */
extern  int uniq_ML;

/**
 *  \brief 0 = BP; 1=any mit GC; 2=any mit AU-parameter
 * 
 *  If set to 1 or 2: fold sequences from an artificial alphabet ABCD..., where A
 *  pairs B, C pairs D, etc. using either GC (1) or AU parameters (2);
 *  default is 0, you probably don't want to change it.
 */
extern int  energy_set;

/**
 *  \brief do backtracking, i.e. compute secondary structures or base pair probabilities
 * 
 *  If 0, do not calculate pair probabilities in pf_fold(); this is about
 *  twice as fast. Default is 1.
 */
extern int    do_backtrack;

/**
 *  \brief A backtrack array marker for inverse_fold()
 * 
 *  If set to 'C': force (1,N) to be paired,
 *  'M' fold as if the sequence were inside a multi-loop. Otherwise ('F') the
 *  usual mfe structure is computed.
 */
extern char backtrack_type;

/**
 *  \brief contains allowed non standard base pairs
 * 
 *  Lists additional base pairs that will be allowed to form in addition to
 *  GC, CG, AU, UA, GU and UG. Nonstandard base pairs are given a stacking
 *  energy of 0.
 */
extern char *nonstandards;

/**
 *  \brief Maximum allowed base pair span
 *
 *  A value of -1 indicates no restriction for distant base pairs.
 */
extern int max_bp_span;

/**
 *  \brief use old alifold energies (with gaps)
 */
extern int oldAliEn;

/**
 *  \brief use ribosum matrices
 */
extern int ribo;            

extern double cv_fact;

extern double nc_fact;

/** \brief if nonzero use logarithmic ML energy in energy_of_struct  */
extern  int logML;

/* END deprecated global variables: */

#endif

/**
 *  \brief The data structure that contains the complete model details used throughout the calculations
 *
 */
typedef struct vrna_md_t{
  double  temperature;      /**<  \brief  The temperature used to scale the thermodynamic parameters */
  double  betaScale;        /**<  \brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
  int     dangles;          /**<  \brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
                                  \note   Some function do not implement all dangle model but only a subset of
                                          (0,1,2,3). Read the documentaion of the particular recurrences or
                                          energy evaluation function for information about the provided dangle
                                          model.
                            */
  int     special_hp;       /**<  \brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int     noLP;             /**<  \brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int     noGU;             /**<  \brief  Do not allow GU pairs */
  int     noGUclosure;      /**<  \brief  Do not allow loops to be closed by GU pair */
  int     logML;            /**<  \brief  Use logarithmic scaling for multi loops */
  int     circ;             /**<  \brief  Assume RNA to be circular instead of linear */
  int     gquad;            /**<  \brief  Include G-quadruplexes in structure prediction */
  int     canonicalBPonly;  /**<  \brief  remove non-canonical bp's from constraint structures  */
  int     uniq_ML;          /**<  \brief  Flag to ensure unique multibranch loop decomposition during folding */
  int     energy_set;       /**<  \brief  Specifies the energy set that defines set of compatible base pairs */
  int     backtrack;        /**<  \brief  Specifies whether or not secondary structures should be backtraced */
  char    backtrack_type;   /**<  \brief  Specifies in which matrix to backtrack */
  int     compute_bpp;      /**<  \brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
  char    nonstandards[33]; /**<  \brief  contains allowed non standard bases */
  int     max_bp_span;      /**<  \brief  maximum allowed base pair span */

  int     min_loop_size;    /**<  \brief  Minimum size of hairpin loops
                              
                                  \note The default value for this field is #TURN, however, it may
                                  be 0 in cofolding context.
                            */

  int     oldAliEn;         /**<  \brief  Use old alifold energy model */
  int     ribo;             /**<  \brief  Use ribosum scoring table in alifold energy model */
  double  cv_fact;          /**<  \brief  Covariance scaling factor for consensus structure prediction */
  double  nc_fact;
  double  sfact;            /**<  \brief  Scaling factor for partition function scaling */
  int     rtype[8];
  short   alias[MAXALPHA+1];
  int     pair[MAXALPHA+1][MAXALPHA+1];
} vrna_md_t;


/**
 * \brief Set default model details
 *
 *  Use this function if you wish to initialize a #vrna_md_t data structure with
 *  its default values, i.e. the global model settings
 *
 *  \see
 *
 *  \param md A pointer to the data structure that is about to be initialized
 */
void vrna_md_set_default(vrna_md_t *md);


void vrna_md_set_nonstandards(vrna_md_t *md, const char *ns);


void vrna_md_set_dangles(vrna_md_t *md, int d);


int vrna_md_get_dangles(vrna_md_t *md);


void vrna_md_set_temperature(vrna_md_t *md, double T);


double vrna_md_get_temperature(vrna_md_t *md);


void vrna_md_set_special_hp(vrna_md_t *md, int shp);


int vrna_md_get_special_hp(vrna_md_t *md);


void vrna_md_set_gquad(vrna_md_t *md, int g);


int vrna_md_get_gquad(vrna_md_t *md);


void vrna_md_set_nolp(vrna_md_t *md, int nolp);


int vrna_md_get_nolp(vrna_md_t *md);


void vrna_md_set_betascale(vrna_md_t *md, double b);


double vrna_md_get_betascale(vrna_md_t *md);

/**
 *  \brief Update the model details
 */
void vrna_md_update(vrna_md_t *md);

/**
 * \brief Set default model details
 *
 *  Use this function if you wish to initialize a #vrna_md_t data structure with
 *  its default values, i.e. the global model settings as provided by the deprecated
 *  global variables.
 *
 *  \deprecated This function will vanish as soon as backward compatibility of
 *              RNAlib is dropped (expected in version 3).
 *              Use vrna_md_set_default instead!
 *
 *  \param md A pointer to the data structure that is about to be initialized
 */
void vrna_md_set_globals(vrna_md_t *md);

/**
 * @}
 */

#endif
