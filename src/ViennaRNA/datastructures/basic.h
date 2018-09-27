#ifndef VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H
#define VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H

/**
 *  @file     ViennaRNA/datastructures/basic.h
 *  @ingroup  data_structures
 *  @brief    Various data structures and pre-processor macros
 */

/**
 *  @addtogroup   data_structures
 *  @{
 *
 *  @brief All datastructures and typedefs shared among the ViennaRNA Package can be found here
 *
 */

/* below are several convenience typedef's we use throughout the ViennaRNA library */

/** @brief Typename for the base pair repesenting data structure #vrna_basepair_s */
typedef struct vrna_basepair_s vrna_basepair_t;

/** @brief Typename for the base pair list repesenting data structure #vrna_elem_prob_s */
typedef struct vrna_elem_prob_s vrna_plist_t;

/** @brief Typename for the base pair stack repesenting data structure #vrna_bp_stack_s */
typedef struct vrna_bp_stack_s vrna_bp_stack_t;

/** @brief Typename for data structure #vrna_cpair_s */
typedef struct vrna_cpair_s vrna_cpair_t;

/** @brief Typename for stack of partial structures #vrna_sect_s */
typedef struct vrna_sect_s vrna_sect_t;

typedef struct vrna_data_linear_s vrna_data_lin_t;

typedef struct vrna_color_s vrna_color_t;

/** @brief Typename for floating point number in partition function computations */
#ifdef  USE_FLOAT_PF
typedef float FLT_OR_DBL;
#else
typedef double FLT_OR_DBL;
#endif


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* the following typedefs are for backward compatibility only */

/**
 *  @brief Old typename of #vrna_basepair_s
 *  @deprecated Use #vrna_basepair_t instead!
 */
typedef struct vrna_basepair_s PAIR;

/**
 *  @brief Old typename of #vrna_elem_prob_s
 *  @deprecated Use #vrna_ep_t or #vrna_elem_prob_s instead!
 */
typedef struct vrna_elem_prob_s plist;
/**
 *  @brief Old typename of #vrna_cpair_s
 *  @deprecated Use #vrna_cpair_t instead!
 */
typedef struct vrna_cpair_s cpair;

/**
 *  @brief Old typename of #vrna_sect_s
 *  @deprecated Use #vrna_sect_t instead!
 */
typedef struct vrna_sect_s sect;

/**
 *  @brief Old typename of #vrna_bp_stack_s
 *  @deprecated Use #vrna_bp_stack_t instead!
 */
typedef struct vrna_bp_stack_s bondT;

#endif

#include <ViennaRNA/params/constants.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/dp_matrices.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/grammar.h>
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/utils/structures.h"

/*
 * ############################################################
 * Here are the type definitions of various datastructures
 * shared among the Vienna RNA Package
 * ############################################################
 */

/**
 *  @brief  Base pair data structure used in subopt.c
 */
struct vrna_basepair_s {
  int i;
  int j;
};

/**
 *  @brief this datastructure is used as input parameter in functions of PS_dot.c
 */
struct vrna_cpair_s {
  int   i, j, mfe;
  float p, hue, sat;
  int   type;
};

struct vrna_color_s {
  float hue;
  float sat;
  float bri;
};

struct vrna_data_linear_s {
  unsigned int  position;
  float         value;
  vrna_color_t  color;
};


/**
 *  @brief  Stack of partial structures for backtracking
 */
struct vrna_sect_s {
  int i;
  int j;
  int ml;
};

/**
 *  @brief  Base pair stack element
 */
struct vrna_bp_stack_s {
  unsigned int  i;
  unsigned int  j;
};


/*
 * ############################################################
 * RNAup data structures
 * ############################################################
 */

/**
 *  @brief contributions to p_u
 */
typedef struct pu_contrib {
  double  **H;    /**<  @brief  hairpin loops */
  double  **I;    /**<  @brief  interior loops */
  double  **M;    /**<  @brief  multi loops */
  double  **E;    /**<  @brief  exterior loop */
  int     length; /**<  @brief  length of the input sequence */
  int     w;      /**<  @brief  longest unpaired region */
} pu_contrib;

/**
 *  @brief  interaction data structure for RNAup
 */
typedef struct interact {
  double  *Pi;      /**<  @brief  probabilities of interaction */
  double  *Gi;      /**<  @brief  free energies of interaction */
  double  Gikjl;    /**<  @brief  full free energy for interaction between [k,i] k<i
                     *            in longer seq and [j,l] j<l in shorter seq */
  double  Gikjl_wo; /**<  @brief  Gikjl without contributions for prob_unpaired */
  int     i;        /**<  @brief  k<i in longer seq */
  int     k;        /**<  @brief  k<i in longer seq */
  int     j;        /**<  @brief  j<l in shorter seq */
  int     l;        /**<  @brief  j<l in shorter seq */
  int     length;   /**<  @brief  length of longer sequence */
} interact;

/**
 *  @brief  Collection of all free_energy of beeing unpaired values for output
 */
typedef struct pu_out {
  int     len;        /**<  @brief  sequence length */
  int     u_vals;     /**<  @brief  number of different -u values */
  int     contribs;   /**<  @brief  [-c "SHIME"] */
  char    **header;   /**<  @brief  header line */
  double  **u_values; /**<  @brief  (the -u values * [-c "SHIME"]) * seq len */
} pu_out;

/**
 *  @brief  constraints for cofolding
 */
typedef struct constrain {
  int   *indx;
  char  *ptype;
} constrain;

/*
 * ############################################################
 * RNAduplex data structures
 * ############################################################
 */

/**
 *  @brief  Data structure for RNAduplex
 */
typedef struct {
  int     i;
  int     j;
  int     end;
  char    *structure;
  double  energy;
  double  energy_backtrack;
  double  opening_backtrack_x;
  double  opening_backtrack_y;
  int     offset;
  double  dG1;
  double  dG2;
  double  ddG;
  int     tb;
  int     te;
  int     qb;
  int     qe;
} duplexT;

/*
 * ############################################################
 * RNAsnoop data structures
 * ############################################################
 */

/**
 *  @brief  Data structure for RNAsnoop (fold energy list)
 */
typedef struct node {
  int         k;
  int         energy;
  struct node *next;
} folden;

/**
 *  @brief  Data structure for RNAsnoop
 */
typedef struct {
  int   i;
  int   j;
  int   u;
  char  *structure;
  float energy;
  float Duplex_El;
  float Duplex_Er;
  float Loop_E;
  float Loop_D;
  float pscd;
  float psct;
  float pscg;
  float Duplex_Ol;
  float Duplex_Or;
  float Duplex_Ot;
  float fullStemEnergy;
} snoopT;


/*
 * ############################################################
 * PKplex data structures
 * ############################################################
 */

/**
 *  @brief  Data structure used in RNApkplex
 */
typedef struct dupVar {
  int     i;
  int     j;
  int     end;
  char    *pk_helix;
  char    *structure;
  double  energy;
  int     offset;
  double  dG1;
  double  dG2;
  double  ddG;
  int     tb;
  int     te;
  int     qb;
  int     qe;
  int     inactive;
  int     processed;
} dupVar;

/**
 *  @brief  Dummy symbol to check whether the library was build using C11/C++11 features
 *
 *  By default, several data structures of our new v3.0 API use C11/C++11 features, such
 *  as unnamed unions, unnamed structs. However, these features can be deactivated at
 *  compile time to allow building the library and executables with compilers that do not
 *  support these features.
 *
 *  Now, the problem arises that once our static library is compiled and a third-party
 *  application is supposed to link against it, it needs to know, at compile time, how to
 *  correctly address particular data structures. This is usually implicitely taken care of
 *  through the API exposed in our header files. Unfortunately, we had some preprocessor directives
 *  in our header files that changed the API depending on the capabilities of the compiler
 *  the third-party application is build with. This in turn prohibited the use of an RNAlib
 *  compiled without C11/C++11 support in a program that compiles/links with enabled C11/C++11
 *  support and vice-versa.
 *
 *  Therefore, we introduce this dummy symbol which can be used to check, whether the
 *  static library was build with C11/C++11 features.
 *
 *  @note If the symbol is present, the library was build with enabled C11/C++11 features support
 *  and no action is required. However, if the symbol is missing in RNAlib >= 2.2.9, programs
 *  that link to RNAlib must define a pre-processor identifier @em VRNA_DISABLE_C11_FEATURES before
 *  including any ViennaRNA Package header file, for instance by adding a @em CPPFLAG
 *  @code
 * CPPFLAGS+=-DVRNA_DISABLE_C11_FEATURES
 *  @endcode
 *
 *  @since v2.2.9
 */
#ifndef VRNA_DISABLE_C11_FEATURES
void vrna_C11_features(void);


#endif

/**
 * @}
 */


#endif
