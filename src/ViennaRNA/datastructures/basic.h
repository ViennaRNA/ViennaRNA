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
 */

/* below are several convenience typedef's we use throughout the ViennaRNA library */

/** @brief Typename for the base pair list repesenting data structure #vrna_elem_prob_s */
typedef struct vrna_elem_prob_s vrna_plist_t;

/** @brief Typename for data structure #vrna_cpair_s */
typedef struct vrna_cpair_s vrna_cpair_t;

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

/** @brief Typename for base pair element
 *  @deprecated Use vrna_bp_t instead!
 */
typedef struct {
  int i;
  int j;
} vrna_basepair_t;


/** @brief Typename for the base pair stack element */
typedef struct {
  unsigned int  i;
  unsigned int  j;
} vrna_bp_stack_t;


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
typedef vrna_bp_stack_t bondT;

#endif

/**
 *  @}
 */


/**
 *  @addtogroup data_structures_bt
 *  @{
 */

/**
 *  @brief  The backtrack stack data structure
 *
 *  @see  vrna_bts_init(), vrna_bts_free(), vrna_bts_push(),
 *        vrna_bts_top(), vrna_bts_pop()
 */
typedef struct vrna_bt_stack_s  *vrna_bts_t;


/**
 *  @brief  The basepair stack data structure
 *
 *  @see  vrna_bps_init(), vrna_bps_free(), vrna_bps_push(),
 *        vrna_bps_top(), vrna_bps_pop(), vrna_bps_at()
 */
typedef struct vrna_bp_stack_s  *vrna_bps_t;


/** @brief Typename for the base pair repesenting data structure #vrna_basepair_s */
typedef struct vrna_basepair_s  vrna_bp_t;

/** @brief Typename for stack of partial structures #vrna_sect_s */
typedef struct vrna_sect_s vrna_sect_t;

/**
 *  @}
 */

#include <ViennaRNA/params/constants.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/datastructures/dp_matrices.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/structures/problist.h"
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/structures/dotbracket.h"
#include "ViennaRNA/structures/utils.h"
#include "ViennaRNA/structures/metrics.h"
#include "ViennaRNA/structures/helix.h"

/*
 * ############################################################
 * Here are the type definitions of various datastructures
 * shared among the Vienna RNA Package
 * ############################################################
 */

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
 *  @addtogroup   data_structures
 *  @{
 */

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
  double  **I;    /**<  @brief  internal loops */
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
 *  @}
 */



/**
 *  @addtogroup data_structures_bt
 *  @{
 */


/**
 *  @brief  Stack of partial structures for backtracking
 */
struct vrna_sect_s {
  int           i;
  int           j;
  unsigned int  ml;
};


/**
 *  @brief  Base pair data structure used in subopt.c
 */
struct vrna_basepair_s {
  unsigned int i;
  unsigned int j;
  unsigned int L;
  unsigned int l[3];
};

/**
 *  @brief  Get an initialized backtrack stack
 *
 *  This function yields an initialized backtracking stack
 *  that holds all elements that need to be further evaluated.
 *  The individual elements stored in the stack are of type
 *  #vrna_sect_t and store the sequence delimiters and corresponding
 *  backtrack DP matrix flag.
 *
 *  @note Memory for the stack must be released via the vrna_bts_free()
 *        function.
 *
 *  @see  #vrna_bts_t, vrna_bts_free(), vrna_bts_push(),
 *        vrna_bts_top(),vrna_bts_pop(), vrna_bts_size(),
 *
 *  @param  n   The initial size of the backtrack stack
 *  @return     An initialized backtrack stack, ready for usage in backtracking functions
 */
vrna_bts_t
vrna_bts_init(size_t  n);

/**
 *  @brief  Release memory occupied by a backtrack stack
 *
 *  @param  bts   The backtrack stack that should be free'd
 */
void
vrna_bts_free(vrna_bts_t bts);


/**
 *  @brief  Push a new interval onto the backtrack stack
 *
 *  This function pushes a new sequence interval for backtracking
 *  onto the backstracking stack @p bts.
 *
 *  @param  bts     The backtrack stack
 *  @param  element The sequence interval and corresponding DP matrix flag
 *  @return         The size of the backtrack stack after pushing the new interval
 */
size_t
vrna_bts_push(vrna_bts_t  bts,
              vrna_sect_t element);


/**
 *  @brief  Retrieve the top element of the backtrack stack
 *
 *  Retrieves the last element put onto the stack, or a zero'd
 *  out #vrna_sect_t structure. The latter is returned on error
 *  or when topping an empty stack.
 *
 *  @param  bts   The backtrack stack
 *  @return       The top element of the backtrack stack, or a zero'd out #vrna_sect_t
 */
vrna_sect_t
vrna_bts_top(vrna_bts_t bts);


/**
 *  @brief  Pop last element of from backtrack stack
 *
 *  Retrieves and removes the last element put onto the stack, or a zero'd
 *  out #vrna_sect_t structure. The latter is returned on error
 *  or when topping an empty stack.
 *
 *  @param  bts   The backtrack stack
 *  @return       The top element of the backtrack stack, or a zero'd out #vrna_sect_t
 */
vrna_sect_t
vrna_bts_pop(vrna_bts_t bts);

/**
 *  @brief  Get the size of the backtrack stack
 *
 *  @param  bts   The backtrack stack
 *  @return       The size of the backtracking stack
 */
size_t
vrna_bts_size(vrna_bts_t bts);


/**
 *  @brief  Get an initialized base pair stack
 *
 *  Base pair stacks are used in the backtracking procedure to store
 *  all base pairs and structural elements that have been identified
 *  so far. Thos function returns an initialized backtracking stack
 *  with initial size @p n. Individual elements stored in this stack
 *  are of type #vrna_bp_t.
 *
 *  @note Memory for the stack must be released via the vrna_bps_free()
 *        function.
 *
 *  @param  n   The initial size of the base pair stack
 *  @return     An initialized base pair stack
 */
vrna_bps_t
vrna_bps_init(size_t  n);


/**
 *  @brief  Release memory of a base pair stack
 *
 *  @param  bps   The base pair stack to be free'd
 */
void
vrna_bps_free(vrna_bps_t bps);


/**
 *  @brief  Put a new base pair element on top of the stack
 *
 *  @param  bps   The base pair stack
 *  @param  pair  The base pair to be put onto the stack
 *  @return       The size of the base pair stack after pushing the base pair
 */
size_t
vrna_bps_push(vrna_bps_t  bps,
              vrna_bp_t   pair);


/**
 *  @brief  Retrieve the top element of the base pair stack
 *
 *  Retrieves the last element put onto the stack, or a zero'd
 *  out #vrna_bp_t structure. The latter is returned on error
 *  or when topping an empty stack.
 *
 *  @param  bps   The base pair stack
 *  @return       The top element of the base pair stack, or a zero'd out #vrna_bp_t
 */
vrna_bp_t
vrna_bps_top(vrna_bps_t bps);


/**
 *  @brief  Pop last element of from base pair stack
 *
 *  Retrieves and removes the last element put onto the stack, or a zero'd
 *  out #vrna_bp_t structure. The latter is returned on error or when
 *  topping an empty stack.
 *
 *  @param  bps   The base pair stack
 *  @return       The top element of the backtrack stack, or a zero'd out #vrna_bp_t
 */
vrna_bp_t
vrna_bps_pop(vrna_bps_t bps);


/**
 *  @brief  Retrieve the n'th element of the base pair stack
 *
 *  Retrieves the n'th element counted from the bottom of the stack (0-based),
 *  or a zero'd out #vrna_bp_t structure. The latter is returned on error
 *  or when @p n is outside the size of the stack.
 *
 *  @param  bps   The base pair stack
 *  @param  n     The position within the stack
 *  @return       The n'th element of the base pair stack, or a zero'd out #vrna_bp_t
 */
vrna_bp_t
vrna_bps_at(vrna_bps_t  bps,
            size_t      n);

/**
 *  @brief  Get the size of the base pair stack
 *
 *  @param  bps   The base pair stack
 *  @return       The size of the base pair stack
 */
size_t
vrna_bps_size(vrna_bps_t bps);

/**
 *  @}
 */

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
