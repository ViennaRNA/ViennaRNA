#ifndef VIENNA_RNA_PACKAGE_MELTING_H
#define VIENNA_RNA_PACKAGE_MELTING_H

#include <stdio.h>

#include <ViennaRNA/datastructures/basic.h>

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
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
 *
 *  @file heat_capacity.h
 *  @ingroup  thermodynamics
 *
 *  @brief Compute heat capacity for an RNA
 *
 *  This file includes the interface to all functions related to predicting the heat capacity for an RNA.
 */


/**
 *  @addtogroup  thermodynamics
 *  @{
 */


/**
 *  @brief  The callback for heat capacity predictions
 *
 *  @callback
 *  @parblock
 *  This function will be called for each evaluated temperature in the heat capacity prediction.
 *  @endparblock
 *
 *  @see vrna_heat_capacity_cb()
 *
 *  @param temp           The current temperature this results corresponds to in &deg;C
 *  @param heat_capacity  The heat capacity in Kcal/(Mol * K)
 *  @param data           Some arbitrary data pointer passed through by the function executing the callback
 */
typedef void (*vrna_heat_capacity_f)(float  temp,
                                           float  heat_capacity,
                                           void   *data);

DEPRECATED(typedef void (vrna_heat_capacity_callback)(float  temp,
                                           float  heat_capacity,
                                           void   *data),
           "Use vrna_heat_capacity_f instead!");


/**
 *  @brief  A single result from heat capacity computations
 *
 *  This is a convenience typedef for #vrna_heat_capacity_s, i.e. results as obtained from vrna_heat_capacity()
 */
typedef struct vrna_heat_capacity_s vrna_heat_capacity_t;


/**
 *  @brief  A single result from heat capacity computations
 *
 *  @see vrna_heat_capacity()
 */
struct vrna_heat_capacity_s {
  float temperature;    /**< @brief   The temperature in &deg;C */
  float heat_capacity;  /**< @brief   The specific heat at this temperature in Kcal/(Mol * K) */
};


/**
 *  @name Basic heat capacity function interface
 *  @{
 */

/**
 *  @brief  Compute the specific heat for an RNA
 *
 *  This function computes an RNAs specific heat in a given temperature range
 *  from the partition function by numeric differentiation. The result is returned
 *  as a list of pairs of temperature in &deg;C and specific heat in Kcal/(Mol*K).
 *
 *  Users can specify the temperature range for the computation from @p T_min to
 *  @p T_max, as well as the increment step size @p T_increment. The latter also determines
 *  how many times the partition function is computed. Finally, the parameter @p mpoints
 *  determines how smooth the curve should be. The algorithm itself fits a parabola
 *  to @f$ 2 \cdot mpoints + 1 @f$ data points to calculate 2nd derivatives. Increasing this
 *  parameter produces a smoother curve.
 *
 *  @see  vrna_heat_capacity_cb(), vrna_heat_capacity_t, vrna_heat_capacity_s
 *
 *  @param  fc            The #vrna_fold_compound_t with the RNA sequence to analyze
 *  @param  T_min         Lowest temperature in &deg;C
 *  @param  T_max         Highest temperature in &deg;C
 *  @param  T_increment   Stepsize for temperature incrementation in &deg;C (a reasonable choice might be 1&deg;C)
 *  @param  mpoints       The number of interpolation points to calculate 2nd derivative (a reasonable choice might be 2, min: 1, max: 100)
 *  @return               A list of pairs of temperatures and corresponding heat capacity or @em NULL upon any failure.
 *                        The last entry of the list is indicated by a @b temperature field set to a value smaller than @p T_min
 */
vrna_heat_capacity_t *
vrna_heat_capacity(vrna_fold_compound_t *fc,
                   float                T_min,
                   float                T_max,
                   float                T_increment,
                   unsigned int         mpoints);


/**
 *  @brief  Compute the specific heat for an RNA  (callback variant)
 *
 *  Similar to vrna_heat_capacity(), this function computes an RNAs specific heat in
 *  a given temperature range from the partition function by numeric differentiation.
 *  Instead of returning a list of temperature/specific heat pairs, however, this
 *  function returns the individual results through a callback mechanism. The provided
 *  function will be called for each result and passed the corresponding temperature
 *  and specific heat values along with the arbitrary data as provided through the
 *  @p data pointer argument.
 *
 *  Users can specify the temperature range for the computation from @p T_min to
 *  @p T_max, as well as the increment step size @p T_increment. The latter also determines
 *  how many times the partition function is computed. Finally, the parameter @p mpoints
 *  determines how smooth the curve should be. The algorithm itself fits a parabola
 *  to @f$ 2 \cdot mpoints + 1 @f$ data points to calculate 2nd derivatives. Increasing this
 *  parameter produces a smoother curve.
 *
 *  @see  vrna_heat_capacity(), vrna_heat_capacity_f
 *
 *  @param  fc            The #vrna_fold_compound_t with the RNA sequence to analyze
 *  @param  T_min         Lowest temperature in &deg;C
 *  @param  T_max         Highest temperature in &deg;C
 *  @param  T_increment   Stepsize for temperature incrementation in &deg;C (a reasonable choice might be 1&deg;C)
 *  @param  mpoints       The number of interpolation points to calculate 2nd derivative (a reasonable choice might be 2, min: 1, max: 100)
 *  @param  cb            The user-defined callback function that receives the individual results
 *  @param  data          An arbitrary data structure that will be passed to the callback in conjunction with the results
 *  @return               Returns 0 upon failure, and non-zero otherwise
 */
int
vrna_heat_capacity_cb(vrna_fold_compound_t        *fc,
                      float                       T_min,
                      float                       T_max,
                      float                       T_increment,
                      unsigned int                mpoints,
                      vrna_heat_capacity_f cb,
                      void                        *data);


/* End basic interface */
/**@}*/

/**
 *  @name Simplified heat capacity computation
 *  @{
 */

/**
 *  @brief  Compute the specific heat for an RNA (simplified variant)
 *
 *  Similar to vrna_heat_capacity(), this function computes an RNAs specific heat
 *  in a given temperature range from the partition function by numeric differentiation.
 *  This simplified version, however, only requires the RNA sequence as input instead of
 *  a vrna_fold_compound_t data structure. The result is returned as a list of pairs of
 *  temperature in &deg;C and specific heat in Kcal/(Mol*K).
 *
 *  Users can specify the temperature range for the computation from @p T_min to
 *  @p T_max, as well as the increment step size @p T_increment. The latter also determines
 *  how many times the partition function is computed. Finally, the parameter @p mpoints
 *  determines how smooth the curve should be. The algorithm itself fits a parabola
 *  to @f$ 2 \cdot mpoints + 1 @f$ data points to calculate 2nd derivatives. Increasing this
 *  parameter produces a smoother curve.
 *
 *  @see  vrna_heat_capacity(), vrna_heat_capacity_cb(), vrna_heat_capacity_t, vrna_heat_capacity_s
 *
 *  @param  sequence      The RNA sequence input (must be uppercase)
 *  @param  T_min         Lowest temperature in &deg;C
 *  @param  T_max         Highest temperature in &deg;C
 *  @param  T_increment   Stepsize for temperature incrementation in &deg;C (a reasonable choice might be 1&deg;C)
 *  @param  mpoints       The number of interpolation points to calculate 2nd derivative (a reasonable choice might be 2, min: 1, max: 100)
 *  @return               A list of pairs of temperatures and corresponding heat capacity or @em NULL upon any failure.
 *                        The last entry of the list is indicated by a @b temperature field set to a value smaller than @p T_min
 */
vrna_heat_capacity_t *
vrna_heat_capacity_simple(const char    *sequence,
                          float         T_min,
                          float         T_max,
                          float         T_increment,
                          unsigned int  mpoints);

/* End basic interface */
/**@}*/

/* End thermodynamics */
/**@}*/

#endif
