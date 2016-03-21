#ifndef VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H
#define VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H

#include "data_structures.h"

/**
 *  @file perturbation_fold.h
 *  @brief Find a vector of perturbation energies that minimizes the discripancies between predicted and observed pairing probabilities and the amount of neccessary adjustments
 *  @ingroup perturbation
 */

/**
 *  @addtogroup perturbation
 *  @brief Find a vector of perturbation energies that minimizes the discripancies between predicted and observed pairing probabilities and the amount of neccessary adjustments
 */

/**
 * @brief Use the sum of squared aberrations as objective function
 *
 * @f$ F(\vec\epsilon) = \sum_{i = 1}^n{ \frac{\epsilon_i^2}{\tau^2} } + \sum_{i = 1}^n{ \frac{(p_i(\vec\epsilon) - q_i)^2}{\sigma^2} } \to min @f$
 *
 * @ingroup perturbation
 */
#define VRNA_OBJECTIVE_FUNCTION_QUADRATIC 0

/**
 * @brief Use the sum of absolute aberrations as objective function
 *
 * @f$ F(\vec\epsilon) = \sum_{i = 1}^n{ \frac{|\epsilon_i|}{\tau^2} } + \sum_{i = 1}^n{ \frac{|p_i(\vec\epsilon) - q_i|}{\sigma^2} } \to min @f$
 *
 * @ingroup perturbation
 */
#define VRNA_OBJECTIVE_FUNCTION_ABSOLUTE 1

/**
 * @brief Use a custom implementation of the gradient descent algorithm to minimize the objective function
 *
 * @ingroup perturbation
 */
#define VRNA_MINIMIZER_DEFAULT 0

/**
 * @brief Use the GNU Scientific Library implementation of the Fletcher-Reeves conjugate gradient algorithm to minimize the objective function
 *
 * Please note that this algorithm can only be used when the GNU Scientific Library is available on your system
 *
 * @ingroup perturbation
 */
#define VRNA_MINIMIZER_CONJUGATE_FR 1

/**
 * @brief Use the GNU Scientific Library implementation of the Polak-Ribiere conjugate gradient algorithm to minimize the objective function
 *
 * Please note that this algorithm can only be used when the GNU Scientific Library is available on your system
 *
 * @ingroup perturbation
 */
#define VRNA_MINIMIZER_CONJUGATE_PR 2

/**
 * @brief Use the GNU Scientific Library implementation of the vector Broyden-Fletcher-Goldfarb-Shanno algorithm to minimize the objective function
 *
 * Please note that this algorithm can only be used when the GNU Scientific Library is available on your system
 *
 * @ingroup perturbation
 */
#define VRNA_MINIMIZER_VECTOR_BFGS 3

/**
 * @brief Use the GNU Scientific Library implementation of the vector Broyden-Fletcher-Goldfarb-Shanno algorithm to minimize the objective function
 *
 * Please note that this algorithm can only be used when the GNU Scientific Library is available on your system
 *
 * @ingroup perturbation
 */
#define VRNA_MINIMIZER_VECTOR_BFGS2 4

/**
 * @brief Use the GNU Scientific Library implementation of the steepest descent algorithm to minimize the objective function
 *
 * Please note that this algorithm can only be used when the GNU Scientific Library is available on your system
 *
 * @ingroup perturbation
 */
#define VRNA_MINIMIZER_STEEPEST_DESCENT 5

/**
 * @brief Callback for following the progress of the minimization process
 *
 * @param iteration The number of the current iteration
 * @param score     The score of the objective function
 * @param epsilon   The perturbation vector yielding the reported score
 *
 * @ingroup perturbation
 */
typedef void (*progress_callback)(int iteration, double score, double *epsilon);

/**
 *  @brief Find a vector of perturbation energies that minimizes the discripancies between predicted and observed pairing probabilities and the amount of neccessary adjustments
 *
 *  Use an iterative minimization algorithm to find a vector of perturbation energies whose incorporation as soft constraints shifts the predicted
 *  pairing probabilities closer to the experimentally observed probabilities.
 *  The algorithm aims to minimize an objective function that penalizes discripancies between predicted and observed pairing probabilities and energy model adjustments,
 *  i.e. an appropriate vector of perturbation energies satisfies
 *  @f[
 *  F(\vec\epsilon) = \sum_{\mu}{ \frac{\epsilon_{\mu}^2}{\tau^2} } + \sum_{i =
 *  1}^n{ \frac{(p_i(\vec\epsilon) - q_i)^2}{\sigma^2} } \to \min.
 *  @f]
 *
 *  An initialized fold compound and an array containing the observed probability for each nucleotide to be unbound are required as input data.
 *  The parameters objective_function, sigma_squared and tau_squared are responsible for adjusting the aim of the objective function.
 *  Dependend on which type of objective function is selected, either squared or absolute aberrations are contributing to the objective function.
 *  The ratio of the parameters sigma_squared and tau_squared can be used to adjust the algorithm to find a solution either close to the thermodynamic prediction
 *  (sigma_squared >> tau_squared) or close to the experimental data (tau_squared >> sigma_squared).
 *  The minimization can be performed by makeing use of a custom gradient descent implementation or using one of the minimizing algorithms provided by the GNU Scientific Library.
 *  All algorithms require the evaluation of the gradient of the objective function, which includes the evaluation of conditional pairing probabilites.
 *  Since an exact evaluation is expensive, the probabilities can also be estimated from sampling by setting an appropriate sample size.
 *  The found vector of perturbation energies will be stored in the array epsilon.
 *  The progress of the minimization process can be tracked by implementing and passing a callback function.
 *
 *  @see For further details we refere to @cite washietl:2012.
 *  @ingroup perturbation
 *
 *  @param vc                 Pointer to a fold compound
 *  @param q_prob_unpaired    Pointer to an array containing the probability to be unpaired for each nucleotide
 *  @param objective_function The type of objective function to be used (VRNA_OBJECTIVE_FUNCTION_QUADRATIC / VRNA_OBJECTIVE_FUNCTION_LINEAR)
 *  @param sigma_squared      A factor used for weighting the objective function.
 *                            More weight on this factor will lead to a solution close to the null vector.
 *  @param tau_squared        A factor used for weighting the objective function.
 *                            More weight on this factor will lead to a solution close to the data provided in q_prob_unpaired.
 *  @param algorithm          The minimization algorithm (VRNA_MINIMIZER_*)
 *  @param sample_size        The number of sampled sequences used for estimating the pairing probabilities. A value <= 0 will lead to an exact evaluation.
 *  @param epsilon            A pointer to an array used for storing the calculated vector of perturbation energies
 *  @param callback           A pointer to a callback function used for reporting the current minimization progress
 *
 */
void vrna_sc_minimize_pertubation(vrna_fold_compound_t *vc,
                                  const double *q_prob_unpaired,
                                  int objective_function,
                                  double sigma_squared,
                                  double tau_squared,
                                  int algorithm,
                                  int sample_size,
                                  double *epsilon,
                                  double initialStepSize,
                                  double minStepSize,
                                  double minImprovement,
                                  double minimizerTolerance,
                                  progress_callback callback);

#endif
