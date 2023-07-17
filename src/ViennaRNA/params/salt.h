#ifndef VIENNA_RNA_PACKAGE_LOOPS_SALT_H
#define VIENNA_RNA_PACKAGE_LOOPS_SALT_H

/**
 *  @file     params/salt.h
 *  @ingroup  energy_parameters_salt
 *  @brief    Functions to compute salt correction
 */

/**
 *  @addtogroup energy_parameters_salt
 *  @{
 */

#include <math.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/model.h"


/**
 *  @brief Get salt correction for a loop at a given salt concentration and temperature
 *
 *  @param L            backbone number in loop
 *  @param salt         salt concentration (M)
 *  @param T            absolute temperature (K)
 *  @param backbonelen  Backbone Length, phosphate-to-phosphate distance (typically 6 for RNA, 6.76 for DNA)
 *  
 *  @return      Salt correction for loop in dcal/mol
 */
double
vrna_salt_loop(int L, double salt, double T, double backbonelen);


/**
 *  @brief Get salt correction for a loop at a given salt concentration and temperature
 *
 *  This functions is same as vrna_salt_loop but returns rounded salt correction in integer
 *
 *  @see vrna_salt_loop
 *
 *  @param L            backbone number in loop
 *  @param salt         salt concentration (M)
 *  @param T            absolute temperature (K)
 *  @param backbonelen  Backbone Length, phosphate-to-phosphate distance (typically 6 for RNA, 6.76 for DNA)
 *  
 *  @return      Rounded salt correction for loop in dcal/mol
 */
int
vrna_salt_loop_int(int L, double salt, double T, double backbonelen);


/**
 *  @brief Get salt correction for a stack at a given salt concentration and temperature
 *
 *  @param salt   salt concentration (M)
 *  @param T      absolute temperature (K)
 *  @param hrise  Helical Rise (typically 2.8 for RNA, 3.4 for DNA)
 *  
 *  @return      Rounded salt correction for stack in dcal/mol
 */
int
vrna_salt_stack(double salt, double T, double hrise);


/**
 *  @brief Fit linear function to loop salt correction
 *
 *  For a given range of loop size (backbone number), we perform a linear fitting
 *  on loop salt correction
 *
 *  @f[ \text{Loop correction} \approx m \cdot L + b. @f]
 *
 *  @see vrna_salt_loop()
 *
 *  @param saltLoop  List of loop salt correction of size from 1
 *  @param lower     Define the size lower bound for fitting
 *  @param upper     Define the size upper bound for fitting
 *  @param m         pointer to store the parameter m in fitting result
 *  @param b         pointer to store the parameter b in fitting result
 */
void
vrna_salt_ml(double saltLoop[], int lower, int upper, int *m, int *b);


/**
 *  @brief Get salt correction for duplex initialization at a given salt concentration
 *
 *  @param md   Model details data structure that specfifies salt concentration in buffer (M)
 *  @return     Rounded correction for duplex initialization in dcal/mol
 */
int
vrna_salt_duplex_init(vrna_md_t *md);

/**
 *  @}
 */

#endif
