#ifndef VIENNA_RNA_PACKAGE_PLOT_UTILS_H
#define VIENNA_RNA_PACKAGE_PLOT_UTILS_H

/**
 *  @file     plot_utils.h
 *  @ingroup  plotting_utils
 *  @brief    Various utilities to assist in plotting secondary structures and consensus structures
 */

/**
 *  @addtogroup  annotation_utils
 *  @{
 */


/**
 *  @brief  Produce covariance annotation for an alignment given a secondary structure
 *
 */
char **
vrna_annotate_bp_covar(const char **alignment,
                       const char *structure,
                       vrna_md_t  *md);


/**
 *  @brief  Produce covariance annotation for an alignment given a set of base pairs
 *
 */
vrna_cpair_t *
vrna_annotate_pr_covar(const char **alignment,
                       vrna_ep_t  *pl,
                       vrna_ep_t  *mfel,
                       double     threshold,
                       vrna_md_t  *md);


/**
 * @}
 */

#endif
