#ifndef VIENNA_RNA_PACKAGE_CONVERT_EPARS_H
#define VIENNA_RNA_PACKAGE_CONVERT_EPARS_H

/**
 *  @file convert_epars.h
 *  @ingroup energy_parameters
 *  @brief Functions and definitions for energy parameter file format conversion
 */

/**
 *  @addtogroup energy_parameters_convert
 *  @brief      Convert energy parameter files into the latest format
 *
 *  To preserve some backward compatibility the RNAlib also provides
 *  functions to convert energy parameter files from the format used
 *  in version 1.4-1.8 into the new format used since version 2.0
 *
 *  @{
 *  @ingroup energy_parameters_convert
 */

/** Flag to indicate printing of a complete parameter set */
#define VRNA_CONVERT_OUTPUT_ALL           1U
/** Flag to indicate printing of hairpin contributions */
#define VRNA_CONVERT_OUTPUT_HP            2U
/** Flag to indicate printing of base pair stack contributions */
#define VRNA_CONVERT_OUTPUT_STACK         4U
/** Flag to indicate printing of  hairpin mismatch contribution  */
#define VRNA_CONVERT_OUTPUT_MM_HP         8U
/** Flag to indicate printing of  interior loop mismatch contribution  */
#define VRNA_CONVERT_OUTPUT_MM_INT        16U
/** Flag to indicate printing of  1:n interior loop mismatch contribution  */
#define VRNA_CONVERT_OUTPUT_MM_INT_1N     32U
/** Flag to indicate printing of  2:3 interior loop mismatch contribution  */
#define VRNA_CONVERT_OUTPUT_MM_INT_23     64U
/** Flag to indicate printing of  multi loop mismatch contribution  */
#define VRNA_CONVERT_OUTPUT_MM_MULTI      128U
/** Flag to indicate printing of  exterior loop mismatch contribution  */
#define VRNA_CONVERT_OUTPUT_MM_EXT        256U
/** Flag to indicate printing of  5' dangle conctribution  */
#define VRNA_CONVERT_OUTPUT_DANGLE5       512U
/** Flag to indicate printing of  3' dangle contribution  */
#define VRNA_CONVERT_OUTPUT_DANGLE3       1024U
/** Flag to indicate printing of  1:1 interior loop contribution  */
#define VRNA_CONVERT_OUTPUT_INT_11        2048U
/** Flag to indicate printing of  2:1 interior loop contribution */
#define VRNA_CONVERT_OUTPUT_INT_21        4096U
/** Flag to indicate printing of  2:2 interior loop contribution  */
#define VRNA_CONVERT_OUTPUT_INT_22        8192U
/** Flag to indicate printing of  bulge loop contribution  */
#define VRNA_CONVERT_OUTPUT_BULGE         16384U
/** Flag to indicate printing of  interior loop contribution  */
#define VRNA_CONVERT_OUTPUT_INT           32768U
/** Flag to indicate printing of  multi loop contribution  */
#define VRNA_CONVERT_OUTPUT_ML            65536U
/** Flag to indicate printing of  misc contributions (such as terminalAU)  */
#define VRNA_CONVERT_OUTPUT_MISC          131072U
/** Flag to indicate printing of  special hairpin contributions (tri-, tetra-, hexa-loops)  */
#define VRNA_CONVERT_OUTPUT_SPECIAL_HP    262144U
/** Flag to indicate printing of  given parameters only\n\note This option overrides all other output options, except #VRNA_CONVERT_OUTPUT_DUMP ! */
#define VRNA_CONVERT_OUTPUT_VANILLA       524288U
/** Flag to indicate printing of  interior loop asymmetry contribution  */
#define VRNA_CONVERT_OUTPUT_NINIO         1048576U
/** Flag to indicate dumping the energy contributions from the library instead of an input file  */
#define VRNA_CONVERT_OUTPUT_DUMP          2097152U

/**
 *  Convert/dump a Vienna 1.8.4 formatted energy parameter file
 * 
 *  The options argument allows one to control the different output modes.\n
 *  Currently available options are:\n
 *  #VRNA_CONVERT_OUTPUT_ALL, #VRNA_CONVERT_OUTPUT_HP, #VRNA_CONVERT_OUTPUT_STACK\n
 *  #VRNA_CONVERT_OUTPUT_MM_HP, #VRNA_CONVERT_OUTPUT_MM_INT, #VRNA_CONVERT_OUTPUT_MM_INT_1N\n
 *  #VRNA_CONVERT_OUTPUT_MM_INT_23, #VRNA_CONVERT_OUTPUT_MM_MULTI, #VRNA_CONVERT_OUTPUT_MM_EXT\n
 *  #VRNA_CONVERT_OUTPUT_DANGLE5, #VRNA_CONVERT_OUTPUT_DANGLE3, #VRNA_CONVERT_OUTPUT_INT_11\n
 *  #VRNA_CONVERT_OUTPUT_INT_21, #VRNA_CONVERT_OUTPUT_INT_22, #VRNA_CONVERT_OUTPUT_BULGE\n
 *  #VRNA_CONVERT_OUTPUT_INT, #VRNA_CONVERT_OUTPUT_ML, #VRNA_CONVERT_OUTPUT_MISC\n
 *  #VRNA_CONVERT_OUTPUT_SPECIAL_HP, #VRNA_CONVERT_OUTPUT_VANILLA, #VRNA_CONVERT_OUTPUT_NINIO\n
 *  #VRNA_CONVERT_OUTPUT_DUMP
 * 
 *  The defined options are fine for bitwise compare- and assignment-operations,
 *  e. g.: pass a collection of options as a single value like this:
 *  @verbatim convert_parameter_file(ifile, ofile, option_1 | option_2 | option_n) @endverbatim
 * 
 *  @param iname    The input file name (If NULL input is read from stdin)
 *  @param oname    The output file name (If NULL output is written to stdout)
 *  @param options  The options (as described above)
 */
void convert_parameter_file(const char *iname,
                            const char *oname,
                            unsigned int options);

/**
 *  @}
 */
#endif
