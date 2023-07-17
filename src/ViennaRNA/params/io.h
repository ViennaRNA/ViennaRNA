#ifndef VIENNA_RNA_PACKAGE_PARAMS_IO_H
#define VIENNA_RNA_PACKAGE_PARAMS_IO_H

#ifdef VRNA_WARN_DEPRECATED
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
 *  @file     ViennaRNA/params/io.h
 *  @ingroup  energy_parameters
 *  @brief    Read and write energy parameter files
 */

/**
 *  @addtogroup energy_parameters_rw
 *  @{
 */


/**
 *  @brief  Default Energy Parameter File format
 *
 *  @see    vrna_params_load(), vrna_params_load_from_string(),
 *          vrna_params_save()
 */
#define VRNA_PARAMETER_FORMAT_DEFAULT     0


/**
 *  @brief Load energy parameters from a file
 *
 *  @see  vrna_params_load_from_string(), vrna_params_save(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @param  fname   The path to the file containing the energy parameters
 *  @param  options File format bit-mask (usually #VRNA_PARAMETER_FORMAT_DEFAULT)
 *  @return         Non-zero on success, 0 on failure
 */
int
vrna_params_load(const char   fname[],
                 unsigned int options);


/**
 *  @brief Save energy parameters to a file
 *
 *  @see vrna_params_load()
 *
 *  @param fname  A filename (path) for the file where the current energy parameters will be written to
 *  @param  options File format bit-mask (usually #VRNA_PARAMETER_FORMAT_DEFAULT)
 *  @return         Non-zero on success, 0 on failure
 */
int
vrna_params_save(const char   fname[],
                 unsigned int options);


/**
 *  @brief  Load energy paramters from string
 *
 *  The string must follow the default energy parameter file convention!
 *  The optional @p name argument allows one to specify a name for the
 *  parameter set which is stored internally.
 *
 *  @see  vrna_params_load(), vrna_params_save(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @param    string  A 0-terminated string containing energy parameters
 *  @param    name    A name for the parameter set in @p string (Maybe @p NULL)
 *  @param    options File format bit-mask (usually #VRNA_PARAMETER_FORMAT_DEFAULT)
 *  @return           Non-zero on success, 0 on failure
 */
int
vrna_params_load_from_string(const char   *string,
                             const char   *name,
                             unsigned int options);


/**
 *  @brief  Load default RNA energy parameter set
 *
 *  This is a convenience function to load the Turner 2004 RNA free
 *  energy parameters. It's the same as calling vrna_params_load_RNA_Turner2004()
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_defaults(void);


/**
 *  @brief  Load Turner 2004 RNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of RNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_defaults(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_RNA_Turner2004(void);


/**
 *  @brief  Load Turner 1999 RNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of RNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_RNA_Turner1999(void);


/**
 *  @brief  Load Andronsecu 2007 RNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of RNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_defaults(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_RNA_Andronescu2007(void);


/**
 *  @brief  Load Langdon 2018 RNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of RNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_RNA_Langdon2018(void);


/**
 *  @brief  Load Misc Special Hairpin RNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of RNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_defaults(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_RNA_misc_special_hairpins(void);


/**
 *  @brief  Load Mathews 2004 DNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of DNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_defaults(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_DNA_Mathews2004(void);


/**
 *  @brief  Load Mathews 1999 DNA energy parameter set
 *
 *  @warning  This function also resets the default geometric parameters
 *            as stored in #vrna_md_t to those of DNA. Only subsequently
 *            initialized #vrna_md_t structures will be affected by this
 *            change.
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_defaults()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int
vrna_params_load_DNA_Mathews1999(void);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief
 *
 */
enum parset {
  UNKNOWN= -1, QUIT,
  S, S_H, HP, HP_H, B, B_H, IL, IL_H, MMH, MMH_H, MMI, MMI_H,
  MMI1N, MMI1N_H, MMI23, MMI23_H, MMM, MMM_H, MME, MME_H, D5, D5_H, D3, D3_H,
  INT11, INT11_H, INT21, INT21_H, INT22, INT22_H, ML, TL,
  TRI, HEX, NIN, MISC
};


/**
 *  @brief Get the file name of the parameter file that was most recently loaded
 *
 *  @return The file name of the last parameter file, or NULL if parameters are still at defaults
 */
const char *
last_parameter_file(void);


/**
 *  @brief Read energy parameters from a file
 *
 *  @deprecated   Use vrna_params_load() instead!
 *  @param fname  The path to the file containing the energy parameters
 */
DEPRECATED(void
           read_parameter_file(const char fname[]),
           "Use vrna_params_load() instead!");


/**
 *  @brief Write energy parameters to a file
 *
 *  @deprecated   Use vrna_params_save() instead!
 *  @param fname  A filename (path) for the file where the current energy parameters will be written to
 */
DEPRECATED(void
           write_parameter_file(const char fname[]),
           "Use vrna_params_save() instead!");


/**
 *  @brief
 *
 */
enum  parset
gettype(const char *ident);


/**
 *  @brief
 *
 */
char *
settype(enum parset s);


/**
 *  @}
 */

#endif

#endif
