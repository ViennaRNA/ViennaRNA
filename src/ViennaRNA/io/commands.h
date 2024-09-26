#ifndef VIENNA_RNA_PACKAGE_IO_COMMANDS_H
#define VIENNA_RNA_PACKAGE_IO_COMMANDS_H

/**
 *  @file     ViennaRNA/io/commands.h
 *  @ingroup  utils, file_utils, command_files
 *  @brief    Parse and apply different commands that alter the behavior of
 *  secondary structure prediction and evaluation
 */

/**
 *  @addtogroup  command_files
 *  @{
 *  @brief  Functions to parse and interpret the content of @ref constraint-formats-file
 */

/** @brief A data structure that contains commands */
typedef struct vrna_command_s *vrna_cmd_t;


#include <ViennaRNA/fold_compound.h>

/**
 * @brief Command parse/apply flag indicating hard constraints
 *
 * @see   #vrna_cmd_t, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_HC      1U
/**
 * @brief Command parse/apply flag indicating soft constraints
 *
 * @see   #vrna_cmd_t, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_SC      2U
/**
 * @brief Command parse/apply flag indicating unstructured domains
 *
 * @see   #vrna_cmd_t, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_UD      4U
/**
 * @brief Command parse/apply flag indicating structured domains
 *
 * @see   #vrna_cmd_t, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_SD      8U
/**
 * @brief Command parse/apply flag indicating default set of commands
 *
 * @see   #vrna_cmd_t, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_DEFAULTS (VRNA_CMD_PARSE_HC \
                                 | VRNA_CMD_PARSE_SC \
                                 | VRNA_CMD_PARSE_UD \
                                 | VRNA_CMD_PARSE_SD \
                                 )

#define VRNA_CMD_PARSE_SILENT   16U

/**
 *  @brief Extract a list of commands from a command file
 *
 *  Read a list of commands specified in the input file
 *  and return them as list of abstract commands
 *
 *  @see  vrna_commands_apply(), vrna_file_commands_apply(),
 *        vrna_commands_free()
 *
 *  @param    filename  The filename
 *  @param    options   Options to limit the type of commands read from the file
 *  @return             A list of abstract commands
 */
vrna_cmd_t
vrna_file_commands_read(const char    *filename,
                        unsigned int  options);


/**
 *  @brief Apply a list of commands from a command file
 *
 *  This function is a shortcut to directly parse a commands file
 *  and apply all successfully parsed commands to a #vrna_fold_compound_t
 *  data structure. It is the same as:
 *  @snippet commands.c Applying commands from file
 *
 *  @param    fc        The #vrna_fold_compound_t the command list will be applied to
 *  @param    filename  The filename
 *  @param    options   Options to limit the type of commands read from the file
 *  @return             The number of commands successfully applied
 */
int
vrna_file_commands_apply(vrna_fold_compound_t *fc,
                         const char           *filename,
                         unsigned int         options);


/**
 *  @brief Apply a list of commands to a #vrna_fold_compound_t
 *
 *  @param    fc        The #vrna_fold_compound_t the command list will be applied to
 *  @param    commands  The commands to apply
 *  @param    options   Options to limit the type of commands read from the file
 *  @return             The number of commands successfully applied
 */
int
vrna_commands_apply(vrna_fold_compound_t  *fc,
                    vrna_cmd_t            commands,
                    unsigned int          options);


/**
 *  @brief Free memory occupied by a list of commands
 *
 *  Release memory occupied by a list of commands
 *  @param  commands  A pointer to a list of commands
 */
void
vrna_commands_free(vrna_cmd_t commands);


/**
 * @}
 */

#endif
