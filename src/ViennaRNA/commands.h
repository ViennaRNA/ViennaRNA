#ifndef VIENNA_RNA_PACKAGE_COMMANDS_H
#define VIENNA_RNA_PACKAGE_COMMANDS_H

/**
 *  @file     commands.h
 *  @ingroup  file_utils
 *  @brief    Parse and apply different commands that alter the behavior of
 *  secondary structure prediction and evaluation
 */

/**
 *  @{
 *  @ingroup  file_utils
 */

/** @brief Typename for the command repesenting data structure #vrna_command_s */
typedef struct vrna_command_s vrna_cmd_t;


#include <ViennaRNA/data_structures.h>

/**
 * @brief Command parse/apply flag indicating hard constraints
 * @see   #vrna_command_s, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_HC      1U
/**
 * @brief Command parse/apply flag indicating soft constraints
 * @see   #vrna_command_s, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_SC      2U
/**
 * @brief Command parse/apply flag indicating unstructured domains
 * @see   #vrna_command_s, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_UD      4U
/**
 * @brief Command parse/apply flag indicating structured domains
 * @see   #vrna_command_s, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_SD      8U
/**
 * @brief Command parse/apply flag indicating default set of commands
 * @see   #vrna_command_s, vrna_file_commands_read(), vrna_file_commands_apply(), vrna_commands_apply()
 */
#define VRNA_CMD_PARSE_DEFAULTS (   VRNA_CMD_PARSE_HC \
                                  | VRNA_CMD_PARSE_SC \
                                  | VRNA_CMD_PARSE_UD \
                                  | VRNA_CMD_PARSE_SD \
                                )

/**
 *  @brief Types of commands within a list of #vrna_command_s structures
 */
typedef enum {
  VRNA_CMD_ERROR=-1,
  VRNA_CMD_LAST=0,
  VRNA_CMD_HC,
  VRNA_CMD_SC,
  VRNA_CMD_MOTIF,
  VRNA_CMD_UD,
  VRNA_CMD_SD
} vrna_command_e;

/**
 *  @brief List element for commands ready for application to a #vrna_fold_compound_t
 *  @see vrna_file_commands_read(), vrna_commands_apply(), vrna_commands_free()
 */
struct vrna_command_s {
  vrna_command_e  type;
  void *data;
};

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
vrna_cmd_t *vrna_file_commands_read(const char *filename,
                                    unsigned int options);

/**
 *  @brief Apply a list of commands from a command file
 *
 *  This function is a shortcut to directly parse a commands file
 *  and apply all successfully parsed commands to a #vrna_fold_compound_t
 *  data structure. It is the same as:
 *  @snippet commands.c Applying commands from file
 *
 *  @param    vc        The #vrna_fold_compound_t the command list will be applied to
 *  @param    filename  The filename
 *  @param    options   Options to limit the type of commands read from the file
 *  @return             The number of commands successfully applied
 */
int vrna_file_commands_apply( vrna_fold_compound_t *vc,
                              const char *filename,
                              unsigned int options);

/**
 *  @brief Apply a list of commands to a #vrna_fold_compound_t
 *
 *  @param    vc        The #vrna_fold_compound_t the command list will be applied to
 *  @param    commands  The list of commands to apply
 *  @param    options   Options to limit the type of commands read from the file
 *  @return             The number of commands successfully applied
 */
int vrna_commands_apply(vrna_fold_compound_t *vc,
                        vrna_cmd_t *commands,
                        unsigned int options);

/**
 *  @brief Free memory occupied by a list of commands
 *
 *  Release memory occupied by a list of commands
 *  @param  commands  A pointer to a list of commands
 */
void vrna_commands_free( vrna_cmd_t *commands);

/**
 * @}
 */

#endif
