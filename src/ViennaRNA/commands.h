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

typedef enum {
  VRNA_CMD_ERROR=-1,
  VRNA_CMD_HC,
  VRNA_CMD_SC,
  VRNA_CMD_MOTIF,
  VRNA_CMD_UD,
  VRNA_CMD_SD
} vrna_command_e;

struct vrna_command_s {
  vrna_command_e  type;
  void *data;
};

void vrna_file_commands_read( const char *filename,
                              unsigned int options);

/**
 * @}
 */

#endif
