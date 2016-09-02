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

#define VRNA_CMD_PARSE_HC      1U
#define VRNA_CMD_PARSE_SC      2U
#define VRNA_CMD_PARSE_UD      4U
#define VRNA_CMD_PARSE_SD      8U
#define VRNA_CMD_PARSE_DEFAULTS (   VRNA_CMD_PARSE_HC \
                                  | VRNA_CMD_PARSE_SC \
                                  | VRNA_CMD_PARSE_UD \
                                  | VRNA_CMD_PARSE_SD \
                                )


typedef enum {
  VRNA_CMD_ERROR=-1,
  VRNA_CMD_LAST=0,
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

vrna_cmd_t *vrna_file_commands_read(const char *filename,
                                    unsigned int options);

int vrna_commands_apply(vrna_fold_compound_t *vc,
                        vrna_cmd_t *commands);

void vrna_commands_free( vrna_cmd_t *commands);

/**
 * @}
 */

#endif
