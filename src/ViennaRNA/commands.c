/*
    commands.c

    Various functions dealing with parsing and application of commands

    (c) 2016 Ronny Lorenz

    ViennaRNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/commands.h"


typedef void *(command_parser_function)(const char *line);

/* number of commands we currently know and are able to interpret */
#define NUM_COMMANDS  6

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

char known_commands[NUM_COMMANDS][3] = {
  "UD",
  "SD",
  "P",
  "F",
  "W",
  "E"
};

vrna_command_e known_command_types[NUM_COMMANDS] = {
  VRNA_CMD_UD,
  VRNA_CMD_SD,
  VRNA_CMD_HC,
  VRNA_CMD_HC,
  VRNA_CMD_HC,
  VRNA_CMD_SC
};

command_parser_function *parsers[NUM_COMMANDS] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE vrna_cmd_t
parse_command(const char *line){

  vrna_cmd_t cmd;
  int i;
  char command[3];

  command[0] = '\0';
  if(sscanf(line, "%2c", command) == 1){
    command[2] = '\0'; /* just a precaution */
    for(i = 0; i < 6; i++){
      if(!strncmp(known_commands[i], command, 2))
        break;
    }
  }

  printf("line: %s\n", line);
  if(i < 6){ /* command is known, so lets try to process it */
    cmd.type = known_command_types[i];
    cmd.data = (parsers[i]) ? parsers[i](line) : NULL;
  } else {
    cmd.type = VRNA_CMD_ERROR;
    cmd.data = NULL;
  }

  return cmd;
}

PUBLIC void
vrna_file_commands_read(const char *filename,
                        unsigned int options){

  FILE *fp;
  char *line;
  vrna_cmd_t cmd;

  line = NULL;

  if(!(fp = fopen(filename, "r"))){
    vrna_message_warning("Hard Constraints File could not be opened!");
    return;
  }

  /* let's go through the file line by line and parse the commands */
  while((line=vrna_read_line(fp))){
    switch(*line){
      /* skip comment lines */
      case '#': case '%': case ';': case '/': case '*': case ' ': case '\0':
        free(line);
        continue;
      default:
        cmd = parse_command((const char *)line);
        break;
    }

    if(cmd.type == VRNA_CMD_UD){ /* command is an instruction for unstructured domains feature */
      printf("got UD command: %s\n", line);
    }

    free(line);
  }

  free(line);
}
