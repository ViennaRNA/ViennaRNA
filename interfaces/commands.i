/**********************************************/
/* BEGIN interface for commands               */
/**********************************************/

%inline %{
  typedef enum {
    CMD_ERROR = VRNA_CMD_ERROR,
    CMD_LAST  = VRNA_CMD_LAST,
    CMD_HC    = VRNA_CMD_HC,
    CMD_SC    = VRNA_CMD_SC,
    CMD_MOTIF = VRNA_CMD_MOTIF,
    CMD_UD    = VRNA_CMD_UD,
    CMD_SD    = VRNA_CMD_SD
  } command_e;
%}

%rename (cmd) vrna_cmd_t;

typedef struct { vrna_command_e type; } vrna_cmd_t;

namespace std {
  %template(CmdVector) std::vector<vrna_cmd_t>;
};

%nodefaultctor vrna_cmd_t;

%extend vrna_cmd_t {
  ~vrna_cmd_t(){
    /* convert vector back into array */
    vrna_cmd_t *cmds = (vrna_cmd_t *)vrna_alloc(sizeof(vrna_cmd_t) * 2);
    cmds[0].type = $self->type;
    cmds[0].data = $self->data;
    cmds[1].type = VRNA_CMD_LAST;
    vrna_commands_free(cmds);
  }

}

%rename (file_commands_read)  my_file_commands_read;
%rename (commands_free)       my_commands_free;

%{

  std::vector<vrna_cmd_t> my_file_commands_read(std::string filename, unsigned int options = VRNA_CMD_PARSE_DEFAULTS){
    int                     i;
    std::vector<vrna_cmd_t> cmd_list;
    vrna_cmd_t              *commands;

    commands = vrna_file_commands_read(filename.c_str(), options);

    for(i = 0; commands[i].type != VRNA_CMD_LAST; i++){
      cmd_list.push_back(commands[i]);
    }

    cmd_list.push_back(commands[i]);

    free(commands);

    return cmd_list;
  }

%}

std::vector<vrna_cmd_t> my_file_commands_read(std::string filename, unsigned int options = VRNA_CMD_PARSE_DEFAULTS);

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {

  int commands_apply(std::vector<vrna_cmd_t> commands, unsigned int options = VRNA_CMD_PARSE_DEFAULTS){

    return vrna_commands_apply($self, &(commands[0]), options);
  }

  int file_commands_apply(std::string filename, unsigned int options = VRNA_CMD_PARSE_DEFAULTS){

    return vrna_file_commands_apply($self, filename.c_str(), options);
  }

}


%constant unsigned int CMD_PARSE_DEFAULTS = VRNA_CMD_PARSE_DEFAULTS;
%constant unsigned int CMD_PARSE_HC       = VRNA_CMD_PARSE_HC;
%constant unsigned int CMD_PARSE_SC       = VRNA_CMD_PARSE_SC;
%constant unsigned int CMD_PARSE_SD       = VRNA_CMD_PARSE_SD;
%constant unsigned int CMD_PARSE_UD       = VRNA_CMD_PARSE_UD;

%include <ViennaRNA/commands.h>
