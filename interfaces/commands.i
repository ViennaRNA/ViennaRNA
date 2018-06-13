/**********************************************/
/* BEGIN interface for commands               */
/**********************************************/

%ignore vrna_cmd_t;
%ignore vrna_command_s;

%rename (cmd)  vrna_command_s;

typedef struct {} vrna_command_s;

%nodefaultdtor vrna_command_s;

%extend vrna_command_s {
  vrna_command_s() {
    vrna_command_s *c = NULL;
    return c;
  }
  ~vrna_command_s(){
    vrna_commands_free($self);
  }
}


%rename (file_commands_read)  my_file_commands_read;
%rename (commands_free)       my_commands_free;

%{
  struct vrna_command_s *
  my_file_commands_read(std::string   filename,
                        unsigned int  options = VRNA_CMD_PARSE_DEFAULTS)
  {
    int i;

    return vrna_file_commands_read(filename.c_str(),
                                   options);
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_file_commands_read;
%feature("kwargs") my_file_commands_read;
#endif

struct vrna_command_s *
my_file_commands_read(std::string   filename,
                      unsigned int  options = VRNA_CMD_PARSE_DEFAULTS);

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") commands_apply;
%feature("kwargs") commands_apply;
#endif

  int
  commands_apply(struct vrna_command_s  *commands,
                 unsigned int           options = VRNA_CMD_PARSE_DEFAULTS)
  {
    return vrna_commands_apply($self,
                               commands,
                               options);
  }

#ifdef SWIGPYTHON
%feature("autodoc") file_commands_apply;
%feature("kwargs") file_commands_apply;
#endif

  int
  file_commands_apply(std::string   filename,
                      unsigned int  options = VRNA_CMD_PARSE_DEFAULTS)
  {
    return vrna_file_commands_apply($self,
                                    filename.c_str(),
                                    options);
  }
}


%constant unsigned int CMD_PARSE_DEFAULTS = VRNA_CMD_PARSE_DEFAULTS;
%constant unsigned int CMD_PARSE_HC       = VRNA_CMD_PARSE_HC;
%constant unsigned int CMD_PARSE_SC       = VRNA_CMD_PARSE_SC;
%constant unsigned int CMD_PARSE_SD       = VRNA_CMD_PARSE_SD;
%constant unsigned int CMD_PARSE_UD       = VRNA_CMD_PARSE_UD;

%include <ViennaRNA/commands.h>
