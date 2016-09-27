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

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE vrna_cmd_t parse_command(const char *line, int line_number, const char *filename);

PRIVATE void *parse_ud_command( const char *line);

PRIVATE void *parse_constraint_force(const char *line);
PRIVATE void *parse_constraint_prohibit(const char *line);
PRIVATE void *parse_constraint_con(const char *line);
PRIVATE void *parse_constraint_allow(const char *line);
PRIVATE void *parse_constraint_energy(const char *line);
PRIVATE void *parse_constraint(const char *line, char command);

PRIVATE int parse_constraints_line( const char *line,
                                    char command,
                                    int *i,
                                    int *j,
                                    int *k,
                                    int *l,
                                    char *loop,
                                    char *orientation,
                                    float *e);


PRIVATE int apply_hard_constraint(vrna_fold_compound_t *vc,
                                  void *constraint);

PRIVATE int apply_soft_constraint(vrna_fold_compound_t *vc,
                                  void *constraint);

PRIVATE int apply_ud(vrna_fold_compound_t *vc, void *data);

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

typedef struct {
  char                    cmd[3];
  vrna_command_e          type;
  command_parser_function *parser;
} parsable;


/* number of commands we currently know and are able to interpret */
#define NUM_COMMANDS  7

/* set of known parsable commands */
parsable known_commands[NUM_COMMANDS] = {
  /* cmd , type , parser */
  { "UD", VRNA_CMD_UD, parse_ud_command },          /* unstructured domain */
  { "SD", VRNA_CMD_SD, NULL },                      /* structured domain */
  { "P",  VRNA_CMD_HC, parse_constraint_prohibit }, /* prohibit base pairing */
  { "F",  VRNA_CMD_HC, parse_constraint_force },    /* force base pairing */
  { "C",  VRNA_CMD_HC, parse_constraint_con },      /* remove conflicting pairs/force nucleotide in loop context */
  { "A",  VRNA_CMD_HC, parse_constraint_allow },    /* allow (non-canonical) pairs */
  { "E",  VRNA_CMD_SC, parse_constraint_energy }    /* soft constraint */
};

typedef struct {
  int   i;
  int   j;
  int   k;
  int   l;
  int   size;
  char  loop;
  char  orientation;
  float e;
  char  command;
} constraint_struct;


typedef struct {
  char          *motif;
  float         motif_en;
  unsigned int  loop_type;
} ud_struct;

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC int
vrna_file_commands_apply( vrna_fold_compound_t *vc,
                          const char *filename,
                          unsigned int options){

  /** [Applying commands from file] */
  int         r;
  vrna_cmd_t  *cmds;

  cmds  = vrna_file_commands_read(filename, options);
  r     = vrna_commands_apply(vc, cmds, options);

  vrna_commands_free(cmds);

  return r;
  /** [Applying commands from file] */
}

PUBLIC vrna_cmd_t *
vrna_file_commands_read(const char *filename,
                        unsigned int options){
  FILE        *fp;
  char        *line;
  int         num_commands, max_commands, line_number, valid;
  vrna_cmd_t  cmd, *output;

  line_number   = 0;
  num_commands  = 0;
  max_commands  = 15;
  line          = NULL;

  if(!(fp = fopen(filename, "r"))){
    vrna_message_warning("Command File could not be opened!");
    return NULL;
  }

  output = (vrna_cmd_t *)vrna_alloc(sizeof(vrna_cmd_t) * max_commands);

  /* let's go through the file line by line and parse the commands */
  while((line=vrna_read_line(fp))){
    line_number++;
    switch(*line){
      /* skip comment lines */
      case '#': case '%': case ';': case '/': case '*': case ' ': case '\0':
        free(line);
        continue;
      default:
        cmd = parse_command((const char *)line, line_number, filename);
        break;
    }

    free(line);

    if(cmd.type == VRNA_CMD_LAST){ /* end of command list */
      break;
    } else { /* check whether command is valid in user-defined context */
      valid = 0;
      switch(cmd.type){
        case VRNA_CMD_HC: valid = options & VRNA_CMD_PARSE_HC;
                          break;

        case VRNA_CMD_SC: valid = options & VRNA_CMD_PARSE_SC;
                          break;

        case VRNA_CMD_UD: valid = options & VRNA_CMD_PARSE_UD;
                          break;

        case VRNA_CMD_SD: valid = options & VRNA_CMD_PARSE_SD;
                          break;

        default:          break;
      }
      
      if(valid){ /* add command to list */
        output[num_commands++] = cmd;

        /* increase length of command list if necessary */
        if(num_commands == max_commands){
          max_commands *= 1.2;
          output = (vrna_cmd_t *)vrna_realloc(output, sizeof(vrna_cmd_t) * max_commands);
        }
      }
    }
  }

  /* mark end of command list */
  output = (vrna_cmd_t *)vrna_realloc(output, sizeof(vrna_cmd_t) * (num_commands + 1));
  output[num_commands].type = VRNA_CMD_LAST;
  output[num_commands].data = NULL;

  /* cleanup */
  free(line);

  return output;
}


PUBLIC int
vrna_commands_apply(vrna_fold_compound_t *vc,
                    vrna_cmd_t *commands,
                    unsigned int options){

  int         r = 0;
  vrna_cmd_t  *ptr;

  if(vc && commands){
    for(ptr = commands; ptr->type != VRNA_CMD_LAST; ptr++){
      switch(ptr->type){
        case VRNA_CMD_HC:   if(options & VRNA_CMD_PARSE_HC)
                              r += apply_hard_constraint(vc, ptr->data);
                            break;

        case VRNA_CMD_SC:   if(options & VRNA_CMD_PARSE_SC)
                              r += apply_soft_constraint(vc, ptr->data);
                            break;

        case VRNA_CMD_UD:   if(options & VRNA_CMD_PARSE_UD)
                              r += apply_ud(vc, ptr->data);
                            break;

        default:            /* do nothing */
                            break;
      }
    }
  }
  
  return r;
}

PUBLIC void
vrna_commands_free( vrna_cmd_t *commands){

  vrna_cmd_t *ptr;

  if(commands){
    for(ptr = commands; ptr->type != VRNA_CMD_LAST; ptr++){
      switch(ptr->type){
        case VRNA_CMD_UD: {
                            ud_struct *d = (ud_struct *)ptr->data;
                            free(d->motif);
                            free(ptr->data);
                          }
                          break;

        default:          free(ptr->data);
                          break;
      }
    }
    free(commands);
  }
}


PRIVATE int
apply_ud(vrna_fold_compound_t *vc, void *data){


  ud_struct *d = (ud_struct *)data;
  vrna_ud_add_motif(vc, d->motif, d->motif_en, d->loop_type);
  
  return 1;
}


PRIVATE int
apply_hard_constraint(vrna_fold_compound_t *vc,
                      void *data){

  int               i, j, k, l, h, cnt1, cnt2, cnt3;
  int               num_hc_up, max_hc_up;
  vrna_hc_up_t      *hc_up;
  char              t, orientation;
  constraint_struct *constraint = (constraint_struct *)data;

  i           = constraint->i;
  j           = constraint->j;
  k           = constraint->k;
  l           = constraint->l;
  t           = constraint->loop;
  orientation = constraint->orientation;
  h           = constraint->size;

  /* actually apply constraints */
  if(h == 0){ /* range mode (prohibit pairs only) */
    for(cnt1 = i; cnt1 <= j; cnt1++)
      for(cnt2 = MAX2(cnt1 + 1, k); cnt2 <= l; cnt2++){
        vrna_hc_add_bp(vc, cnt1, cnt2, t);
      }
  } else {

    /* we'll collect hard constraints for unpairedness */
    num_hc_up = 0;
    max_hc_up = 15;
    hc_up     = vrna_alloc(sizeof(vrna_hc_up_t) * max_hc_up);

    for(cnt1 = i; cnt1 <= j; cnt1++)
      for(cnt2 = k; cnt2 <= l; cnt2++)
        for(cnt3 = h; cnt3 != 0; cnt3--){
          if(cnt2 == 0){  /* enforce unpairedness of nucleotide */
            /* just store this constraint, we'll apply it later */
            hc_up[num_hc_up].position = cnt1 + (cnt3 - 1);
            hc_up[num_hc_up].options  = t;
            num_hc_up++;
            if(num_hc_up == max_hc_up){ /* increase size of hc_up if necessary */
              max_hc_up  = 1.2 * max_hc_up;
              hc_up      = (vrna_hc_up_t *)vrna_realloc(hc_up, sizeof(vrna_hc_up_t) * max_hc_up);
            }
          } else if((i == j) && (j == k) && (k == l)){  /* enforce pairedness of nucleotide */
            int d = 0;
            if(orientation != '\0')
              d = (orientation == 'U') ? -1 : 1;
            vrna_hc_add_bp_nonspecific(vc, cnt1 + (cnt3 - 1), d, t);
          } else {  /* enforce / prohibit base pair */
            vrna_hc_add_bp(vc, cnt1 + (cnt3 - 1), cnt2 - (cnt3 - 1), t);
          }
        }

    /* add hard constraints for unpairedness */
    if(num_hc_up > 0){
      hc_up[num_hc_up].position = 0; /* mark end of list */
      vrna_hc_add_up_batch(vc, hc_up);
    }
    free(hc_up);

  }

  return 1;
}


PRIVATE int
apply_soft_constraint(vrna_fold_compound_t *vc,
                      void *data){

  int               i, j, k, l, h, cnt1, cnt2, cnt3;
  float             e;
  constraint_struct *constraint = (constraint_struct *)data;

  i           = constraint->i;
  j           = constraint->j;
  k           = constraint->k;
  l           = constraint->l;
  h           = constraint->size;
  e           = constraint->e;

  for(cnt1 = i; cnt1 <= j; cnt1++)
    for(cnt2 = k; cnt2 <= l; cnt2++)
      for(cnt3 = h; cnt3 != 0; cnt3--){
        if((cnt2 == 0) || ((i == j) && (j == k) && (k == l))){  /* enforce nucleotide constraint */
          vrna_sc_add_up(vc, cnt1 + (cnt3 - 1), e, VRNA_OPTION_DEFAULT);
        } else {  /* enforce base pair constraint */
          vrna_sc_add_bp(vc, cnt1 + (cnt3 - 1), cnt2 - (cnt3 - 1), e, VRNA_OPTION_DEFAULT);
        }
      }

  return 1;
}


PRIVATE vrna_cmd_t
parse_command(const char *line, int line_number, const char *filename){

  vrna_cmd_t cmd;
  int i, r;
  char command[3];

  command[0] = '\0';
  i = NUM_COMMANDS;

  r = sscanf(line, "%2c", command);
  if(r == 1){
    command[2] = '\0'; /* just a precaution */
    for(i = 0; i < NUM_COMMANDS; i++){
      if(!strncmp(known_commands[i].cmd, command, strlen(known_commands[i].cmd)))
        break;
    }
  }

  if(i < NUM_COMMANDS){ /* command is known, so lets try to process it */
    cmd.data = (known_commands[i].parser) ? known_commands[i].parser(line) : NULL;
    if(cmd.data)
      cmd.type = known_commands[i].type;
    else {
      vrna_message_warning("Ignoring invalid command in file \"%s\":\nline %d: %s", filename, line_number, line);
      cmd.type = VRNA_CMD_ERROR;
    }
  } else {
    vrna_message_warning("Ignoring unknown command in file \"%s\":\nline %d: %s", filename, line_number, line);
    cmd.type = VRNA_CMD_ERROR;
    cmd.data = NULL;
  }

  return cmd;
}


PRIVATE void *
parse_ud_command( const char *line){

  int           ret, entries_seen, max_entries, pos, pp;
  char          *ptr, *buffer;
  float         e;
  unsigned int  loop_type;
  ud_struct     *data;

  buffer        = (char *)vrna_alloc(sizeof(char) * (strlen(line) + 1));
  data          = (ud_struct *)vrna_alloc(sizeof(ud_struct));
  data->motif   = NULL;
  ret           = 0;  /* error indicator */
  entries_seen  = 0;  /* entries seen so far */
  max_entries   = 3;  /* expected number of entries */
  pos           = 2;  /* position relative to start of line */
  pp            = 0;

  while(!ret && (entries_seen < max_entries) && (sscanf(line+pos,"%s%n", buffer, &pp) == 1)){
    pos += pp;
    switch(entries_seen){
      case 0:       /* motif in IUPAC format */
                    data->motif = strdup(buffer);
                    break;

      case 1:       /* motif energy in kcal/mol */
                    if(sscanf(buffer, "%g", &e) == 1){
                      data->motif_en = e;
                    } else {
                      ret = 1;
                    }
                    break;

      case 2:       /* motif loop type */
                    loop_type = 0;
                    for(ptr=buffer; *ptr != '\0'; ptr++){
                      switch(*ptr){
                        case 'A'  : loop_type |= VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS;
                                    break;
                        case 'E'  : loop_type |= VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP;
                                    break;
                        case 'H'  : loop_type |= VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP;
                                    break;
                        case 'I'  : loop_type |= VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP;
                                    break;
                        case 'M'  : loop_type |= VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP;
                                    break;
                        default:    ret = 1;
                                    break;
                      }
                      if(ret)
                        break;
                    }
                    data->loop_type = loop_type;
                    break;
    }
    entries_seen++;
  }

  free(buffer);

  if(ret){
    free(data->motif);
    free(data);
    return NULL;
  }

  if(data->loop_type == 0){
    data->loop_type = VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS;
    vrna_message_warning("");
  }
  return (void *)data;
}

PRIVATE void *
parse_constraint_force(const char *line){

  return parse_constraint(line, 'F');
}

PRIVATE void *
parse_constraint_prohibit(const char *line){

  return parse_constraint(line, 'P');
}

PRIVATE void *
parse_constraint_con(const char *line){

  return parse_constraint(line, 'C');
}

PRIVATE void *
parse_constraint_allow(const char *line){

  return parse_constraint(line, 'A');
}

PRIVATE void *
parse_constraint_energy(const char *line){

  return parse_constraint(line, 'E');
}


PRIVATE void *
parse_constraint( const char *line,
                  char command){

  int               ret, i, j, k, l, h, valid;
  char              loop, orientation;
  float             e;
  constraint_struct *output;

  output = NULL;

  i = j = k = l = -1;
  orientation = '\0'; /* no orientation */
  e = 0.;

  ret = parse_constraints_line(line + 1, command, &i, &j, &k, &l, &loop, &orientation, &e);

  if(ret == 0){

    /* do something with the constraint we've just read */

    h = 1; /* helix length for pairs, or number of unpaired nucleotides */

    /* check indices */
    valid = 0;
    if(i > 0){
      if(j == -1){ /* i and range [k:l] */
        if((k > 0) && (l > 0)){
          if((k < l) && (i < k) && (orientation == '\0')){
            j     = i;
            valid = 1;
          }
        }
      } else if(k <= 0){ /* range [i:j] and l */
        if((i < j) && (j < l) && (orientation == '\0')){
          k     = l;
          valid = 1;
        }
      } else if(l <= 0){ /* helix of size k starting with pair (i,j), or segment [i:i+k-1] */
        if(i != j){
          if((j == 0) || (((j - i + 1) > 2*k) && (orientation == '\0'))){
            h     = k;
            k = l = j;
            j     = i;
            valid = 1;
          }
        }
      } else if((i < j) && (k < l) && (i <= k) && (j <= l) && (orientation == '\0')){  /* range [i:j] and [k:l] */
        if(command == 'P'){ /* we only allow this for 'prohibit pairing between two ranges' */
          h     = 0;
          valid = 1;
        }
      }
    }

    if(valid){
      /* nucleotide constraint? */
      if((k == 0) && (l == 0) && (i == j) && (h > 0)){
        /* set correct loop type context */
        switch(command){
          case 'P': break;
          case 'A': /* this case allows particular nucleotides to form non-canonical pairs */
                    loop |= VRNA_CONSTRAINT_CONTEXT_NO_REMOVE; /* do not remove possibility to stay unpaired */
                    /* fall through */
          case 'F': /* set i == j == k == l */
                    k = l = i;
                    break;
          case 'E': loop = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS; /* soft constraints are always applied for all loops */
                    break;
          case 'C': loop |= VRNA_CONSTRAINT_CONTEXT_ENFORCE;  /* enforce context dependency */
                    break;
          default:  break;
        }
      } else { /* base pair constraint */
        /* set correct loop type context */
        switch(command){
          case 'P': loop = ~loop; /* prohibit */
                    loop &= VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                    loop |= VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;  /* since we prohibit pairs, we do not want to remove incompatible pairs */
                    break;
          case 'F': loop |= VRNA_CONSTRAINT_CONTEXT_ENFORCE;  /* enforce */
                    break;
          case 'E': loop = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;  /* soft constraints are always applied for all loops */
                    break;
          case 'C': break;  /* remove conflicting pairs only */
          case 'A': loop |= VRNA_CONSTRAINT_CONTEXT_NO_REMOVE; /* since we allow pairs, we do not want to remove incompatible pairs */
                    break;
          default:  break;
        }
      }

      output = (constraint_struct *)vrna_alloc(sizeof(constraint_struct));
      output->command     = command;
      output->i           = i;
      output->j           = j;
      output->k           = k;
      output->l           = l;
      output->size        = h;
      output->loop        = loop;
      output->orientation = orientation;
      output->e           = e;
    }
  }

  return (void *)output;
}

PRIVATE int
parse_constraints_line( const char *line,
                        char command,
                        int *i,
                        int *j,
                        int *k,
                        int *l,
                        char *loop,
                        char *orientation,
                        float *e){

  int v1, v2;
  int ret = 0;
  int range_mode = 0;
  int pos = 0;
  int max_entries = 5;
  int entries_seen = 0;
  int pp;
  float energy;
  char buf[256], buf2[10], *c, tmp_loop;

  switch(command){
    case 'A':   /* fall through */
    case 'F':   /* fall through */
    case 'P':   max_entries = 5;
                break;
    case 'C':   /* fall through */
    case 'E':   max_entries = 4;
                break;
    default:    ret = 1;  /* error */
                break;
  }

  /* default to all loop types */
  *loop     = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  tmp_loop  = (char)0;

  /* now lets scan the entire line for content */
  while(!ret && (entries_seen < max_entries) && (sscanf(line+pos,"%15s%n", &buf[0], &pp) == 1)){
    pos += pp;
    switch(entries_seen){
      case 0: /* must be i, or range */
              if(sscanf(buf, "%d-%d%n", &v1, &v2, &pp) == 2){
                if(pp == strlen(buf)){
                  *i = v1;
                  *j = v2;
                  range_mode = 1;
                  --max_entries; /* no orientation allowed now */
                  break;
                }
              } else if(sscanf(buf, "%d%n", &v1, &pp) == 1){
                if(pp == strlen(buf)){
                  *i = v1;
                  break;
                }
              }
              ret = 1;
              break;
      case 1: /* must be j, or range */
              if(sscanf(buf, "%d-%d%n", &v1, &v2, &pp) == 2){
                if(pp == strlen(buf)){
                  *k = v1;
                  *l = v2;
                  if(!range_mode)
                    --max_entries; /* no orientation allowed now */
                  range_mode = 1;
                  break;
                }
              } else if(range_mode){
                if(sscanf(buf, "%d%n", &v1, &pp) == 1){
                  if(pp == strlen(buf)){
                    *l = v1;
                    break;
                  }
                }
              } else if(sscanf(buf, "%d%n", &v1, &pp) == 1){
                if(pp == strlen(buf)){
                  *j = v1;
                  break;
                }
              }
              ret = 1;
              break;
      case 2: /* skip if in range_mode */
              if(!range_mode){
                /* must be k */
                if(sscanf(buf, "%d%n", &v1, &pp) == 1){
                  if(pp == strlen(buf)){
                    *k = v1;
                    break;
                  }
                }
                ret = 1;
                break;
              } else {
                --max_entries;
                /* fall through */
              }
      case 3: 
              if(command == 'E'){ /* must be pseudo energy */
                if(sscanf(buf, "%g%n", &energy, &pp) == 1){
                  if(pp == strlen(buf)){
                    *e = energy;
                    break;
                  }
                }
              } else { /*  must be loop type, or orientation */
                if(sscanf(buf, "%8s%n", &buf2[0], &pp) == 1){
                  buf2[8] = '\0';
                  if(pp == strlen(buf)){
                    for(c = &(buf2[0]); (*c != '\0') && (!ret); c++){
                      switch(*c){
                        case 'E': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_EXT_LOOP;
                                  break;
                        case 'H': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_HP_LOOP;
                                  break;
                        case 'I': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_INT_LOOP;
                                  break;
                        case 'i': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
                                  break;
                        case 'M': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_MB_LOOP;
                                  break;
                        case 'm': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;
                                  break;
                        case 'A': tmp_loop |= VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                                  break;
                        case 'U': case 'D':
                                  *orientation = *c;
                                  break;
                        default:  ret = 1;
                      }
                    }
                    if(tmp_loop)
                      *loop = tmp_loop;

                    break;
                  }
                }
              }
              ret = 1;
              break;
      case 4: /* must be orientation */
              if(!(sscanf(buf, "%c", orientation) == 1)){
                ret = 1;
              }
              break;
    }
    ++entries_seen;
  }

  return ret;
}
