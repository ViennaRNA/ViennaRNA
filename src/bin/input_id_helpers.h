#ifndef VRNA_INPUT_ID_HELPERS
#define VRNA_INPUT_ID_HELPERS

typedef struct id_data *dataset_id;

void ID_number_increase(long int    *num,
                        const char  *name);


#define ID_generate(ID, \
                    ID_default, \
                    ID_auto_switch, \
                    ID_prefix, \
                    ID_delimiter, \
                    ID_digits, \
                    ID_number, \
                    ID_full)  ({ \
    if ((!ID_auto_switch) && (ID_default) && (ID_default[0] != '\0')) { \
      if (ID_full) {  /* use full default ID */ \
        ID = strdup(ID_default); \
      } else {        /* use first word from default ID */ \
        ID = (char *)vrna_alloc(sizeof(char) * (strlen(ID_default) + 1)); \
        (void)sscanf(ID_default, "%s", ID);  \
        ID = (char *)vrna_realloc(ID, sizeof(char) * (strlen(ID) + 1)); \
      } \
    } else if (ID_auto_switch) {  /* we have nuffin', Jon Snow (...so we simply generate an ID) */ \
      ID = vrna_strdup_printf("%s%s%0*ld", ID_prefix, ID_delimiter, ID_digits, ID_number); \
      free(ID_default);           /* reset default ID */ \
      ID_default = strdup(ID); \
    } \
  })


dataset_id
init_id_data(const char *name,
             const char *default_prefix,
             const char *default_delimiter,
             int        default_digits,
             long       default_start);


void
free_id_data(dataset_id dat);


void
set_auto_id(dataset_id  dat,
            int         sw);


void
set_id_prefix(dataset_id  dat,
              const char  *prefix);


void
set_id_delim(dataset_id dat,
             const char *delim);


void
set_id_digits(dataset_id  dat,
              int         digits);


void
set_id_start(dataset_id dat,
             long int   start);


const char *
get_id_delim(struct id_data *dat);


int
get_auto_id(struct id_data *dat);


void
set_next_id(char            **ID_default,
            struct id_data  *dat);


long int
get_current_id(struct id_data *dat);


char *
fileprefix_from_id(const char     *id,
                   struct id_data *dat,
                   int            ID_full);


char *
fileprefix_from_id_alifold(const char     *id,
                           struct id_data *dat,
                           int            continuous_names);


#endif
