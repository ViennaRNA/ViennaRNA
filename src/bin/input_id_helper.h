#define ID_number_increase(num, name) ({ \
    if (num == LONG_MAX) { \
      vrna_message_warning(name " ID number overflow, beginning with 1 (again)!"); \
      num = 1; \
    } else { num++; } \
  })


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
