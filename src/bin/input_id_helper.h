#define ID_number_increase(num, name) ({ \
    if (num == LONG_MAX) { \
      vrna_message_warning(name " ID number overflow, beginning with 1 (again)!"); \
      num = 1; \
    } else { num++; } \
  })


#define ID_generate(ID, \
                    fname, \
                    ID_auto_switch, \
                    ID_prefix, \
                    ID_delimiter, \
                    ID_digits, \
                    ID_number)  ({ \
    if ((!ID_auto_switch) && (fname[0] != '\0')) { \
      ID = strdup(fname);         /* we've got an ID from somewhere, so we use it */ \
    } else if (ID_auto_switch) {  /* we have nuffin', Jon Snow (...so we simply generate an ID) */ \
      ID = vrna_strdup_printf("%s%s%0*ld", ID_prefix, ID_delimiter, ID_digits, ID_number); \
    } \
  })
