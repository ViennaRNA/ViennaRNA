
/* deactivate ANSI colors in TTY output if we compile for windows */
#ifndef WITHOUT_TTY_COLORS
# ifdef _WIN32
#   define WITHOUT_TTY_COLORS
# endif
#endif

#ifndef WITHOUT_TTY_COLORS

#define ANSI_COLOR_BRIGHT     "\x1b[1m"
#define ANSI_COLOR_UNDERLINE  "\x1b[4m"
#define ANSI_COLOR_RED        "\x1b[31m"
#define ANSI_COLOR_GREEN      "\x1b[32m"
#define ANSI_COLOR_YELLOW     "\x1b[33m"
#define ANSI_COLOR_BLUE       "\x1b[34m"
#define ANSI_COLOR_MAGENTA    "\x1b[35m"
#define ANSI_COLOR_CYAN       "\x1b[36m"
#define ANSI_COLOR_RED_B      "\x1b[1;31m"
#define ANSI_COLOR_GREEN_B    "\x1b[1;32m"
#define ANSI_COLOR_YELLOW_B   "\x1b[1;33m"
#define ANSI_COLOR_BLUE_B     "\x1b[1;34m"
#define ANSI_COLOR_MAGENTA_B  "\x1b[1;35m"
#define ANSI_COLOR_CYAN_B     "\x1b[1;36m"
#define ANSI_COLOR_RESET      "\x1b[0m"

static void
print_fasta_header( FILE *fp,
                    const char *head){

  if(head){
    if(isatty(fileno(fp))){
      fprintf(fp, ANSI_COLOR_YELLOW ">%s" ANSI_COLOR_RESET "\n", head);
    } else {
      fprintf(fp, ">%s\n", head);
    }
  }
}

static void
print_structure(FILE *fp,
                const char *structure,
                const char *data){

  if(structure){
    if(data){
      if(isatty(fileno(fp))){
        fprintf(fp, "%s" ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET "\n", structure, data);
      } else {
        fprintf(fp, "%s%s\n", structure, data);
      }
    } else {
      fprintf(fp, "%s\n", structure);
    }
  } else {
    if(data){
      if(isatty(fileno(fp))){
        fprintf(fp, ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET "\n", data);
      } else {
        fprintf(fp, "%s\n", data);
      }
    }
  }
}

static void
print_table(FILE *fp,
            const char *head,
            const char *line){

  if(head){
    if(isatty(fileno(fp))){
      fprintf(fp, ANSI_COLOR_UNDERLINE ANSI_COLOR_BRIGHT "%s" ANSI_COLOR_RESET "\n", head);
    } else {
      fprintf(fp, "%s\n", head);
    }
  }
  if(line){
    if(isatty(fileno(fp))){
      fprintf(fp, ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET "\n", line);
    } else {
      fprintf(fp, "%s\n", line);
    }
  }
}

static void
print_comment(FILE *fp,
              const char *line){

  if(line){
    if(isatty(fileno(fp))){
      fprintf(fp, ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", line);
    } else {
      fprintf(fp, "%s\n", line);
    }
  }
}

#else

static void
print_fasta_header( FILE *fp,
                    const char *head){

  if(head){
    fprintf(fp, ">%s\n", head);
  }
}

static void
print_structure(FILE *fp,
                const char *structure,
                const char *data){

  if(structure){
    if(data){
      fprintf(fp, "%s%s\n", structure, data);
    } else {
      fprintf(fp, "%s\n", structure);
    }
  } else {
    if(data){
      fprintf(fp, "%s\n", data);
    }
  }
}

static void
print_table(FILE *fp,
            const char *head,
            const char *line){

  if(head){
    fprintf(fp, "%s\n", head);
  }
  if(line){
    fprintf(fp, "%s\n", line);
  }
}

static void
print_comment(FILE *fp,
              const char *line){

  if(line){
    fprintf(fp, "%s\n", line);
  }
}

#endif

