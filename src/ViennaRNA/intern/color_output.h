#ifndef VRNA_INTERN_COLOR_OUTPUT_H
#define VRNA_INTERN_COLOR_OUTPUT_H

/* deactivate ANSI colors in TTY output if we compile for windows */
#ifndef VRNA_WITHOUT_TTY_COLORS
# ifdef _WIN32
#   define VRNA_WITHOUT_TTY_COLORS
# endif
#endif

#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

#ifndef VRNA_WITHOUT_TTY_COLORS

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

static INLINE void
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

static INLINE void
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


static INLINE void
print_eval_sd_corr(FILE *fp){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_BRIGHT "Correcting for presence of structured domains" ANSI_COLOR_RESET "\n");
  } else {
    fprintf(fp, "Correcting for presence of structured domains\n");
  }
}


static INLINE void
print_eval_ext_loop(FILE  *fp,
                    int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_CYAN "External loop" ANSI_COLOR_RESET
                "                           : "
                ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n", energy);
  } else {
    fprintf(fp, "External loop"
                "                           : "
                "%5d\n", energy);
  }
}

static INLINE void
print_eval_hp_loop( FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_CYAN "Hairpin  loop" ANSI_COLOR_RESET
                " (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                "              : "
                ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                i, j,
                si, sj,
                energy);
  } else {
    fprintf(fp, "Hairpin  loop"
                " (%3d,%3d) %c%c              : "
                "%5d\n",
                i, j,
                si, sj,
                energy);
  }
}


static INLINE void
print_eval_hp_loop_revert(FILE  *fp,
                          int   i,
                          int   j,
                          char  si,
                          char  sj,
                          int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_MAGENTA "Hairpin  loop" ANSI_COLOR_RESET
                " (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                "              : "
                ANSI_COLOR_RED "%5d" ANSI_COLOR_RESET "\n",
                i, j,
                si, sj,
                -energy);
  } else {
    fprintf(fp, "Hairpin  loop"
                " (%3d,%3d) %c%c              : "
                "%5d\n",
                i, j,
                si, sj,
                -energy);
  }
}


static INLINE void
print_eval_int_loop(FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   k,
                    int   l,
                    char  sk,
                    char  sl,
                    int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_CYAN "Interior loop" ANSI_COLOR_RESET
                " (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                "; (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                ": "
                ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                i, j,
                si, sj,
                k, l,
                sk, sl,
                energy);
  } else {
    fprintf(fp, "Interior loop"
                " (%3d,%3d) "
                "%c%c"
                "; (%3d,%3d) "
                "%c%c"
                ": "
                "%5d\n",
                i, j,
                si, sj,
                k, l,
                sk, sl,
                energy);
  }
}


static INLINE void
print_eval_int_loop_revert(FILE  *fp,
                          int   i,
                          int   j,
                          char  si,
                          char  sj,
                          int   k,
                          int   l,
                          char  sk,
                          char  sl,
                          int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_MAGENTA "Interior loop" ANSI_COLOR_RESET
                " (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                "; (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                ": "
                ANSI_COLOR_RED "%5d" ANSI_COLOR_RESET "\n",
                i, j,
                si, sj,
                k, l,
                sk, sl,
                -energy);
  } else {
    fprintf(fp, "Interior loop"
                " (%3d,%3d) "
                "%c%c"
                "; (%3d,%3d) "
                "%c%c"
                ": "
                "%5d\n",
                i, j,
                si, sj,
                k, l,
                sk, sl,
                -energy);
  }
}


static INLINE void
print_eval_mb_loop( FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_CYAN "Multi    loop" ANSI_COLOR_RESET
                " (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                "              : "
                ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                i, j,
                si, sj,
                energy);
  } else {
    fprintf(fp, "Multi    loop"
                " (%3d,%3d) %c%c              : "
                "%5d\n",
                i, j,
                si, sj,
                energy);
  }
}


static INLINE void
print_eval_mb_loop_revert(FILE  *fp,
                          int   i,
                          int   j,
                          char  si,
                          char  sj,
                          int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_MAGENTA "Multi    loop" ANSI_COLOR_RESET
                " (%3d,%3d) "
                ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                "              : "
                ANSI_COLOR_RED "%5d" ANSI_COLOR_RESET "\n",
                i, j,
                si, sj,
                -energy);
  } else {
    fprintf(fp, "Multi    loop"
                " (%3d,%3d) %c%c              : "
                "%5d\n",
                i, j,
                si, sj,
                -energy);
  }
}


static INLINE void
print_eval_gquad(FILE  *fp,
                 int   i,
                 int   L,
                 int   l[3],
                 int   energy){

  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_CYAN "G-Quadruplex " ANSI_COLOR_RESET
              " (%3d,%3d) "
              ANSI_COLOR_BRIGHT "L%d  " ANSI_COLOR_RESET
              "(%2d,%2d,%2d)  : "
              ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
              i, i + 4*L + l[0] + l[1] + l[2] - 1,
              L, l[0], l[1], l[2],
              energy);
  } else {
    fprintf(fp, "G-Quadruplex "
              " (%3d,%3d) "
              "L%d  "
              "(%2d,%2d,%2d)  : "
              "%5d\n",
              i, i + 4*L + l[0] + l[1] + l[2] - 1,
              L, l[0], l[1], l[2],
              energy);
  }
}


static INLINE void
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

static INLINE void
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

static INLINE void
print_fasta_header( FILE *fp,
                    const char *head){

  if(head){
    fprintf(fp, ">%s\n", head);
  }
}

static INLINE void
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


static INLINE void
print_eval_sd_corr(FILE *fp){

  fprintf(fp, "Correcting for presence of structured domains\n");
}


static INLINE void
print_eval_ext_loop(FILE  *fp,
                    int   energy){

  fprintf(fp, "External loop"
              "                           : "
              "%5d\n", energy);
}


static INLINE void
print_eval_hp_loop( FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   energy){

  fprintf(fp, "Hairpin  loop"
              " (%3d,%3d) %c%c              : "
              "%5d\n",
              i, j,
              si, sj,
              energy);
}


static INLINE void
print_eval_hp_loop_revert( FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   energy){

  print_eval_hp_loop(fp, i, j, si, sj, -energy);
}


static INLINE void
print_eval_int_loop(FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   k,
                    int   l,
                    char  sk,
                    char  sl,
                    int   energy){

  fprintf(fp, "Interior loop"
              " (%3d,%3d) "
              "%c%c"
              "; (%3d,%3d) "
              "%c%c"
              ": "
              "%5d\n",
              i, j,
              si, sj,
              k, l,
              sk, sl,
              energy);
}


static INLINE void
print_eval_int_loop_revert(FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   k,
                    int   l,
                    char  sk,
                    char  sl,
                    int   energy){

  print_eval_int_loop(fp, i, j, si, sj, k, l, sk, sl, -energy);
}


static INLINE void
print_eval_mb_loop( FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   energy){

  fprintf(fp, "Multi    loop"
              " (%3d,%3d) %c%c              : "
              "%5d\n",
              i, j,
              si, sj,
              energy);
}


static INLINE void
print_eval_mb_loop_revert( FILE  *fp,
                    int   i,
                    int   j,
                    char  si,
                    char  sj,
                    int   energy){

  print_eval_mb_loop(fp, i, j, si, sj, -energy);
}


static INLINE void
print_eval_gquad(FILE  *fp,
                 int   i,
                 int   L,
                 int   l[3],
                 int   energy){

  fprintf(fp, "G-Quadruplex "
              " (%3d,%3d) "
              "L%d  "
              "(%2d,%2d,%2d)  : "
              "%5d\n",
              i, i + 4*L + l[0] + l[1] + l[2] - 1,
              L, l[0], l[1], l[2],
              energy);
}


static INLINE void
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

static INLINE void
print_comment(FILE *fp,
              const char *line){

  if(line){
    fprintf(fp, "%s\n", line);
  }
}

#endif

#endif
