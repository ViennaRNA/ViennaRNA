/*
 *           Compute duplex structure of two RNA strands
 *
 *                         c Ivo L Hofacker
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/ali_plex.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/plex.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/read_epars.h"
#include "RNAplex_cmdl.h"


clock_t
BeginTimer()
{
  /* timer declaration */
  clock_t Begin;    /* initialize Begin */

  Begin = clock();  /* start the timer */
  return Begin;
}


clock_t
EndTimer(clock_t begin)
{
  clock_t End;

  End = clock();    /* stop the timer */
  return End;
}


/* --------------------end include timer */
extern int subopt_sorted;
/* static int print_struc(duplexT const *dup); */
static int **average_accessibility_target(char      **names,
                                          char      **ALN,
                                          int       number,
                                          char      *access,
                                          double    verhaeltnis,
                                          const int alignment_length,
                                          int       binaries,
                                          int       fast);


/* static int ** average_accessibility_query(char **names, char **ALN, int number, char *access, double verhaeltnis); */
static int get_max_u(const char *s,
                     char       delim);


static int **read_plfold_i(char       *fname,
                           const int  beg,
                           const int  end,
                           double     verhaeltnis,
                           const int  length,
                           int        fast);


/**
 * Compute the conditional per nucleotide directly
 */
/* static int ** read_plfold_j(char *fname, const int beg, const int end,double verhaeltnis); */
static int **read_plfold_i_bin(char       *fname,
                               const int  beg,
                               const int  end,
                               double     verhaeltnis,
                               const int  length,
                               int        fast);


/* Compute and pass opening energies in case of f=2*/
static int get_sequence_length_from_alignment(char *sequence);


/* take as argument a list of hit from an alignment interaction */
/* Accessibility can currently not be provided */
static void aliprint_struct(FILE      *Result,  /* result file */
                            FILE      *Target,  /* target alignment */
                            FILE      *Query,
                            const int WindowsLength);


static void linear_fit(int  *a,
                       int  *b,
                       int  *c,
                       int  *d);                        /* get linear fits for Iopn, Iextension, Bopn, Bextension, I^1open and I^1extension */


/**
 * Compute Tm based on silvana's parameters (1999)
 */
double probcompute_silvana_98(char    *s1,
                              double  k_concentration,
                              double  tris_concentration,
                              double  mg_concentration,
                              double  na_concentration,
                              double  probe_concentration);


double probcompute_silvana_04(char    *s1,
                              double  k_concentration,
                              double  tris_concentration,
                              double  mg_concentration,
                              double  na_concentration,
                              double  probe_concentration);


double     probcompute_xia_98(char    *s1,
                              double  na_concentration,
                              double  probe_concentration);


double     probcompute_sug_95(char    *s1,
                              double  na_concentration,
                              double  probe_concentration);


double probcompute_newparameters(char   *s1,
                                 double k_concentration,
                                 double tris_concentration,
                                 double mg_concentration,
                                 double na_concentration,
                                 double probe_concentration);


/*@unused@*/
static int convert_plfold_i(char *fname);/* convert test accessibility into bin accessibility. */


static char scale[] = "....,....1....,....2....,....3....,....4"
                      "....,....5....,....6....,....7....,....8";

/*--------------------------------------------------------------------------*/

int
main(int  argc,
     char *argv[])
{
  struct        RNAplex_args_info args_info;

#define MAX_NUM_NAMES    500
  char                            *temp1[MAX_NUM_NAMES], *temp2[MAX_NUM_NAMES], *AS1[MAX_NUM_NAMES], *AS2[MAX_NUM_NAMES], *names1[MAX_NUM_NAMES], *names2[MAX_NUM_NAMES];
  char                            *s1     = NULL, *s2 = NULL, *line, *cstruc = NULL, *structure = NULL;
  char                            *line_q = NULL, *line_t = NULL;
  char                            *tname  = NULL;
  char                            *qname  = NULL;
  char                            *access = NULL;
  char                            fname[FILENAME_MAX_LENGTH];
  char                            *ParamFile  = NULL;
  char                            *ns_bases   = NULL, *c;

  FILE                            *Result = NULL; /* file containing the results */
  FILE                            *mRNA   = NULL, *sRNA = NULL;

  char                            *Resultfile       = NULL;
  int                             fold_constrained  = 0;
  int                             i, l, sym;

  /* double kT, sfact=1.07; */
  int                             istty, delta = -INF;
  int                             noconv            = 0;
  int                             extension_cost    = 0;
  int                             deltaz            = 0;
  int                             alignment_length  = 40;
  int                             fast              = 0;
  int                             redraw            = 0;
  int                             binaries          = 0;
  int                             convert           = 0;
  /**
   * Defines how many nucleotides has to be added at the begining and end of the target and query sequence in order to generate the structure figure
   */
  int                             WindowsLength   = 0;
  int                             alignment_mode  = 0;
  /**
   * Define scaling factor for the accessibility. Default 1
   */
  double                          verhaeltnis = 1;
  /**
   * Probe Tm computation
   */
  double                          probe_concentration = 0.05;
  double                          na_concentration    = 1;
  double                          mg_concentration    = 0;
  double                          k_concentration     = 0;
  double                          tris_concentration  = 0;
  int                             probe_mode          = 0;
  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAplex_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /*temperature*/
  if (args_info.temp_given)
    temperature = args_info.temp_arg;

  /*query file*/
  if (args_info.query_given)
    qname = strdup(args_info.query_arg);

  /*target file*/
  if (args_info.target_given)
    tname = strdup(args_info.target_arg);

  /*interaction_length*/
  alignment_length = args_info.interaction_length_arg;
  /*extension_cost*/
  extension_cost = args_info.extension_cost_arg;
  /*duplex_distance*/
  deltaz = args_info.duplex_distance_arg;
  /*energy_threshold*/
  delta = (int)(100 * args_info.energy_threshold_arg);
  /*fast_folding*/
  fast = args_info.fast_folding_arg;
  /*accessibility*/
  if (args_info.accessibility_dir_given)
    access = strdup(args_info.accessibility_dir_arg);

  /*produce ps arg*/
  if (args_info.produce_ps_given) {
    Resultfile  = strdup(args_info.produce_ps_arg);
    redraw      = 1;
  }

  /*WindowLength*/
  WindowsLength = args_info.WindowLength_arg;
  /*scale_accessibility_arg*/
  verhaeltnis = args_info.scale_accessibility_arg;
  /*constraint_flag*/
  if (args_info.constraint_given)
    fold_constrained = 1;

  /*paramFile*/
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  /*binary*/
  if (args_info.binary_given)
    binaries = 1;

  /*convert_to_bin*/
  if (args_info.convert_to_bin_given)
    convert = 1;

  /*alignment_mode*/
  if (args_info.alignment_mode_given)
    alignment_mode = 1;

  /*Probe mode*/
  if (args_info.probe_mode_given)
    probe_mode = 1;

  /*sodium concentration*/
  na_concentration = args_info.na_concentration_arg;
  /*magnesium concentration*/
  mg_concentration = args_info.mg_concentration_arg;
  /*potassium concentration*/
  k_concentration = args_info.k_concentration_arg;
  /*tris concentration*/
  tris_concentration = args_info.tris_concentration_arg;
  /*probe concentration*/
  probe_concentration = args_info.probe_concentration_arg;

  /*Probe mode Salt concentration*/
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards  = vrna_alloc(33);
    c             = ns_bases;
    i             = sym = 0;
    if (*c == '-') {
      sym = 1;
      c++;
    }

    while (*c != '\0') {
      if (*c != ',') {
        nonstandards[i++] = *c++;
        nonstandards[i++] = *c;
        if ((sym) && (*c != *(c - 1))) {
          nonstandards[i++] = *c;
          nonstandards[i++] = *(c - 1);
        }
      }

      c++;
    }
  }

  int il_a, il_b, b_a, b_b;
  linear_fit(&il_a, &il_b, &b_a, &b_b);
  /**
   * check if we have two input files
   */


  /**
   * Here we test if the user wants to convert a bunch of text opening energy files
   * into binary
   */
  if (probe_mode) {
    if (qname || tname) {
      vrna_message_error("No query/target file allowed in Tm probe mode\nPlease pipe your input into RNAplex\n");
      /* get sequence */
    } else {
      /* fix temperature to 37C */
      temperature = 37;
      printf("Probe mode\n");
      char          *id_s1  = NULL;
      char          *s1     = NULL;
      vrna_param_t  *P      = NULL;
      vrna_md_t     md;
      set_model_details(&md);

      if ((!P) || (fabs(P->temperature - temperature) > 1e-6)) {
        update_fold_params();
        P = vrna_params(&md);
        make_pair_matrix();
      }

      /*Initialize parameter */
      printf("Concentration K:%3.3f TNP:%3.3f Mg:%3.3f Na:%3.3f probe:%3.3f\n\n", k_concentration, tris_concentration, mg_concentration, na_concentration, probe_concentration);
      printf("%100s %7s %7s %7s %7s %7s\n", "sequence", "DDSL98", "DDSL04", "DRSU95", "RRXI98", "CURRENT");
      do {
        istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
        if ((line = vrna_read_line(stdin)) == NULL)
          break;

        /* skip empty lines, comment lines, name lines */
        while ((*line == '*') || (*line == '\0') || (*line == '>')) {
          printf("%s\n", line);
          if (*line == '>') {
            id_s1 = (char *)vrna_alloc(strlen(line) + 2);
            (void)sscanf(line, "%s", id_s1);
            memmove(id_s1, id_s1 + 1, strlen(id_s1));
          }

          free(line);
          if ((line = vrna_read_line(stdin)) == NULL) {
            free(id_s1);
            break;
          }
        }
        if ((line == NULL) || (strcmp(line, "@") == 0))
          break;

        s1 = (char *)vrna_alloc(strlen(line) + 1);
        strcpy(s1, line);
        /*compute duplex/entropy energy for the reverse complement*/;
        double Tm;
        Tm = probcompute_silvana_98(line, k_concentration, tris_concentration, mg_concentration, na_concentration, probe_concentration);
        printf("%100s  %*.2f ", s1, 6, Tm);
        Tm = probcompute_silvana_04(line, k_concentration, tris_concentration, mg_concentration, na_concentration, probe_concentration);
        printf("%*.2f  ", 6, Tm);
        Tm = probcompute_sug_95(line, na_concentration, probe_concentration);
        printf("%*.2f  ", 6, Tm);
        Tm = probcompute_xia_98(line, na_concentration, probe_concentration);
        printf("%*.2f  ", 6, Tm);
        Tm = probcompute_newparameters(line, k_concentration, tris_concentration, mg_concentration, na_concentration, probe_concentration);
        printf("%*.2f\n", 6, Tm);
      } while (1);
    }

    RNAplex_cmdline_parser_free(&args_info);
    return 0;
  }

  if (convert && access) {
    char          pattern[8];
    strcpy(pattern, "_openen");
    DIR           *dfd;
    struct dirent *dir;
    dfd = opendir(access);
    /**
     * Check if the directory where the opening energy files are located exists
     */
    if (dfd) {
      while ((dir = readdir(dfd)) != NULL) {
        int   strPos      = 0;
        int   PatternPos  = 0;
        char  name[128];
        strcpy(name, access);
        strcat(name, "/");
        strcat(name, dir->d_name);
        while (name[strPos] != '\0') {
          if (name[strPos] == pattern[0]) {
            /**
             * Check if the files we are looking ends in openen and not openen_bin,
             * in order to avoid that RNAplex converts bin-file to bin-file
             */
            while (name[strPos + PatternPos] == pattern[PatternPos] && pattern[PatternPos] != '\0' && name[strPos + PatternPos] != '\0')
              PatternPos++;
            if (name[strPos + PatternPos] == '\0' && pattern[PatternPos] == '\0')
              /**
               * convert_plfold_i is the function that makes the hardwork
               */
              convert_plfold_i(name);
            else
              PatternPos = 0;
          }

          strPos++;
        }
      }
    }

    closedir(dfd);
    RNAplex_cmdline_parser_free(&args_info);
    return 0;
  }

  /**
   * Here we check if the user wants to produce PS formatted structure files
   * from existing RNAplex dot-parenthesis-formated results. Depending on
   * the kind of input, Alignments or single sequence we will produce
   * either a color annotated alignment or RNAfold-like structure,
   * respectively.
   */
  if (Resultfile) {
    /**
     * Check single sequence case.
     */
    if (!(alignment_mode) && (tname && qname)) {
      mRNA = fopen(tname, "r");
      if (mRNA == NULL) {
        printf("%s: Wrong target file name\n", tname);
        RNAplex_cmdline_parser_free(&args_info);
        return 0;
      }

      sRNA = fopen(qname, "r");
      if (sRNA == NULL) {
        printf("%s: Wrong query file name\n", qname);
        RNAplex_cmdline_parser_free(&args_info);
        return 0;
      }

      vrna_message_error("Sorry not implemented yet");
    }/**
      * We have no single sequence case. Check if we have alignments.
      */
    else if ((alignment_mode) && (tname && qname)) {
      mRNA = fopen(tname, "r");
      if (mRNA == NULL) {
        printf("%s: Wrong target file name\n", tname);
        RNAplex_cmdline_parser_free(&args_info);
        return 0;
      }

      sRNA = fopen(qname, "r");
      if (sRNA == NULL) {
        printf("%s: Wrong query file name\n", qname);
        RNAplex_cmdline_parser_free(&args_info);
        return 0;
      }

      Result = fopen(Resultfile, "r");
      if (sRNA == NULL) {
        printf("%s: Wrong query file name\n", qname);
        RNAplex_cmdline_parser_free(&args_info);
        return 0;
      }

      aliprint_struct(Result, mRNA, sRNA, (const int)WindowsLength);
      RNAplex_cmdline_parser_free(&args_info);
      return 0;
    } else {
      /**
       * User was not able to input either two alignments or two single sequence files
       */
      printf("Please enter either two alignments or single sequence files\n");
      RNAplex_cmdline_parser_print_help();
    }
  }

#if 0
  update_fold_params();
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards  = vrna_alloc(33);
    c             = ns_bases;
    i             = sym = 0;
    if (*c == '-') {
      sym = 1;
      c++;
    }

    while (*c != '\0') {
      if (*c != ',') {
        nonstandards[i++] = *c++;
        nonstandards[i++] = *c;
        if ((sym) && (*c != *(c - 1))) {
          nonstandards[i++] = *c;
          nonstandards[i++] = *(c - 1);
        }
      }

      c++;
    }
  }

  int il_a, il_b, b_a, b_b;
  linear_fit(&il_a, &il_b, &b_a, &b_b);
#endif
  /**
   * check if we have two input files
   */
  if ((qname == NULL && tname) || (qname && tname == NULL)) {
    RNAplex_cmdline_parser_print_help();
  } else if (qname && tname && !(alignment_mode)) {
    /*free allocated memory of commandline parser*/
    RNAplex_cmdline_parser_free(&args_info);

    if (!fold_constrained) {
      if (access) {
        char *id_s1 = NULL;
        mRNA = fopen(tname, "r");
        if (mRNA == NULL) {
          printf("%s: Wrong target file name\n", tname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        sRNA = fopen(qname, "r");
        if (sRNA == NULL) {
          printf("%s: Wrong quert file name\n", qname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        do {
          /* main loop: continue until end of file */
          if ((line_t = vrna_read_line(mRNA)) == NULL)
            break;

          /*parse line, get id for further accessibility fetching*/
          while ((*line_t == '*') || (*line_t == '\0') || (*line_t == '>')) {
            if (*line_t == '>') {
              if (id_s1) {
                free(id_s1);
                id_s1 = NULL;
              }

              id_s1 = (char *)vrna_alloc(strlen(line_t) + 2);
              (void)sscanf(line_t, "%s", id_s1);
              memmove(id_s1, id_s1 + 1, strlen(id_s1));
              free(line_t);
              line_t = NULL;
            }

            if (line_t)
              free(line_t);

            if ((line_t = vrna_read_line(mRNA)) == NULL)
              break;
          }
          /*append N's to the sequence in order to avoid boundary checking*/
          if ((line_t == NULL) || (strcmp(line_t, "@") == 0)) {
            if (id_s1) {
              free(id_s1);
              id_s1 = NULL;
            }

            break;
          }

          s1 = (char *)vrna_alloc(strlen(line_t) + 1 + 20);
          strcpy(s1, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
          strcat(s1, line_t);
          free(line_t);
          strcat(s1, "NNNNNNNNNN\0");
          int s1_len;
          s1_len = strlen(s1);
          for (l = 0; l < s1_len; l++) {
            s1[l] = toupper(s1[l]);
            if (!noconv && s1[l] == 'T')
              s1[l] = 'U';
          }

          /*read accessibility*/
          int   **access_s1;
          char  *file_s1;
          file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(id_s1) + strlen(access) + 20));
          strcpy(file_s1, access);
          strcat(file_s1, "/");
          strcat(file_s1, id_s1);
          strcat(file_s1, "_openen");
          if (!binaries) {
            access_s1 = read_plfold_i(file_s1, 1, s1_len, verhaeltnis, alignment_length, fast);
          } else {
            strcat(file_s1, "_bin");
            access_s1 = read_plfold_i_bin(file_s1, 1, s1_len, verhaeltnis, alignment_length, fast);
          }

          if (access_s1 == NULL) {
            printf("Accessibility file %s not found or corrupt, look at next target RNA\n", file_s1);
            free(file_s1);
            free(s1);
            free(id_s1);
            id_s1 = NULL;
            continue;
          }

          do {
            char *id_s2 = NULL;
            if ((line_q = vrna_read_line(sRNA)) == NULL)
              break;

            while ((*line_q == '*') || (*line_q == '\0') || (*line_q == '>')) {
              if (*line_q == '>') {
                if (id_s2) {
                  free(id_s2);
                  id_s2 = NULL;
                }

                id_s2 = (char *)vrna_alloc(strlen(line_q) + 2);
                (void)sscanf(line_q, "%s", id_s2);
                memmove(id_s2, id_s2 + 1, strlen(id_s2));
                free(line_q);
                line_q = NULL;
              }

              if (line_q)
                free(line_q);

              if ((line_q = vrna_read_line(sRNA)) == NULL)
                break;
            }
            if ((line_q == NULL) || (strcmp(line_q, "@") == 0)) {
              if (id_s2) {
                free(id_s2);
                id_s2 = NULL;
              }

              break;
            }

            if (!id_s2) {
              free(line_q);
              continue;
            }

            s2 = (char *)vrna_alloc(strlen(line_q) + 1 + 20);
            strcpy(s2, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
            strcat(s2, line_q);
            free(line_q);
            strcat(s2, "NNNNNNNNNN\0");
            int s2_len;
            s2_len = strlen(s2);
            for (l = 0; l < s2_len; l++) {
              s2[l] = toupper(s2[l]);
              if (!noconv && s2[l] == 'T')
                s2[l] = 'U';
            }
            int   **access_s2;
            char  *file_s2;
            file_s2 = (char *)vrna_alloc(sizeof(char) * (strlen(id_s2) + strlen(access) + 20));
            strcpy(file_s2, access);
            strcat(file_s2, "/");
            strcat(file_s2, id_s2);
            strcat(file_s2, "_openen");
            if (!binaries) {
              access_s2 = read_plfold_i(file_s2, 1, s2_len, verhaeltnis, alignment_length, fast);
            } else {
              strcat(file_s2, "_bin");
              access_s2 = read_plfold_i_bin(file_s2, 1, s2_len, verhaeltnis, alignment_length, fast);
            }

            if (access_s2 == NULL) {
              printf("Accessibility file %s not found, look at next target RNA\n", file_s2);
              free(file_s2);
              free(s2);
              free(id_s2);
              id_s2 = NULL;
              continue;
            }

            printf(">%s\n>%s\n", id_s1, id_s2);
            double  begin = BeginTimer();
            Lduplexfold_XS(s1, s2, (const int **)access_s1, (const int **)access_s2, delta, alignment_length, deltaz, fast, il_a, il_b, b_a, b_b);
            float   elapTicks;
            float   elapMilli;
            elapTicks = (EndTimer(begin) - begin);
            elapMilli = elapTicks / 1000;
            /*             printf("Millisecond %.2f\n",elapMilli); */
            free(id_s2);
            free(file_s2);
            free(s2);
            id_s2 = NULL;
            i     = access_s2[0][0];
            while (--i > -1)
              free(access_s2[i]);
            free(access_s2);
          } while (1);
          free(id_s1);
          id_s1 = NULL;
          free(file_s1);
          free(s1);
          rewind(sRNA);
          i = access_s1[0][0];
          while (--i > -1)
            free(access_s1[i]);
          free(access_s1);
        } while (1);
        fclose(mRNA);
        fclose(sRNA);
      } else if (access == NULL) {
        /* t and q are defined, but no accessibility is provided */
        mRNA = fopen(tname, "r");
        if (mRNA == NULL) {
          printf("%s: Wrong target file name\n", tname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        sRNA = fopen(qname, "r");
        if (sRNA == NULL) {
          printf("%s: Wrong query file name\n", qname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        do {
          /* main loop: continue until end of file */
          char *id_s1 = NULL; /* header of the target file  */
          if ((line_t = vrna_read_line(mRNA)) == NULL)
            break;

          /*parse line, get id for further accessibility fetching*/
          while ((*line_t == '*') || (*line_t == '\0') || (*line_t == '>')) {
            if (*line_t == '>') {
              if (id_s1) {
                /* in case we have two header the one after the other */
                free(id_s1);     /* free the old header, a put the new one instead */
                id_s1 = NULL;
              }

              id_s1 = (char *)vrna_alloc(strlen(line_t) + 2);
              (void)sscanf(line_t, "%s", id_s1);
              memmove(id_s1, id_s1 + 1, strlen(id_s1));
              free(line_t);
              line_t = NULL;
            }

            if (line_t)
              free(line_t);

            if ((line_t = vrna_read_line(mRNA)) == NULL)
              break;
          }
          /*append N's to the sequence in order to avoid boundary checking*/
          if ((line_t == NULL) || (strcmp(line_t, "@") == 0)) {
            if (id_s1) {
              free(id_s1);
              id_s1 = NULL;
            }

            break;
          }

          s1 = (char *)vrna_alloc(strlen(line_t) + 1 + 20);
          strcpy(s1, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
          strcat(s1, line_t);
          free(line_t);
          strcat(s1, "NNNNNNNNNN\0");
          int s1_len;
          s1_len = strlen(s1);
          for (l = 0; l < s1_len; l++) {
            s1[l] = toupper(s1[l]);
            if (!noconv && s1[l] == 'T')
              s1[l] = 'U';
          }
          do {
            /*read sRNA files*/
            char *id_s2 = NULL;
            if ((line_q = vrna_read_line(sRNA)) == NULL)
              break;

            while ((*line_q == '*') || (*line_q == '\0') || (*line_q == '>')) {
              if (*line_q == '>') {
                if (id_s2) {
                  free(id_s2);
                  id_s2 = NULL;
                }

                id_s2 = (char *)vrna_alloc(strlen(line_q) + 2);
                (void)sscanf(line_q, "%s", id_s2);
                memmove(id_s2, id_s2 + 1, strlen(id_s2));
                free(line_q);
                line_q = NULL;
              }

              if (line_q)
                free(line_q);

              if ((line_q = vrna_read_line(sRNA)) == NULL)
                break;
            }
            if ((line_q == NULL) || (strcmp(line_q, "@") == 0)) {
              if (id_s2)
                free(id_s2);

              break;
            }

            if (!id_s2) {
              free(line_q);
              continue;
            }

            s2 = (char *)vrna_alloc(strlen(line_q) + 1 + 20);
            strcpy(s2, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
            strcat(s2, line_q);
            free(line_q);
            strcat(s2, "NNNNNNNNNN\0");
            int s2_len;
            s2_len = strlen(s2);
            for (l = 0; l < s2_len; l++) {
              s2[l] = toupper(s2[l]);
              if (!noconv && s2[l] == 'T')
                s2[l] = 'U';
            }
            printf(">%s\n>%s\n", id_s1, id_s2);
            /*             double begin = BeginTimer(); */
            Lduplexfold(s1, s2, delta, extension_cost, alignment_length, deltaz, fast, il_a, il_b, b_a, b_b);
            /* float elapTicks; */
            /* float elapMilli; */
            /* elapTicks = (EndTimer(begin) - begin); */
            /* elapMilli = elapTicks/1000; */
            /* printf("Millisecond %.2f\n",elapMilli); */
            /* printf("\n"); */
            free(id_s2);
            id_s2 = NULL;
            free(s2);
          } while (1);
          free(id_s1);
          id_s1 = NULL;
          free(s1);
          rewind(sRNA);
        } while (1);
        fclose(mRNA);
        fclose(sRNA);
      }
    } else {
      if (access) {
        char *id_s1 = NULL;
        mRNA = fopen(tname, "r");
        if (mRNA == NULL) {
          printf("%s: Wrong target file name\n", tname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        sRNA = fopen(qname, "r");
        if (sRNA == NULL) {
          printf("%s: Wrong query file name\n", qname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        do {
          /* main loop: continue until end of file */
          if ((line_t = vrna_read_line(mRNA)) == NULL)
            break;

          /*parse line, get id for further accessibility fetching*/
          while ((*line_t == '*') || (*line_t == '\0') || (*line_t == '>')) {
            if (*line_t == '>') {
              if (id_s1) {
                /* in case we have two header the one after the other */
                free(id_s1);     /* free the old header, a put the new one instead */
                id_s1 = NULL;
              }

              id_s1 = (char *)vrna_alloc(strlen(line_t) + 2);
              (void)sscanf(line_t, "%s", id_s1);
              memmove(id_s1, id_s1 + 1, strlen(id_s1));
              free(line_t);
              line_t = NULL;
            }

            if (line_t)
              free(line_t);

            if ((line_t = vrna_read_line(mRNA)) == NULL)
              break;
          }
          /*append N's to the sequence in order to avoid boundary checking*/
          if ((line_t == NULL) || (strcmp(line_t, "@") == 0)) {
            if (id_s1) {
              free(id_s1);
              id_s1 = NULL;
            }

            break;
          }

          s1 = (char *)vrna_alloc(strlen(line_t) + 1 + 20);
          strcpy(s1, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
          strcat(s1, line_t);
          free(line_t);
          strcat(s1, "NNNNNNNNNN\0");
          int s1_len;
          s1_len = strlen(s1);
          for (l = 0; l < s1_len; l++) {
            s1[l] = toupper(s1[l]);
            if (!noconv && s1[l] == 'T')
              s1[l] = 'U';
          }

          /*read accessibility*/
          int   **access_s1;
          char  *file_s1;
          file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(id_s1) + strlen(access) + 20));
          strcpy(file_s1, access);
          strcat(file_s1, "/");
          strcat(file_s1, id_s1);
          strcat(file_s1, "_openen");
          if (!binaries) {
            access_s1 = read_plfold_i(file_s1, 1, s1_len, verhaeltnis, alignment_length, fast);
          } else {
            strcat(file_s1, "_bin");
            access_s1 = read_plfold_i_bin(file_s1, 1, s1_len, verhaeltnis, alignment_length, fast);
          }

          if (access_s1 == NULL) {
            printf("Accessibility file %s not found, look at next target RNA\n", file_s1);
            free(file_s1);
            free(s1);
            free(id_s1);
            id_s1 = NULL;
            continue;
          }

          do {
            char *id_s2 = NULL;
            /*read sRNA files*/
            if ((line_q = vrna_read_line(sRNA)) == NULL)
              break;

            while ((*line_q == '*') || (*line_q == '\0') || (*line_q == '>')) {
              if (*line_q == '>') {
                if (id_s2) {
                  free(id_s2);
                  id_s2 = NULL;
                }

                id_s2 = (char *)vrna_alloc(strlen(line_q) + 2);
                (void)sscanf(line_q, "%s", id_s2);
                memmove(id_s2, id_s2 + 1, strlen(id_s2));
                free(line_q);
                line_q = NULL;
              }

              if (line_q)
                free(line_q);

              if ((line_q = vrna_read_line(sRNA)) == NULL)
                break;
            }
            if ((line_q == NULL) || (strcmp(line_q, "@") == 0)) {
              if (id_s2)
                free(id_s2);

              break;
            }

            /* if ((line_t ==NULL) || (strcmp(line_t, "@") == 0)) break; */
            s2 = (char *)vrna_alloc(strlen(line_q) + 1 + 20);
            strcpy(s2, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
            strcat(s2, line_q);
            free(line_q);
            strcat(s2, "NNNNNNNNNN\0");
            int s2_len;
            s2_len = strlen(s2);
            for (l = 0; l < s2_len; l++) {
              s2[l] = toupper(s2[l]);
              if (!noconv && s2[l] == 'T')
                s2[l] = 'U';
            }
            structure = (char *)vrna_alloc((unsigned)s2_len + 1);
            cstruc    = vrna_read_line(sRNA);
            if (cstruc != NULL) {
              int dn3 = strlen(cstruc) - (s2_len - 20);
              strcpy(structure, "..........");
              strncat(structure, cstruc, s2_len - 20);
              if (dn3 >= 0) {
                strcat(structure, "..........\0");
              } else {
                while (dn3++)
                  strcat(structure, ".");
                strcat(structure, "\0");
              }

              free(cstruc);
            } else {
              vrna_message_warning("constraints missing");
            }

            int   a = strchr(structure, '|') - structure;
            int   b = strrchr(structure, '|') - structure;
            if (alignment_length < b - a + 1)
              vrna_message_error("Maximal duplex length (-l option) is smaller than constraint on the structures\n. Please adjust the -l option accordingly\n");

            int   **access_s2;
            char  *file_s2;
            file_s2 = (char *)vrna_alloc(sizeof(char) * (strlen(id_s2) + strlen(access) + 20));
            strcpy(file_s2, access);
            strcat(file_s2, "/");
            strcat(file_s2, id_s2);
            strcat(file_s2, "_openen");
            if (!binaries) {
              access_s2 = read_plfold_i(file_s2, 1, s2_len, verhaeltnis, alignment_length, fast);
            } else {
              strcat(file_s2, "_bin");
              access_s2 = read_plfold_i_bin(file_s2, 1, s2_len, verhaeltnis, alignment_length, fast);
            }

            if (access_s2 == NULL) {
              printf("Accessibility file %s not found, look at next target RNA\n", file_s2);
              free(file_s2);
              free(s2);
              free(id_s2);
              continue;
            }

            printf(">%s\n>%s\n", id_s1, id_s2);
            Lduplexfold_CXS(s1, s2, (const int **)access_s1, (const int **)access_s2, delta, alignment_length, deltaz, fast, structure, il_a, il_b, b_a, b_b);/* , target_dead, query_dead); */
            free(id_s2);
            free(file_s2);
            free(s2);
            i = access_s2[0][0];
            while (--i > -1)
              free(access_s2[i]);
            free(access_s2);
            free(structure);
          } while (1);
          free(id_s1);
          id_s1 = NULL;
          free(file_s1);
          free(s1);
          rewind(sRNA);
          i = access_s1[0][0];
          while (--i > -1)
            free(access_s1[i]);
          free(access_s1);
        } while (1);
        fclose(mRNA);
        fclose(sRNA);
      } else if (access == NULL) {
        /* t and q are defined, but no accessibility is provided */
        char *id_s1 = NULL;
        mRNA = fopen(tname, "r");
        if (mRNA == NULL) {
          printf("%s: Wrong target file name\n", tname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        sRNA = fopen(qname, "r");
        if (sRNA == NULL) {
          printf("%s: Wrong query file name\n", qname);
          RNAplex_cmdline_parser_free(&args_info);
          return 0;
        }

        do {
          /* main loop: continue until end of file */
          if ((line_t = vrna_read_line(mRNA)) == NULL)
            break;

          /*parse line, get id for further accessibility fetching*/
          while ((*line_t == '*') || (*line_t == '\0') || (*line_t == '>')) {
            if (*line_t == '>') {
              if (id_s1) {
                /* in case we have two header the one after the other */
                free(id_s1);     /* free the old header, a put the new one instead */
                id_s1 = NULL;
              }

              id_s1 = (char *)vrna_alloc(strlen(line_t) + 2);
              (void)sscanf(line_t, "%s", id_s1);
              memmove(id_s1, id_s1 + 1, strlen(id_s1));
              free(line_t);
              line_t = NULL;
            }

            if (line_t)
              free(line_t);

            if ((line_t = vrna_read_line(mRNA)) == NULL)
              break;
          }
          /*append N's to the sequence in order to avoid boundary checking*/
          /* if ((line_t ==NULL) || (strcmp(line_t, "@") == 0)) break; */
          if ((line_t == NULL) || (strcmp(line_t, "@") == 0)) {
            if (id_s1) {
              free(id_s1);
              id_s1 = NULL;
            }

            break;
          }

          s1 = (char *)vrna_alloc(strlen(line_t) + 1 + 20);
          strcpy(s1, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
          strcat(s1, line_t);
          free(line_t);
          strcat(s1, "NNNNNNNNNN\0");
          int s1_len;
          s1_len = strlen(s1);
          for (l = 0; l < s1_len; l++) {
            s1[l] = toupper(s1[l]);
            if (!noconv && s1[l] == 'T')
              s1[l] = 'U';
          }
          do {
            char *id_s2 = NULL;
            /*read sRNA files*/
            if ((line_q = vrna_read_line(sRNA)) == NULL)
              break;

            while ((*line_q == '*') || (*line_q == '\0') || (*line_q == '>')) {
              if (*line_q == '>') {
                if (id_s2) {
                  free(id_s2);
                  id_s2 = NULL;
                }

                id_s2 = (char *)vrna_alloc(strlen(line_q) + 2);
                (void)sscanf(line_q, "%s", id_s2);
                memmove(id_s2, id_s2 + 1, strlen(id_s2));
                free(line_q);
                line_q = NULL;
              }

              if (line_q)
                free(line_q);

              if ((line_q = vrna_read_line(sRNA)) == NULL)
                break;
            }
            /* if ((line_t ==NULL) || (strcmp(line_t, "@") == 0)) break; */
            if ((line_q == NULL) || (strcmp(line_q, "@") == 0)) {
              if (id_s2)
                free(id_s2);

              break;
            }

            s2 = (char *)vrna_alloc(strlen(line_q) + 1 + 20);
            strcpy(s2, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
            strcat(s2, line_q);
            free(line_q);
            strcat(s2, "NNNNNNNNNN\0");
            int s2_len;
            s2_len = strlen(s2);
            for (l = 0; l < s2_len; l++) {
              s2[l] = toupper(s2[l]);
              if (!noconv && s2[l] == 'T')
                s2[l] = 'U';
            }
            structure = (char *)vrna_alloc((unsigned)s2_len + 1);
            cstruc    = vrna_read_line(sRNA);
            if (cstruc != NULL) {
              int dn3 = strlen(cstruc) - (s2_len - 20);
              strcpy(structure, "..........");
              strncat(structure, cstruc, s2_len - 20);
              if (dn3 >= 0) {
                strcat(structure, "..........\0");
              } else {
                while (dn3++)
                  strcat(structure, ".");
                strcat(structure, "\0");
              }

              free(cstruc);
            } else {
              vrna_message_warning("constraints missing");
            }

            int     a = strchr(structure, '|') - structure;
            int     b = strrchr(structure, '|') - structure;
            if (alignment_length < b - a + 1)
              vrna_message_error("Maximal duplex length (-l option) is smaller than constraint on the structures\n. Please adjust the -l option accordingly\n");

            printf(">%s\n>%s\n", id_s1, id_s2);
            double  begin = BeginTimer();
            Lduplexfold_C(s1, s2, delta, extension_cost, alignment_length, deltaz, fast, structure, il_a, il_b, b_a, b_b);
            float   elapTicks;
            float   elapMilli;
            elapTicks = EndTimer(begin);
            elapMilli = elapTicks / 1000;
            /*             printf("Millisecond %.2f\n",elapMilli); */
            free(id_s2);
            free(s2);
            free(structure);
          } while (1);
          free(id_s1);
          id_s1 = NULL;
          free(s1);
          rewind(sRNA);
        } while (1);
        fclose(mRNA);
        fclose(sRNA);
      }
    }
  } else if (!qname && !tname && !(alignment_mode)) {
    istty = isatty(fileno(stdout)) && isatty(fileno(stdin));

    if ((fold_constrained) && (istty)) {
      printf("Input constraints using the following notation:\n");
      printf("| : paired with another base\n");
      printf(". : no constraint at all\n");
      printf("x : base must not pair\n");
    }

    do {
      /* main loop: continue until end of file */
      /* duplexT mfe, *subopt         */
      char *id_s1 = NULL, *id_s2 = NULL;
      if (istty) {
        printf("\nInput two sequences (one line each); @ to quit\n");
        printf("%s\n", scale);
      }

      fname[0] = '\0';

      if ((line = vrna_read_line(stdin)) == NULL)
        break;

      /* skip empty lines, comment lines, name lines */
      while ((*line == '*') || (*line == '\0') || (*line == '>')) {
        printf("%s\n", line);
        if (*line == '>') {
          id_s1 = (char *)vrna_alloc(strlen(line) + 2);
          (void)sscanf(line, "%s", id_s1);
          memmove(id_s1, id_s1 + 1, strlen(id_s1));
        }

        free(line);
        if ((line = vrna_read_line(stdin)) == NULL) {
          free(id_s1);
          break;
        }
      }
      if ((line == NULL) || (strcmp(line, "@") == 0))
        break;

      s1 = (char *)vrna_alloc(strlen(line) + 1 + 20);
      strcpy(s1, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
      strcat(s1, line);
      free(line);
      strcat(s1, "NNNNNNNNNN\0");
      if ((line = vrna_read_line(stdin)) == NULL)
        break;

      /* skip comment lines and get filenames */

      while ((*line == '*') || (*line == '\0') || (*line == '>')) {
        printf("%s\n", line);
        if (*line == '>') {
          id_s2 = (char *)vrna_alloc(strlen(line) + 2);
          (void)sscanf(line, "%s", id_s2);
          memmove(id_s2, id_s2 + 1, strlen(id_s2));
        }

        free(line);
        if ((line = vrna_read_line(stdin)) == NULL) {
          free(id_s2);
          break;
        }
      }
      if ((line == NULL) || (strcmp(line, "@") == 0))
        break;

      s2 = (char *)vrna_alloc(strlen(line) + 1 + 20);
      strcpy(s2, "NNNNNNNNNN"); /*add NNNNNNNNNN to avoid boundary check*/
      strcat(s2, line);
      free(line);
      strcat(s2, "NNNNNNNNNN\0");

      int n1  = strlen(s1);
      int n2  = strlen(s2);


      structure = (char *)vrna_alloc((unsigned)n2 + 1);
      if (fold_constrained) {
        cstruc = vrna_read_line(stdin);
        if (cstruc != NULL && (cstruc[0] == '>')) {
          int dn3 = strlen(cstruc) - (n2 - 20);
          strcpy(structure, "..........");
          strncat(structure, cstruc, n2 - 20);
          if (dn3 >= 0) {
            strcat(structure, "..........\0");
          } else {
            while (dn3++)
              strcat(structure, ".");
            strcat(structure, "\0");
          }

          free(cstruc);
        } else {
          vrna_message_warning("constraints missing");
        }
      }

      for (l = 0; l < n1; l++) {
        s1[l] = toupper(s1[l]);
        if (!noconv && s1[l] == 'T')
          s1[l] = 'U';
      }
      for (l = 0; l < n2; l++) {
        s2[l] = toupper(s2[l]);
        if (!noconv && s2[l] == 'T')
          s2[l] = 'U';
      }
      if (istty)
        printf("lengths = %d,%d\n", (int)strlen(s1), (int)strlen(s2));

      if (alignment_length == 0)
        alignment_length = 40;

      if (access == NULL) {
        if (!fold_constrained) {
          Lduplexfold(s1, s2, delta, extension_cost, alignment_length, deltaz, fast, il_a, il_b, b_a, b_b);
        } else {
          int a = strchr(structure, '|') - structure;
          int b = strrchr(structure, '|') - structure;
          if (alignment_length < b - a + 1)
            vrna_message_error("Maximal duplex length (-l option) is smaller than constraint on the structures\n. Please adjust the -l option accordingly\n");

          Lduplexfold_C(s1, s2, delta, extension_cost, alignment_length, deltaz, fast, structure, il_a, il_b, b_a, b_b);
        }
      } else {
        int   **access_s1, **access_s2;
        char  *file_s1, *file_s2;
        int   s1_len, s2_len;
        s1_len  = strlen(s1);
        s2_len  = strlen(s2);
        if (!(id_s1 && id_s2))
          vrna_message_error("The fasta files has no header information..., cant fetch accessibility file\n");

        file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(id_s1) + strlen(access) + 20));
        file_s2 = (char *)vrna_alloc(sizeof(char) * (strlen(id_s2) + strlen(access) + 20));
        strcpy(file_s1, access);
        strcpy(file_s2, access);
        strcat(file_s1, "/");
        strcat(file_s2, "/");
        strcat(file_s1, id_s1);
        strcat(file_s2, id_s2);
        strcat(file_s1, "_openen");
        strcat(file_s2, "_openen");
        if (!binaries) {
          access_s1 = read_plfold_i(file_s1, 1, s1_len, verhaeltnis, alignment_length, fast);
        } else {
          strcat(file_s1, "_bin");
          access_s1 = read_plfold_i_bin(file_s1, 1, s1_len, verhaeltnis, alignment_length, fast);
        }

        if (access_s1 == NULL) {
          free(file_s1);
          free(s1);
          free(s2);
          free(file_s2);
          free(id_s1);
          free(id_s2);
          continue;
        }

        if (!binaries) {
          access_s2 = read_plfold_i(file_s2, 1, s2_len, verhaeltnis, alignment_length, fast);
        } else {
          strcat(file_s2, "_bin");
          access_s2 = read_plfold_i_bin(file_s2, 1, s2_len, verhaeltnis, alignment_length, fast);
        }

        if (access_s2 == NULL) {
          free(access_s1);
          free(file_s1);
          free(s1);
          free(s2);
          free(file_s2);
          free(id_s1);
          free(id_s2);
          continue;
        }

        if (!fold_constrained) {
          Lduplexfold_XS(s1, s2, (const int **)access_s1, (const int **)access_s2, delta, alignment_length, deltaz, fast, il_a, il_b, b_a, b_b);/* , target_dead, query_dead); */
        } else {
          int a = strchr(structure, '|') - structure;
          int b = strrchr(structure, '|') - structure;
          if (alignment_length < b - a + 1)
            vrna_message_error("Maximal duplex length (-l option) is smaller than constraint on the structures\n. Please adjust the -l option accordingly\n");

          Lduplexfold_CXS(s1, s2, (const int **)access_s1, (const int **)access_s2, delta, alignment_length, deltaz, fast, structure, il_a, il_b, b_a, b_b);/* , target_dead, query_dead); */
        }

        i = access_s1[0][0];
        while (--i > -1)
          free(access_s1[i]);
        i = access_s2[0][0];
        while (--i > -1)
          free(access_s2[i]);
        free(access_s1);
        free(access_s2);
        free(file_s1);
        free(file_s2);
        free(id_s1);
        id_s1 = NULL;
        free(id_s2);
        id_s2 = NULL;
      }

      free(s1);
      free(s2);
      s1  = NULL;
      s2  = NULL;
      free(structure);
      printf("\n");
      if (id_s1)
        free(id_s1);

      if (id_s2)
        free(id_s2);
    } while (1);
    if (s1)
      free(s1);

    if (s2)
      free(s2);

    /* if(id_s1){free(id_s1);}if(id_s2){free(id_s2);} */
  } else if (qname && tname && alignment_mode) {
    int n_seq, n_seq2;
    mRNA = fopen(tname, "r");
    if (mRNA == NULL) {
      printf("%s: Wrong target file name\n", tname);
      RNAplex_cmdline_parser_free(&args_info);
      return 0;
    }

    sRNA = fopen(qname, "r");
    if (sRNA == NULL) {
      printf("%s: Wrong query file name\n", qname);
      RNAplex_cmdline_parser_free(&args_info);
      return 0;
    }

    n_seq   = read_clustal(mRNA, temp1, names1);
    n_seq2  = read_clustal(sRNA, temp2, names2);
    fclose(mRNA);
    fclose(sRNA);
    if (n_seq != n_seq2) {
      for (i = 0; temp1[i]; i++) {
        free(temp1[i]);
        free(temp2[i]);
      }
      vrna_message_error("unequal number of seqs in alignments");
    }

    for (i = 0; temp1[i]; i++) {
      AS1[i]  = (char *)vrna_alloc((strlen(temp1[i]) + 21) * sizeof(char));
      AS2[i]  = (char *)vrna_alloc((strlen(temp2[i]) + 21) * sizeof(char));
      strcpy(AS1[i], "NNNNNNNNNN");
      strcpy(AS2[i], "NNNNNNNNNN");
      strcat(AS1[i], temp1[i]);
      strcat(AS2[i], temp2[i]);
      strcat(AS1[i], "NNNNNNNNNN\0");
      strcat(AS2[i], "NNNNNNNNNN\0");
    }
    for (i = 0; temp1[i]; i++) {
      free(temp1[i]);
      free(temp2[i]);
    }
    AS1[n_seq]  = NULL;
    AS2[n_seq]  = NULL;


    if (access == NULL) {
      aliLduplexfold((const char **)AS1, (const char **)AS2, n_seq * delta, extension_cost, alignment_length, deltaz, fast, il_a, il_b, b_a, b_b);
    } else {
      int **target_access = NULL, **query_access = NULL;
      target_access = average_accessibility_target(names1, AS1, n_seq, access, verhaeltnis, alignment_length, binaries, fast); /* get averaged accessibility for alignments */
      query_access  = average_accessibility_target(names2, AS2, n_seq, access, verhaeltnis, alignment_length, binaries, fast);
      if (!(target_access && query_access)) {
        for (i = 0; AS1[i]; i++) {
          free(AS1[i]);
          free(names1[i]);
          free(AS2[i]);
          free(names2[i]);
        }
        RNAplex_cmdline_parser_free(&args_info);
        return 0;
      }

      aliLduplexfold_XS((const char **)AS1, (const char **)AS2, (const int **)target_access, (const int **)query_access, n_seq * delta, alignment_length, deltaz, fast, il_a, il_b, b_a, b_b);
      if (!(target_access == NULL)) {
        i = target_access[0][0];
        while (--i > -1)
          free(target_access[i]);
        free(target_access);
      }

      if (!(query_access == NULL)) {
        i = query_access[0][0];
        while (--i > -1)
          free(query_access[i]);
        free(query_access);
      }
    }

    for (i = 0; AS1[i]; i++) {
      free(AS1[i]);
      free(names1[i]);
      free(AS2[i]);
      free(names2[i]);
    }
  }

  if (access) {
    free(access);
    access = NULL;
  }

  if (qname) {
    free(tname);
    access = NULL;
  }

  if (tname) {
    free(qname);
    access = NULL;
  }

  RNAplex_cmdline_parser_free(&args_info);
  return 0;
}


#if 0
static int
print_struc(duplexT const *dup)
{
  /*this function print out the structure of a hybridization site
   * and return the value from which one can begin to look for the next
   * non-overlappig hybrid*/

  int   l1;
  float energy  = dup->energy * 0.01;
  int   init    = dup->end;

  l1 = strchr(dup->structure, '&') - dup->structure;
  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, init - l1 - 1,
         init - 2, dup->j + l1 + 2 - strlen(dup->structure), dup->j, energy);
  return init - l1 - 1;
}


#endif

static int **
read_plfold_i(char      *fname,
              const int beg,
              const int end,
              double    verhaeltnis,
              const int length,
              int       fast)
{
  double  begin = BeginTimer();
  FILE    *in   = fopen(fname, "r");

  if (in == NULL) {
    vrna_message_warning("File ' %s ' open error", fname);
    return NULL;
  }

  int   i, j;
  int   **access;
  int   offset, temp;
  temp    = 0;
  offset  = 0;
  int   seq_pos;
  int   beg_r, end_r;
  beg_r = beg;
  end_r = end;
  char  tmp[2048] = {
    0x0
  };
  /* char *ptr; */
  if (fgets(tmp, sizeof(tmp), in) == 0)
    perror("Empty File");

  if (strchr(tmp, '>')) {
    vrna_message_warning("file %s is not in RNAplfold format", fname);
    return NULL;
  }

  if (fgets(tmp, sizeof(tmp), in) == 0) {
    vrna_message_warning("No accessibility data");
    return NULL;
  }

  int dim_x;
  dim_x = get_max_u(tmp, '\t');
  if (length > dim_x && fast == 0) {
    printf("Interaction length %d is larger than the length of the largest region %d \nfor which the opening energy was computed (-u parameter of RNAplfold)\n", length, dim_x);
    printf("Please recompute your profiles with a larger -u or set -l to a smaller interaction length\n");
    return NULL;
  }

  access = (int **)vrna_alloc(sizeof(int *) * (dim_x + 2));
  for (i = 0; i < dim_x + 2; i++)
    access[i] = (int *)vrna_alloc(sizeof(int) * (end_r - beg_r + 1)); /* normally +1 */
  for (i = 0; i < end_r - beg_r + 1; i++)
    for (j = 0; j < dim_x + 2; j++)
      access[j][i] = INF;
  access[0][0] = dim_x + 2;
  while (fgets(tmp, sizeof(tmp), in) != 0 && --end_r > 10) {
    /* read a record, before we have --end_r > 10*/
    float n;
    /* int i; */
    /* int u; */
    beg_r--;
    if (beg_r < 1) {
      if (sscanf(tmp, "%d\t%n", &seq_pos, &temp) == 1) {
        offset += temp;
        int count;
        count = 1;
        while (sscanf(tmp + offset, "%f%n", &n, &temp) == 1) {
          offset += temp;
          /* seq_pos - beg allows one to get the accessibility right */
          access[count][seq_pos - beg + 11] = (int)rint(100 * n); /* round the number */
          access[count][seq_pos - beg + 11] *= verhaeltnis;       /* 10 stands here for the number of nucleotides */
          count++;
        }
        offset = 0;
      }
    }
  }
  if (end_r > 20) {
    printf("Accessibility files contains %d less entries than expected based on the sequence length\n", end_r - 20);
    printf("Please recompute your profiles so that profile length and sequence length match\n");
    return NULL;
  }

  fclose(in);
  float elapTicks;
  float elapMilli;
  elapTicks = (EndTimer(begin) - begin);
  elapMilli = elapTicks / 1000;
  return access;
}


static int
convert_plfold_i(char *fname)
{
  int   i;
  FILE  *in = fopen(fname, "r");

  if (in == NULL) {
    vrna_message_warning("File ' %s ' open error", fname);
    return -1;
  }

  char tmp[2048] = {
    0x0
  };
  if (fgets(tmp, sizeof(tmp), in) == 0)
    perror("Empty File");

  if (strchr(tmp, '>')) {
    vrna_message_warning("file %s is not in RNAplfold format", fname);
    return -1;
  }

  if (fgets(tmp, sizeof(tmp), in) == 0) {
    vrna_message_warning("No accessibility data");
    return -1;
  }

  int   u_length;
  u_length = get_max_u(tmp, '\t'); /* get the x dimension */
  char  c;
  int   length = 0;
  while ((c = fgetc(in)) != EOF)
    if (c == '\n')
      length++;

  fclose(in);
  int   **access = read_plfold_i(fname, 1, length + 20, 1, u_length, 2);
  char  *outname;
  outname = (char *)vrna_alloc((strlen(fname) + 5) * sizeof(char));
  strcpy(outname, fname);
  strcat(outname, "_bin");
  FILE  *fp = fopen(outname, "wb");
  int   p[1];
  p[0] = 1000000;
  for (i = 0; i < u_length + 2; i++) {
    fwrite(access[i], sizeof(int), length + 20, fp);
    free(access[i]);
  }
  fseek(fp, 0, SEEK_SET);
  fwrite(&u_length, sizeof(int), 1, fp);
  fwrite(&length, sizeof(int), 1, fp);
  free(outname);
  fclose(fp);
  free(access);
  return 1;
}


static int **
read_plfold_i_bin(char      *fname,
                  const int beg,
                  const int end,
                  double    verhaeltnis,
                  const int length,
                  int       fast)
{
  double  begin = BeginTimer();
  FILE    *fp   = fopen(fname, "rb");
  int     seqlength;

  if (fp == NULL) {
    vrna_message_warning("File ' %s ' open error", fname);
    return NULL;
  }

  int *first_line;
  first_line = (int *)vrna_alloc(sizeof(int) * (end - beg + 1));              /* check length of the line LOOK at read_plfold_i */
  if (!fread(first_line, sizeof(int), (end - beg) + 1, fp)) {
    vrna_message_warning("Problem reading size of profile from '%s'", fname); /* get the value of the u option */
    return NULL;
  }

  int lim_x;
  lim_x     = first_line[0];
  seqlength = first_line[1];                                  /* length of the sequence RNAplfold was ran on. */
  if (length > lim_x && fast == 0) {
    printf("Interaction length %d is larger than the length of the largest region %d \nfor which the opening energy was computed (-u parameter of RNAplfold)\n", length, lim_x);
    printf("Please recompute your profiles with a larger -u or set -l to a smaller interaction length\n");
    return NULL;
  }

  fseek(fp, 0, SEEK_SET);                                   /* set pointer to the begining of the file */
  int       **access;                                       /* here we store the access  values */
  access = (int **)vrna_alloc(sizeof(int *) * (lim_x + 1)); /* !!!! lim_x+1 -> lim_x+2 */
  int       count;
  long int  position;
  for (count = 0; count < lim_x + 1; count++) {
    /* read file */
    access[count] = (int *)vrna_alloc(sizeof(int) * (end - beg + 1)); /* declare array length */
    /* now we should be sure to read the right position */
    /* we first need a seek to the right position */
    /* 0->9 line begin, + begin  */
    /* here we need information about the total line length..., we can get this if this is saved by RNAplfold... */
    position = ftell(fp);
    fseek(fp, (beg - 1) * sizeof(int), SEEK_CUR);                 /* go to the desired position, note the 10 offset */
    position = ftell(fp);
    if (!fread(access[count], sizeof(int), (end - beg) + 1, fp))  /* read the needed number of accessibility values */
      printf("File '%s' is corrupted \n", fname);

    position = ftell(fp);
    fseek(fp, (seqlength - end + 20) * sizeof(int), SEEK_CUR); /* place to the begining of the next file */
    position = ftell(fp);
  }
  fseek(fp, 0, SEEK_SET);
  access[0][0]  = lim_x;
  access[0][0]  += 1; /* !!!!!!!!!!!!!!!!!!!put 2 in case of problem */
  fclose(fp);         /* close file */
  float elapTicks;
  float elapMilli;
  free(first_line);
  elapTicks = (EndTimer(begin) - begin);
  elapMilli = elapTicks / 1000;
  /* printf("read_plfold_i_bin Millisecond %.2f\n",elapMilli); */ /* return time */
  return access; /* return access */
}


static int
get_max_u(const char  *s,
          char        delim)
{
  char  *pch;
  int   total;

  total = 0;
  pch   = strchr(s, delim);
  total++;
  while (pch != NULL) {
    pch = strchr(pch + 1, delim);
    total++;
  }
  return total - 2;
}


static int **
average_accessibility_target(char       **names,
                             char       **ALN,
                             int        number,
                             char       *access,
                             double     verhaeltnis,
                             const int  alignment_length,
                             int        binaries,
                             int        fast)
{
  int           i;
  int           ***master_access  = NULL;           /* contains the accessibility arrays for different */
  int           **average_access  = NULL;           /* contains the averaged accessibility */
  int           aln_size          = strlen(ALN[0]); /* aln size -> define size of the averaged accessibility array */
  int           *index;                             /* contains the index used for navigating inside the alignments */
  long long int begin, end;                         /* contains the begin and end region to read the accessibility */

  index = (int *)vrna_alloc(sizeof(int) * number);
  for (i = 0; i < number; i++)
    index[i] = 1;
  master_access = (int ***)vrna_alloc(sizeof(int **) * number);
  int   dim_x; /* contains the minimal of the maximum u length for all sequences of the alignment */
  dim_x = INF;
  int   u;
  int   location_flag; /* location flag checks if all line contains a location flag or not */
  /* if some line contains location information and other not, return a warning but keep going; */
  /* The location flag will be set to TRUE if the following pattern is found */
  /*  sequence-name/[0-9]+-[0-9]+  */
  /*  |----NAME----||BEG|  |END| */
  char  bla[255];
  if (sscanf(names[0], "%255[^/]/%lld-%lld", bla, &begin, &end) == 3) /* initialize location_flag; */
    location_flag = 1;
  else
    location_flag = 0;

  char *file_s1 = NULL;
  for (i = 0; i < number; i++) {
    /*  be careful!!!! Name should contain all characters from begin till the "/" character */
    /* char *s1; */
    file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(names[i]) + strlen(access) + 20));
    begin   = 1;
    int sequence_length = get_sequence_length_from_alignment(ALN[i]);
    end = sequence_length;
    if (sscanf(names[i], "%255[^/]/%lld-%lld", bla, &begin, &end) == 3) {
      /* if subsequence and range do not have the same length, stop RNAplex. */
      if (end - begin + 1 != sequence_length - 20) {
        printf("Your range %lld %lld in sequence %s does not correspond to the sequence length %d\n", begin, end, names[i], sequence_length - 20);
        printf("Please check your input alignments and rerun RNAplex");
        int a = 0;
        if (master_access != NULL) {
          for (a = 0; a < i; a++) {
            int j;
            j = master_access[a][0][0];
            while (--j > -1)
              free(master_access[a][j]);
            free(master_access[a]);
          }
        }

        free(master_access);
        free(index);
        free(file_s1);
        return average_access;
      }

      end += 20; /* add 20 to the end, in order to take the N's into account */
      if (location_flag == 0) {
        vrna_message_warning("\n!! Line %d in your target alignment contains location information\n"
                             "while line %d did not. PLEASE CHECK your alignments!!\n"
                             "RNAplex will continue the target search.",
                             i + 1, i);
      }

      location_flag = 1;
      strcpy(file_s1, access);
      strcat(file_s1, "/");
      strcat(file_s1, bla);
    } else {
      if (location_flag == 1) {
        vrna_message_warning("\n!! Line %d in your target alignment does not contain location information\n"
                             "while line %d in your target alignment did. PLEASE CHECK your alignments!!\n"
                             "RNAplex will continue the target search.",
                             i + 1, i);
      }

      location_flag = 0;
      strcpy(file_s1, access);
      strcat(file_s1, "/");
      strcat(file_s1, names[i]);
    }

    strcat(file_s1, "_openen");
    if (!binaries) {
      master_access[i] = read_plfold_i(file_s1, begin, end, verhaeltnis, alignment_length, fast); /* read */
    } else {
      strcat(file_s1, "_bin");
      master_access[i] = read_plfold_i_bin(file_s1, begin, end, verhaeltnis, alignment_length, fast); /* read */
    }

    free(file_s1);
    file_s1 = NULL;
    dim_x   = MIN2(dim_x, master_access[i][0][0]);
  }
  average_access = (int **)vrna_alloc(sizeof(int *) * (dim_x));
  for (i = 0; i < dim_x; i++)
    average_access[i] = (int *)vrna_alloc(sizeof(int) * (aln_size + 9));  /* average_access is of the length of the alignment */
  /* while master_access is of the length of the sequences. */
  average_access[0][0] = dim_x;
  for (i = 1; i < aln_size; i++) {
    /* go through all aln position */
    int j;
    int count_j = number;
    for (j = 0; j < number; j++) {
      /* go through all aln members */
      if (!(ALN[j][i - 1] == '-')) {
        /* if we have a gap, skip this position */
        for (u = 1; u < dim_x; u++) {
          /* go through all u */
          if (index[j] < u + 10)
            continue;

          average_access[u][i] += master_access[j][u][index[j]];  /* index[j] position in sequence j */
        }
        index[j]++;                                               /* increase position in sequence */
      } else {
        count_j--;
      }
    }
    if (!(count_j > 0)) {
      printf("Alignment contains only gap at column %d\n exiting program\n", i);
      return NULL;
    }

    for (u = 1; u < dim_x; u++)
      average_access[u][i] = (int)rint(average_access[u][i] / (count_j));
  }
  /* free(index); */
  for (i = 0; i < number; i++) {
    int j;
    j = master_access[i][0][0];
    while (--j > -1)
      free(master_access[i][j]);
    free(master_access[i]);
  }
  free(master_access);
  free(index);
  return average_access;
}


/**
 * aliprint_struct generate two postscript files.
 * The first one is a structure annotated alignment a la RNAalifold.
 * The second file is a structural representation of the interaction
 * with colour annotation of the base-pair conservation.
 */
static void
aliprint_struct(FILE      *Result,        /* result file */
                FILE      *mrna,          /* target alignment */
                FILE      *query,         /* query alignment */
                const int WindowsLength   /* Number of nucleotide around the target sequence */
                )
{
  char  *result = NULL;
  /**
   * Defines arrays for names and sequences of the query and target alignment
   */
  char  *AS1[MAX_NUM_NAMES], *AS2[MAX_NUM_NAMES], *names1[MAX_NUM_NAMES], *names2[MAX_NUM_NAMES];
  int   n_seq, n_seq2, i, s;

  /**
   * Check if target and sequence alignment contains the same number of sequences
   * The user should ensure that the sequence are in the same species order in both the target and query alignment files.
   */
  n_seq   = read_clustal(mrna, AS1, names1);
  n_seq2  = read_clustal(query, AS2, names2);
  if (n_seq != n_seq2) {
    for (i = 0; AS1[i]; i++) {
      free(AS1[i]);
      free(AS2[i]);
    }
    vrna_message_error("unequal number of seqs in alignments");
  }

  /**
   * Here we get the length of the target and query alignments. Important for initiliazing the necessary arrays
   */
  int target_alignment_length, query_alignment_length;
  target_alignment_length = strlen(AS1[0]);
  query_alignment_length  = strlen(AS2[0]);


  /**
   * Initialize files containing sequences
   */
  AS1[n_seq]  = NULL;
  AS2[n_seq]  = NULL;
  do {
    /**
     * Quit if the result file does not exist
     */
    if ((result = vrna_read_line(Result)) == NULL) {
      free(result);
      break;
    }

    /**
     * If we have a line in the result file that looks like a result line, then parse it
     */
    if (strchr(result, '(') && strchr(result, '&') && strchr(result, '(') && strchr(result, ',') && strchr(result, ':') && strchr(result, '-')) {
      char  *structure, *pos;
      /**
       * tbegin,tend,qbegin,qend store the duplex boundaries
       * on both the target and query sequences, respectively
       */
      int   tbegin, tend, qbegin, qend;
      int   length, posi;
      /**
       * sbegin, send are the coordinates of the target alignment slide
       */
      int   sbegin, send;
      /**
       * Check in the line where is the first space.
       * This gives us the total length of the duplex
       */
      pos     = strchr(result, ' ');
      posi    = (int)(pos - result);
      length  = posi;
      /**
       * Copy the structure in the line into variable structure
       */
      structure = (char *)vrna_alloc((length + 1) * sizeof(char));
      sscanf(result, "%s", structure);
      /**
       * Save the coordinates on the target
       */
      sscanf(pos, "%10d,%10d", &tbegin, &tend);
      pos = strchr(pos, ':');
      pos = strchr(pos, ' ');
      /**
       * Save the coordinates on the query
       */
      sscanf(pos, "%10d,%10d", &qbegin, &qend);

      /**
       * To produce the ps files we use the full query alignment
       * and a user defined section of the target coordinates.
       * The size of the target slide is determined by the WindowsLength variable
       */
      sbegin  = MAX2(tbegin - WindowsLength, 0);
      send    = MIN2(tend + WindowsLength, target_alignment_length);
      /**
       * We retrieved the coordinates of the duplex
       * Now we can slide the alignment
       * Target and query will hold the alignments for the target and query sequences
       */
      char  **target;
      char  **query;
      char  **join;
      target  = (char **)vrna_alloc((n_seq + 1) * sizeof(char *));
      query   = (char **)vrna_alloc((n_seq + 1) * sizeof(char *));
      join    = (char **)vrna_alloc((n_seq + 1) * sizeof(char *));
      for (s = 0; s < n_seq; s++) {
        target[s] = (char *)vrna_alloc((send - sbegin + 2) * sizeof(char));
        strncpy(target[s], (AS1[s] + sbegin - 1), (send - sbegin + 1));
        target[s][send - sbegin + 1]  = '\0';
        query[s]                      = (char *)vrna_alloc((query_alignment_length + 1) * sizeof(char));
        strcpy(query[s], AS2[s]);
        query[s][query_alignment_length]  = '\0';
        join[s]                           = (char *)vrna_alloc((query_alignment_length + (send - sbegin) + 3) * sizeof(char));
        strncpy(join[s], (AS1[s] + sbegin - 1), (send - sbegin + 1));
        strcat(join[s], "&");
        strcat(join[s], AS2[s]);
      }
      target[n_seq] = NULL;
      query[n_seq]  = NULL;
      /**
       * Get the consensus structure for the query and target sequence
       * This consensus structure is constrained such that the interaction sites are unbound.
       * We define and set the query and target constraint
       */
      char  *target_constraint;
      char  *query_constraint;
      char  *joint_structure;
      target_constraint = (char *)vrna_alloc((unsigned)(send - sbegin + 2));
      query_constraint  = (char *)vrna_alloc((unsigned)(query_alignment_length + 1));
      joint_structure   = (char *)vrna_alloc((unsigned)(send - sbegin + 3 + query_alignment_length));
      for (i = 0; i < strlen(target[0]); i++) {
        if (i > tbegin - sbegin - 1 && i < tend - sbegin + 1)
          target_constraint[i] = 'x';
        else
          target_constraint[i] = '.';
      }
      for (i = 0; i < strlen(AS2[0]); i++) {
        if (i > qbegin - 2 && i < qend)
          query_constraint[i] = 'x';
        else
          query_constraint[i] = '.';
      }
      /**
       * Now produce structure
       */
      fold_constrained = 1;
      alifold((const char **)target, target_constraint);
      for (i = 0; !(target_constraint[i] == '\0'); i++) {
        if (target_constraint[i] == '(')
          target_constraint[i] = '[';
        else if (target_constraint[i] == ')')
          target_constraint[i] = ']';
      }
      alifold((const char **)query, query_constraint);
      for (i = 0; !(query_constraint[i] == '\0'); i++) {
        if (query_constraint[i] == '(')
          query_constraint[i] = '<';
        else if (query_constraint[i] == ')')
          query_constraint[i] = '>';
      }
      /**
       * Add duplex structure, and produce joint structure
       */
      int count_duplex_structure = 0;
      for (i = tbegin - sbegin; i < tend - sbegin + 1; i++) {
        target_constraint[i] = structure[count_duplex_structure];
        count_duplex_structure++;
      }
      count_duplex_structure++;
      for (i = qbegin - 1; i < qend; i++) {
        query_constraint[i] = structure[count_duplex_structure];
        count_duplex_structure++;
      }
      strncpy(joint_structure, target_constraint, (send - sbegin + 1));
      strcat(joint_structure, "&");
      strcat(joint_structure, query_constraint);
      /**
       * Produce consensus sequence
       */
      char  *temp_target;
      char  *temp_query;
      char  *string;
      temp_target = consensus((const char **)target);
      temp_query  = consensus((const char **)query);
      string      = (char *)vrna_alloc((strlen(temp_target) + strlen(temp_query) + 1) * sizeof(char));
      strcpy(string, temp_target);
      free(temp_target);
      strcat(string, temp_query);
      free(temp_query);
      /**
       * now produce output name, based on the first two names of the alignment
       */
      char *psoutput = vrna_strdup_printf("annaln_%s_%s_%d_%d_%d_%d.ps",
                                          names1[0],
                                          names2[0],
                                          tbegin,
                                          tend,
                                          qbegin,
                                          qend);
      s = 0;
      while (psoutput[s] != '\0') {
        if (psoutput[s] == '\\' || psoutput[s] == '/')
          psoutput[s] = '_';

        s++;
      }
      /**
       * Produce a structure annotated, colorated alignment
       */
      aliPS_color_aln(joint_structure, psoutput, (const char **)join, (const char **)names1);
      free(psoutput);
#if 0
      /**
       * We also need the consensus structure annotated with conservation information
       *
       */

      /**
       * First we change the output file name from annaln -> annstr
       */
      psoutput[3] = 's';
      psoutput[4] = 't';
      psoutput[5] = 'r';


      /**
       * Now we change the different parenthesis into one
       */

      char *joint_structure_parenthesis;
      joint_structure_parenthesis = (char *)vrna_alloc((strlen(joint_structure) + 1) * sizeof(char));
      strcpy(joint_structure_parenthesis, joint_structure);
      for (i = 0; i < strlen(joint_structure); i++) {
        if ((joint_structure[i] == '[') || (joint_structure[i] == '<'))
          joint_structure[i] = '(';
        else if (joint_structure[i] == ']' || joint_structure[i] == '>')
          joint_structure[i] = ')';
      }
      cut_point = send - sbegin + 1;
      /**
       * Remove & from the alignment and the consensus structure
       */
      int counter = 0;
      for (i = 0; i < strlen(joint_structure); i++) {
        if (joint_structure[i] == '&')
          continue;

        joint_structure[counter]              = joint_structure[i];
        joint_structure_parenthesis[counter]  = joint_structure_parenthesis[i];
        for (s = 0; s < n_seq; s++)
          join[s][counter] = join[s][i];
        counter++;
      }
      free(string);
      string = consensus((const char **)join);
      printf("%s\n%s\n%s\n", join[0], string, joint_structure);
      char **A;
      //A=annote(joint_structure_parenthesis,(const char**) join);
      xrna_plex_plot(string, joint_structure_parenthesis, psoutput);
      /**
       * Free stuff...
       */
      free(structure);
      free(result);
      free(psoutput);
      for (i = 0; i < n_seq + 1; i++)
        free(target[i]);
      free(target);
#endif
    }
  } while (1);
  for (i = 0; AS1[i]; i++) {
    free(AS1[i]);
    free(AS2[i]);
    free(names1[i]);
    free(names2[i]);
  }
  fclose(mrna);
  fclose(query);
}


/*
 * All characters not being included into ATCGUN are considered as gap characters.
 */
static int
get_sequence_length_from_alignment(char *sequence)
{
  int counter = 0;

  while (!(*sequence == '\0')) {
    if (*sequence == 'A' || *sequence == 'T' || *sequence == 'C' || *sequence == 'G' || *sequence == 'U' || *sequence == 'N')
      counter++;

    sequence++;
  }
  return counter;
}


static void
linear_fit(int  *il_a,
           int  *il_b,
           int  *b_a,
           int  *b_b)
{
  /*get fit parameters*/
  vrna_md_t     md;
  vrna_param_t  *P = NULL;

  set_model_details(&md);

  if ((!P) || (fabs(P->temperature - temperature) > 1e-6)) {
    update_fold_params();
    P = vrna_params(&md);
    make_pair_matrix();
  }

  int     internal_loop_x[25] = {
    6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
  };
  int     internal_loop_y[25];
  int     i;
  for (i = 0; i < 25; i++) /* initialize the y-value for internal loops */
    internal_loop_y[i] = P->internal_loop[internal_loop_x[i]];
  int     sumx, sumy, sumxx, sumxy;
  sumx = sumy = sumxx = sumxy = 0;
  double  til_a, til_b;
  int     n = 25; /* only take data points till int_loop size <=20;  */
                  /* no measurement exists for larger loop &&  */
                  /* such large loop are not expected in RNA-RNA interactions. */
  for (i = 0; i < n; i++) {
    sumx  += internal_loop_x[i];
    sumy  += internal_loop_y[i];
    sumxx += internal_loop_x[i] * internal_loop_x[i];
    sumxy += internal_loop_x[i] * internal_loop_y[i];
  }
  if ((sumxx - (sumx * sumx) / n) < 1e-6) {
    free(P);
    printf("divisor for internal loop is too small %d\n", (sumxx - (sumx * sumx) / n));
    vrna_message_error("Problem in fitting");
  }

  til_a = (sumxy - (sumx * sumy) / n) / (sumxx - (sumx * sumx) / n);
  til_b = sumy / n - (til_a * sumx / n);
  *il_a = (int)til_a;
  *il_b = (int)til_b;

  /* ###bulge#### */

  n = 5;
  int     bulge_loop_x[5] = {
    2, 3, 4, 5, 6
  };
  int     bulge_loop_y[5];
  for (i = 0; i < n; i++)
    bulge_loop_y[i] = P->bulge[bulge_loop_x[i]];
  sumx = sumy = sumxx = sumxy = 0;
  double  tb_a, tb_b;
  for (i = 0; i < n; i++) {
    sumx  += bulge_loop_x[i];
    sumy  += bulge_loop_y[i];
    sumxx += bulge_loop_x[i] * bulge_loop_x[i];
    sumxy += bulge_loop_x[i] * bulge_loop_y[i];
  }
  tb_a  = (sumxy - (sumx * sumy) / n) / (sumxx - (sumx * sumx) / n);
  tb_b  = sumy / n - (tb_a * sumx / n);
  *b_a  = (int)tb_a;
  *b_b  = (int)tb_b;
  free(P);
  if ((sumxx - (sumx * sumx) / n) < 1e-6) {
    printf("divisor for bulge loop is too small %d\n", (sumxx - (sumx * sumx) / n));
    vrna_message_error("Problem in fitting");
  }
}


double
probcompute_silvana_98(char   *s1,
                       double k_concentration,
                       double tris_concentration,
                       double mg_concentration,
                       double na_concentration,
                       double probe_concentration)
{
  double  d_a               = 3.92 / 100000.0;  /*Parameters from the article of Owczarzy*/
  double  d_b               = 9.11 / 1000000;   /*Parameters from the article of Owczarzy*/
  double  d_c               = 6.26 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_d               = 1.42 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_e               = 4.82 / 10000;     /*Parameters from the article of Owczarzy*/
  double  d_f               = 5.25 / 10000;     /*Parameters from the article of Owczarzy*/
  double  d_g               = 8.31 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_magn_corr_value = 0;
  double  d_fgc             = 0;
  double  dH, dS;

  dS  = 0;
  dH  = 0;
  double  salt_correction;
  double  enthalpy[4][4] =
  { { -7.9, -8.4,  -7.8,  -7.2 },
    { -8.5, -8.0,  -10.6, -7.8 },
    { -8.2, -10.6, -8.0,  -8.4 },
    { -7.2, -8.2,  -8.5,  -7.9 } };
  double  entropy[4][4] = { {
                              -22.2, -22.4, -21.0, -20.4
                            },
                            { -22.7, -19.9, -27.2, -21.0 },
                            { -22.2, -27.2, -19.9, -22.4 },
                            { -21.3, -22.2, -22.7, -22.2 } };
  int     posst;
  posst = 0;
  int     converted = encode_char(toupper(s1[posst])) - 1;
  int     seqlen    = strlen(s1);
  double  Tm;

  /* terminal correction */
  if (s1[0] == 'G' || s1[0] == 'C') {
    dH  += 0.1;
    dS  += -2.8;
    d_fgc++;
  }

  if (s1[0] == 'A' || s1[0] == 'T' || s1[0] == 'U') {
    dH  += 2.3;
    dS  += 4.8;
  }

  if (s1[seqlen - 1] == 'G' || s1[seqlen - 1] == 'C') {
    dH  += 0.1;
    dS  += -2.8;
  }

  if (s1[seqlen - 1] == 'A' || s1[seqlen - 1] == 'T' || s1[seqlen - 1] == 'U') {
    dH  += 2.3;
    dS  += 4.8;
  }

  /* compute dH and dH based on sequence */
  for (posst = 1; posst < seqlen; posst++) {
    if (toupper(s1[posst]) == 'G' || toupper(s1[posst]) == 'C')
      d_fgc++;

    int nextconverted = encode_char(toupper(s1[posst])) - 1;
    dH        += enthalpy[converted][nextconverted];
    dS        += entropy[converted][nextconverted];
    converted = nextconverted;
  }
  d_fgc = d_fgc / ((double)(seqlen));
  if (mg_concentration == 0) {
    dS  += 0.368 * (seqlen - 1) * log(na_concentration);
    Tm  = (1000 * dH) / (dS + 1.987 * log(probe_concentration / 4)) - 273.15;
    return Tm;
  } else {
    double  single_charged  = k_concentration + tris_concentration / 2 + na_concentration;
    double  d_ratio_ions    = sqrt(mg_concentration / single_charged);
    if (single_charged == 0) {
      d_magn_corr_value =
        d_a -
        d_b * log(mg_concentration) +
        d_fgc * (d_c + d_d * log(mg_concentration)) +
        1 / (2 * ((double)seqlen - 1)) *
        (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
    } else {
      if (d_ratio_ions < 0.22) {
        d_magn_corr_value = (4.29 * d_fgc - 3.95) * 1 / 100000 * log(single_charged) + 9.40 * 1 / 1000000 * pow(log(single_charged), 2);
      } else {
        if (d_ratio_ions < 6) {
          d_a = 3.92 / 100000 * (0.843 - 0.352 * sqrt(single_charged) * log(single_charged));
          d_d = 1.42 / 100000 * (1.279 - 4.03 / 1000 * log(single_charged) - 8.03 / 1000 * pow(log(single_charged), 2));
          d_g = 8.31 / 100000 * (0.486 - 0.258 * log(single_charged) + 5.25 / 1000 * pow(log(single_charged), 3));

          d_magn_corr_value = d_a - d_b * log
                                (mg_concentration) + d_fgc *
                              (d_c + d_d * log(mg_concentration)) + 1 / (2 * ((double)seqlen - 1)) * (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
        } else {
          d_magn_corr_value = d_a - d_b * log(mg_concentration) + d_fgc * (d_c + d_d * log(mg_concentration)) + 1 / (2 * ((double)seqlen - 1)) * (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
        }
      }
    }

    double temp_Tm = dH * 1000 / (dS + 1.987 * log(probe_concentration / 4));
    Tm = 1 / (1 / temp_Tm + d_magn_corr_value) - 273.15;
    return Tm;
  }
}


double
probcompute_silvana_04(char   *s1,
                       double k_concentration,
                       double tris_concentration,
                       double mg_concentration,
                       double na_concentration,
                       double probe_concentration)
{
  double  d_a               = 3.92 / 100000.0;  /*Parameters from the article of Owczarzy*/
  double  d_b               = 9.11 / 1000000;   /*Parameters from the article of Owczarzy*/
  double  d_c               = 6.26 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_d               = 1.42 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_e               = 4.82 / 10000;     /*Parameters from the article of Owczarzy*/
  double  d_f               = 5.25 / 10000;     /*Parameters from the article of Owczarzy*/
  double  d_g               = 8.31 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_magn_corr_value = 0;
  double  d_fgc             = 0;
  double  dH, dS;

  dS  = 0;
  dH  = 0;
  double  salt_correction;
  double  enthalpy[4][4] =
  { { -7.6, -8.4, -7.8,  -7.2 },
    { -8.5, -8.0, -10.6, -7.8 },
    { -8.2, -9.8, -8.0,  -8.4 },
    { -7.2, -8.2, -8.5,  -7.6 } };
  double  entropy[4][4] = { {
                              -21.3, -22.4, -21.0, -20.4
                            },
                            { -22.7, -19.9, -27.2, -21.0 },
                            { -22.2, -24.4, -19.9, -22.4 },
                            { -21.3, -22.2, -22.7, -21.3 } };
  int     posst;
  posst = 0;
  int     converted = encode_char(toupper(s1[posst])) - 1;
  int     seqlen    = strlen(s1);
  double  Tm;

  /* terminal correction */
  if (s1[0] == 'G' || s1[0] == 'C') {
    dH  += 0.1;
    dS  += -2.85;
    d_fgc++;
  }

  if (s1[0] == 'A' || s1[0] == 'T' || s1[0] == 'U') {
    dH  += 2.4;
    dS  += 4.1;
  }

  if (s1[seqlen - 1] == 'G' || s1[seqlen - 1] == 'C') {
    dH  += 0.1;
    dS  += -2.85;
  }

  if (s1[seqlen - 1] == 'A' || s1[seqlen - 1] == 'T' || s1[seqlen - 1] == 'U') {
    dH  += 2.4;
    dS  += 4.1;
  }

  /* compute dH and dH based on sequence */
  for (posst = 1; posst < seqlen; posst++) {
    if (toupper(s1[posst]) == 'G' || toupper(s1[posst]) == 'C')
      d_fgc++;

    int nextconverted = encode_char(toupper(s1[posst])) - 1;
    dH        += enthalpy[converted][nextconverted];
    dS        += entropy[converted][nextconverted];
    converted = nextconverted;
  }
  d_fgc = d_fgc / ((double)(seqlen));
  if (mg_concentration == 0) {
    dS  += 0.368 * (seqlen - 1) * log(na_concentration);
    Tm  = (1000 * dH) / (dS + 1.987 * log(probe_concentration / 4)) - 273.15;
    return Tm;
  } else {
    double  single_charged  = k_concentration + tris_concentration / 2 + na_concentration;
    double  d_ratio_ions    = sqrt(mg_concentration / single_charged);
    if (single_charged == 0) {
      d_magn_corr_value =
        d_a -
        d_b * log(mg_concentration) +
        d_fgc * (d_c + d_d * log(mg_concentration)) +
        1 / (2 * ((double)seqlen - 1)) *
        (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
    } else {
      if (d_ratio_ions < 0.22) {
        d_magn_corr_value = (4.29 * d_fgc - 3.95) * 1 / 100000 * log(single_charged) + 9.40 * 1 / 1000000 * pow(log(single_charged), 2);
      } else {
        if (d_ratio_ions < 6) {
          d_a = 3.92 / 100000 * (0.843 - 0.352 * sqrt(single_charged) * log(single_charged));
          d_d = 1.42 / 100000 * (1.279 - 4.03 / 1000 * log(single_charged) - 8.03 / 1000 * pow(log(single_charged), 2));
          d_g = 8.31 / 100000 * (0.486 - 0.258 * log(single_charged) + 5.25 / 1000 * pow(log(single_charged), 3));

          d_magn_corr_value = d_a - d_b * log
                                (mg_concentration) + d_fgc *
                              (d_c + d_d * log(mg_concentration)) + 1 / (2 * ((double)seqlen - 1)) * (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
        } else {
          d_magn_corr_value = d_a - d_b * log(mg_concentration) + d_fgc * (d_c + d_d * log(mg_concentration)) + 1 / (2 * ((double)seqlen - 1)) * (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
        }
      }
    }

    double temp_Tm = dH * 1000 / (dS + 1.987 * log(probe_concentration / 4));
    Tm = 1 / (1 / temp_Tm + d_magn_corr_value) - 273.15;
    return Tm;
  }
}


double
probcompute_xia_98(char   *s1,
                   double na_concentration,
                   double probe_concentration)
{
  double  dH, dS;

  dS  = 0;
  dH  = 0;
  double  salt_correction;
  double  enthalpy[4][4] =
  { { -6.820,  -11.400, -10.480, -9.380  },
    { -10.440, -13.390, -10.640, -10.480 },
    { -12.440, -14.880, -13.390, -11.400 },
    { -7.690,  -12.440, -10.440, -6.820  } };
  double  entropy[4][4] =
  { { -19.0, -29.5, -27.1, -26.7 },
    { -26.9, -32.7, -26.7, -27.1 },
    { -32.5, -36.9, -32.7, -29.5 },
    { -20.5, -32.5, -26.9, -19.0 }, };
  int     posst;
  posst = 0;
  int     converted = encode_char(toupper(s1[posst])) - 1;
  int     seqlen    = strlen(s1);
  double  Tm;

  /* terminal correction */
  if (s1[0] == 'G' || s1[0] == 'C') {
    dH  += 3.61;
    dS  += -1.5;
  }

  if (s1[0] == 'A' || s1[0] == 'T' || s1[0] == 'U') {
    dH  += 7.73;
    dS  += 9.0;
  }

  if (s1[seqlen - 1] == 'G' || s1[seqlen - 1] == 'C') {
    dH  += 3.61;
    dS  += -1.5;
  }

  if (s1[seqlen - 1] == 'A' || s1[seqlen - 1] == 'T' || s1[seqlen - 1] == 'U') {
    dH  += 7.73;
    dS  += 9.0;
  }

  /* compute dH and dH based on sequence */
  for (posst = 1; posst < seqlen; posst++) {
    int nextconverted = encode_char(toupper(s1[posst])) - 1;
    dH        += enthalpy[converted][nextconverted];
    dS        += entropy[converted][nextconverted];
    converted = nextconverted;
  }
  dS  += 0.368 * (seqlen - 1) * log(na_concentration);
  Tm  = (1000 * dH) / (dS + 1.987 * log(probe_concentration / 4)) - 273.15;
  return Tm;
}


double
probcompute_sug_95(char   *s1,
                   double na_concentration,
                   double probe_concentration)
{
  double  dH, dS;

  dS  = 0;
  dH  = 0;
  double  salt_correction;
  double  enthalpy[4][4] =
  { { -11.500, -7.800,  -7.000,  -8.300 },
    { -10.400, -12.800, -16.300, -9.100 },
    { -8.600,  -8.000,  -9.300,  -5.900 },
    { -7.800,  -5.500,  -9.000,  -7.800 } };
  double  entropy[4][4] =
  { { -36.4, -21.6, -19.7, -23.9 },
    { -28.4, -31.9, -47.1, -23.5 },
    { -22.9, -17.1, -23.2, -12.3 },
    { -23.2, -13.5, -26.1, -21.9 } };
  int     posst;
  posst = 0;
  int     converted = encode_char(toupper(s1[posst])) - 1;
  int     seqlen    = strlen(s1);
  double  Tm;

  /*  salt entropy correction von Ahsen 1999 */
  /* terminal correction */
  if (s1[0] == 'G' || s1[0] == 'C') {
    dH  += 0.95;
    dS  += -1.95;
  }

  if (s1[0] == 'A' || s1[0] == 'T' || s1[0] == 'U') {
    dH  += 0.95;
    dS  += -1.95;
  }

  if (s1[seqlen - 1] == 'G' || s1[seqlen - 1] == 'C') {
    dH  += 0.95;
    dS  += -1.95;
  }

  if (s1[seqlen - 1] == 'A' || s1[seqlen - 1] == 'T' || s1[seqlen - 1] == 'U') {
    dH  += 0.95;
    dS  += -1.95;
  }

  /* compute dH and dH based on sequence */
  for (posst = 1; posst < seqlen; posst++) {
    int nextconverted = encode_char(toupper(s1[posst])) - 1;
    dH        += enthalpy[converted][nextconverted];
    dS        += entropy[converted][nextconverted];
    converted = nextconverted;
  }
  /*  salt entropy correction von Ahsen 1999 */
  dS  += 0.368 * (seqlen - 1) * log(na_concentration);
  Tm  = (1000 * dH) / (dS + 1.987 * log(probe_concentration / 4)) - 273.15;
  return Tm;
}


double
probcompute_newparameters(char    *s1,
                          double  k_concentration,
                          double  tris_concentration,
                          double  mg_concentration,
                          double  na_concentration,
                          double  probe_concentration)
{
  /* ////////////////////////////////////////// */
  /* Folding Init */
  /* //////////////////////////////////////////   */
  vrna_md_t     md;
  vrna_param_t  *P = NULL;

  set_model_details(&md);

  if ((!P) || (fabs(P->temperature - temperature) > 1e-6)) {
    update_fold_params();
    P = vrna_params(&md);
    make_pair_matrix();
  }

  /* ///////////////////////////////////////// */
  /* Salt Variable Declaration */
  /* ///////////////////////////////////////// */

  double  d_temp;             /* melting temperature */
  double  d_temp_na;          /*melting temperature in 1M na+ */
  double  d_salt_corr_value = 0.0;
  double  d_magn_corr_value = 0.0;

  double  d_fgc;              /* gc content */
  double  d_conc_monovalents  = k_concentration + na_concentration + tris_concentration / 2;
  double  d_ratio_ions        = sqrt(mg_concentration) / d_conc_monovalents;
  double  d_a                 = 3.92 / 100000.0;  /*Parameters from the article of Owczarzy*/
  double  d_b                 = 9.11 / 1000000;   /*Parameters from the article of Owczarzy*/
  double  d_c                 = 6.26 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_d                 = 1.42 / 100000;    /*Parameters from the article of Owczarzy*/
  double  d_e                 = 4.82 / 10000;     /*Parameters from the article of Owczarzy*/
  double  d_f                 = 5.25 / 10000;     /*Parameters from the article of Owczarzy*/
  double  d_g                 = 8.31 / 100000;    /*Parameters from the article of Owczarzy*/

  /* ////////////////////////////////////// */
  /* Sequence Variable Declaration */
  /* ////////////////////////////////////// */

  int seqlen = strlen(s1);
  int posst;
  posst = 0;
  int converted     = encode_char(toupper(s1[posst]));
  int revconverted  = abs(5 - converted);

  /* ////////////////////////////////////// */
  /* Energy Variable Declaration */
  /* ////////////////////////////////////// */

  double dT, dG, dS, dH;
  dS  = 0;
  dT  = 0;
  dH  = 0;

  /* ////////////////////////////////////// */
  /* Computation */
  /* ////////////////////////////////////// */


  /* GC-Content */
  d_fgc = 0;
  for (posst = 0; posst < seqlen; posst++)
    if (s1[posst] == 'G' || s1[posst] == 'C' || s1[posst] == 'g' || s1[posst] == 'c')
      d_fgc++;

  d_fgc = d_fgc / seqlen;

  /* dH dG dS */
  int type = pair[converted][revconverted];
  /* Add twice the duplex penalty */
  dG  = 2 * P->DuplexInit;
  dH  = 2 * DuplexInitdH;
  dS  = (dH - dG) / (37 + K0) * 10;
  if (type > 2) {
    dG  += P->TerminalAU;
    dH  += TerminalAUdH;
    dS  += (TerminalAUdH - P->TerminalAU) / (37 + K0) * 10;
  }

  for (posst = 1; posst < seqlen; posst++) {
    int nextconverted     = encode_char(toupper(s1[posst]));
    int nextrevconverted  = abs(5 - nextconverted);
    int nexttype          = pair[nextconverted][nextrevconverted];
    dG            += stack37[rtype[type]][nexttype];
    dH            += stackdH[rtype[type]][nexttype];
    dS            += (stackdH[rtype[type]][nexttype] - stack37[rtype[type]][nexttype]) / (37 + K0) * 10;
    converted     = nextconverted;
    revconverted  = nextrevconverted;
    type          = nexttype;
  }
  if (type > 2) {
    dG  += P->TerminalAU;
    dH  += TerminalAUdH;
    dS  += (TerminalAUdH - P->TerminalAU) / (37 + K0) * 10;
  }

  /* add initiation again because of symmetry. */


  if (mg_concentration == 0) {
    dS += 0.368 * (seqlen - 1) * log(na_concentration);
    double Tm;
    Tm = (10 * dH) / (dS + 1.987 * log(probe_concentration / 4)) - 273.15;
    return Tm;
  } else {
    double  single_charged  = k_concentration + tris_concentration / 2 + na_concentration;
    double  d_ratio_ions    = sqrt(mg_concentration / single_charged);
    if (single_charged == 0) {
      d_magn_corr_value =
        d_a -
        d_b * log(mg_concentration) +
        d_fgc * (d_c + d_d * log(mg_concentration)) +
        1 / (2 * ((double)seqlen - 1)) *
        (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
    } else {
      if (d_ratio_ions < 0.22) {
        d_magn_corr_value = (4.29 * d_fgc - 3.95) * 1 / 100000 * log(single_charged) + 9.40 * 1 / 1000000 * pow(log(single_charged), 2);
      } else {
        if (d_ratio_ions < 6) {
          d_a = 3.92 / 100000 * (0.843 - 0.352 * sqrt(single_charged) * log(single_charged));
          d_d = 1.42 / 100000 * (1.279 - 4.03 / 1000 * log(single_charged) - 8.03 / 1000 * pow(log(single_charged), 2));
          d_g = 8.31 / 100000 * (0.486 - 0.258 * log(single_charged) + 5.25 / 1000 * pow(log(single_charged), 3));

          d_magn_corr_value = d_a - d_b * log
                                (mg_concentration) + d_fgc *
                              (d_c + d_d * log(mg_concentration)) + 1 / (2 * ((double)seqlen - 1)) * (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
        } else {
          d_magn_corr_value = d_a - d_b * log(mg_concentration) + d_fgc * (d_c + d_d * log(mg_concentration)) + 1 / (2 * ((double)seqlen - 1)) * (-d_e + d_f * log(mg_concentration) + d_g * pow(log(mg_concentration), 2));
        }
      }
    }

    double  temp_Tm = dH * 10 / (dS + 1.987 * log(probe_concentration / 4));
    int     Tm      = 1 / (1 / temp_Tm + d_magn_corr_value) - 273.15;
    return Tm;
  }
}
