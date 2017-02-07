/*
 *                Ineractive Access to folding Routines
 *
 *                c Ivo L Hofacker
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "ViennaRNA/snofold.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/snoop.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "RNAsnoop_cmdl.h"
static void  aliprint_struc(snoopT      *dup,
                            const char  **s1,
                            const char  **s2,
                            char **,
                            char **,
                            int         count,
                            int);


static void  print_struc(snoopT     *dup,
                         const char *s1,
                         const char *s2,
                         char *,
                         char *,
                         int,
                         int);


/* static void  print_struc_L(snoopT const *dup, const char *s1, const char *s2, char*, char*); */


static void redraw_output(char  *fname,
                          char  *output,
                          int   plfold_up_flag,
                          char  *suffix,
                          char  *access,
                          char  *sname,
                          char  *tname);


static int **read_plfold_i(char       *fname,
                           const int  beg,
                           const int  end);                             /* read plfold profiles */


static int **read_rnaup(char      *fname,
                        const int beg,
                        const int end);


static int get_max_u(const char *s,
                     char       delim);


extern int cut_point;


#define PRIVATE static
#define MAX_NUM_NAMES    500

/*--------------------------------------------------------------------------*/

int
main(int  argc,
     char *argv[])
{
  struct        RNAsnoop_args_info  args_info;
  char                              *string_s, *line_s, *name_s, *temp_s; /*string for snoRNA*/
  char                              *string_t, *line_t, *temp_t;          /*string for target RNA*/
  char                              *structure  = NULL, *cstruc = NULL;
  char                              *sname      = NULL, *tname = NULL /*name of the sno RNa file, mRNA file respectively*/;
  char                              *name       = NULL;
  char                              *access;
  char                              *result_file;
  char                              *output_directory;
  int                               fullStemEnergy;

  output_directory  = NULL;
  result_file       = NULL;
  access            = NULL;
  char                              *suffix;
  suffix = NULL;
  FILE                              *sno, *mrna;
  int                               fast;
  fast = 0;
  /* long int elapTicks; */
  /* clock_t Begin, End; */
  int                               plfold_up_flag = 0;
  int                               nice, i, r, l, length_s,  /* , length_t,  */
                                    penalty,                  /*extension penalty*/
                                    threshloop,               /*energy threshold on loop*/
                                    threshLE,                 /*energy threshold on the S2 part*/
                                    threshRE,                 /*energy threshold on the S1 part*/
                                    threshDE,                 /*Total duplex energy threshold*/
                                    threshTE,                 /*Total duplex + Loop energy*/
                                    threshSE,                 /*Sum of all energies*/
                                    threshD,                  /*lower stem energy*/
                                    half_stem,                /*minimal length of stem*/
                                    max_half_stem,            /*maximal length of stem*/
                                    max_s2 /*maximale position wo stem anfangen darf in 3->5*/,
                                    min_s2 /*minimale position wo stem anfangen darf in 3->5*/,
                                    max_s1 /*minimal position wo stem enden darf in 3->5 */,
                                    min_s1 /*maximale position wo stem enden darf in 3->5*/,
                                    distance /*distance between two subopts*/,
                                    min_d1 /*minimal distance between 5' sno and first duplex*/,
                                    min_d2 /*minimal distance between 3' sno and second duplex*/,
                                    max_asymm /*maximal asymmetry in the stem interior loop*/,
                                    alignment_length /*maximal target RNA alignment length*/,
                                    alignment /*flag to use ali version*/,
                                    delta,
                                    redraw /*if used (option I) allow to redraw command line output into ps files */;

  int noconv = 0;
  string_s            = NULL;
  string_t            = NULL;
  plfold_up_flag      = 0;
  alignment           = 0;
  redraw              = 0;
  nice                = 0;
  fast                = 0;
  snoop_subopt_sorted = 0;
  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAsnoop_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* Redraw results from RNAsnoop */
  if (args_info.constraint_given)
    fold_constrained = 1;

  /* structure constraint */
  if (args_info.produce_ps_given)
    redraw = 1;

  /*output directory for RNAsnoop results */
  output_directory = strdup(args_info.output_directory_arg);
  /* Draw directly nice pictures from RNAsnoop*/
  if (args_info.direct_redraw_given)
    nice = 1;

  /*Accessibility files are coming from RNAup*/
  if (args_info.from_RNAup_given) {
    plfold_up_flag  = 2;
    access          = strdup(args_info.from_RNAup_arg);
  }

  /*Suffix for the name of the accessibility files*/
  suffix = strdup(args_info.suffix_arg);
  /*Accessibility files are coming from RNAup*/
  if (args_info.from_RNAplfold_given) {
    plfold_up_flag  = 1;
    access          = strdup(args_info.from_RNAplfold_arg);
  }

  /*Check if we are in alignment mode*/
  if (args_info.alignment_mode_given)
    alignment = 1;

  /*Name of the stems file*/
  if (args_info.query_given)
    sname = strdup(args_info.query_arg);

  /*Name of the target file*/
  if (args_info.target_given)
    tname = strdup(args_info.target_arg);

  /*Delta for suboptimals*/
  delta = (int)(100 * (args_info.energy_threshold_arg + 0.1));
  /*fast folding and target search*/
  fast = args_info.fast_folding_arg;
  /*Penalty for duplex extension*/
  penalty = args_info.extension_cost_arg;
  /*threshold on loop energy*/
  threshloop = args_info.minimal_loop_energy_arg;
  /*threshold on right duplex*/
  threshRE = args_info.minimal_right_duplex_arg;
  /*threshold on left duplex*/
  threshLE = args_info.minimal_left_duplex_arg;
  /*threshold on minimal duplex*/
  threshDE = args_info.minimal_duplex_arg;
  /*threshold on duplex distance*/
  distance = args_info.duplex_distance_arg;
  /*threshold on minimal stem length*/
  half_stem = args_info.minimal_stem_length_arg;
  /*threshold on maximal stem length*/
  max_half_stem = args_info.maximal_stem_length_arg;
  /*threshold on minimal duplex box length*/
  min_s2 = args_info.minimal_duplex_box_length_arg;
  /*threshold on maximal duplex box length*/
  max_s2 = args_info.maximal_duplex_box_length_arg;
  /*threshold on minimal snoRNA stem loop length*/
  min_s1 = args_info.minimal_snoRNA_stem_loop_length_arg;
  /*thrshold on maximal snoRNA stem loop length*/
  max_s1 = args_info.maximal_snoRNA_stem_loop_length_arg;
  /*threshold on minimal snoRNA duplex length*/
  min_d1 = args_info.minimal_snoRNA_duplex_length_arg;
  /*threshold on minimal snoRNA duplex length*/
  min_d2 = args_info.minimal_snoRNA_duplex_length_arg;
  /*threshold on minimal duplex stem energy*/
  threshTE = args_info.minimal_duplex_stem_energy_arg;
  /*threshold on minimal total energy*/
  threshSE = args_info.minimal_total_energy_arg;
  /*threshold on maximal stem asymmetry*/
  max_asymm = args_info.maximal_stem_asymmetry_arg;
  /*threshold on minimal lower stem energy*/
  threshD = args_info.minimal_lower_stem_energy_arg;
  /*threshold on minimal lower stem energy*/
  alignment_length = args_info.alignmentLength_arg;

  threshloop = MIN2(threshloop, 0);

  /*   if(plfold_up_flag && !fast){ */
  /*     printf("Sorry, no accessibility implementation with the non fast implementation\n"); */
  /*     printf("If you want to include accessibility information please run RNAsnoop with the -f option\n"); */
  /*     printf("If you want to run RNAsnoop with the slow algorithm, please remove run RNAsnoop without -P\n"); */
  /*     exit(0); */
  /*   } */


  if (plfold_up_flag == 2 && suffix == NULL) {
    printf("RNAsnoop needs a suffix (-S u1-to-30) for the RNAup accessibility file in order to localize them\n");
    exit(0);
  }

  if (redraw) {
    if (!alignment)
      redraw_output(result_file, output_directory, plfold_up_flag, suffix, access, NULL, NULL);
    else
      redraw_output(result_file, output_directory, plfold_up_flag, suffix, access, sname, tname);

    /* readfile */
    /* parse lines */
    /* use current ploting function to output the results */
    exit(0);
  }

  /* istty = isatty(fileno(stdout))&&isatty(fileno(stdin)); */

  min_s2  += 5;
  max_s2  += 5;
  min_s1  += 5;
  max_s1  += 5;
  min_d1  += 5;
  min_d2  += 5;


  if (!alignment) {
    if (tname == NULL || sname == NULL)
      RNAsnoop_cmdline_parser_print_help();

    sno = fopen(sname, "r");
    if (sno == NULL) {
      printf("%s: Wrong snoRNA file name\n", sname);
      return 0;
    }

    mrna = fopen(tname, "r");
    if (mrna == NULL) {
      printf("%s: Wrong target file name\n", tname);
      return 0;
    }

    do {
      /* main loop: continue until end of file */
      if ((line_s = vrna_read_line(sno)) == NULL) {
        free(line_s);
        break;
      }

      /* skip comment lines and get filenames */
      while ((*line_s == '*') || (*line_s == '\0') || (*line_s == '>')) {
        if (*line_s == '>')
          name_s = (char *)vrna_alloc(strlen(line_s) + 1);

        (void)sscanf(line_s, "%s", name_s);
        free(line_s);
        if ((line_s = vrna_read_line(sno)) == NULL) {
          free(line_s);
          break;
        }
      }
      if (name_s == NULL) {
        printf("Your snoRNA sequence: \n%s\nhas no header. Please update your fasta file\n", line_s);
        exit(0);
      }

      /*   if ((line ==NULL) || (strcmp(line, "@") == 0)) break; */
      temp_s = (char *)vrna_alloc(strlen(line_s) + 1);
      (void)sscanf(line_s, "%s", temp_s);
      free(line_s);
      length_s = (int)strlen(temp_s);
      for (l = 0; l < length_s; l++) {
        temp_s[l] = toupper(temp_s[l]);
        if (!noconv && temp_s[l] == 'T')
          temp_s[l] = 'U';
      }
      string_s = (char *)vrna_alloc(length_s + 11);
      strcpy(string_s, "NNNNN");
      strcat(string_s + 5, temp_s);
      strcat(string_s + 5 + length_s, "NNNNN\0");
      free(temp_s);
      /* We declare the structure variable here as it will also contains the final stem structure */
      structure = (char *)vrna_alloc((unsigned)length_s + 11);
      if (fold_constrained) {
        cstruc = vrna_read_line(sno);
        if (cstruc != NULL) {
          int dn3 = strlen(cstruc) - (length_s - 10);
          strcpy(structure, ".....");
          strcat(structure, cstruc);
          if (dn3 >= 0) {
            strcat(structure, ".....\0");
          } else {
            while (dn3++)
              strcat(structure, ".");
            strcat(structure, "\0");
          }

          /* Now we fold with constraints the  */
        }
      }

      fullStemEnergy = snofold(string_s, structure, max_asymm, threshloop, min_s2, max_s2, half_stem, max_half_stem);
      do {
        /* main loop for target continue until end of file */
        snoopT  mfe;
        if ((line_t = vrna_read_line(mrna)) == NULL)
          /* free(line_t); */
          break;

        char    *name_t = NULL;

        /* skip comment lines and get filenames */
        while ((*line_t == '*') || (*line_t == '\0') || (*line_t == '>')) {
          if (*line_t == '>') {
            printf("%s\n", name_s);
            name_t = (char *)vrna_alloc(strlen(line_t) + 1);
            (void)sscanf(line_t, "%s", name_t);

            printf("%s\n", name_t);
            free(line_t);
            /*             free(string_t); */
          }

          if ((line_t = vrna_read_line(mrna)) == NULL) {
            free(line_t);
            break;
          }
        }
        if (name_t == NULL) {
          printf("Your target sequence: \n%s\nhas no header. Please update your fasta file\n", line_t);
          exit(0);
        }

        /*   if ((line ==NULL) || (strcmp(line, "@") == 0)) break; */
        temp_t = (char *)vrna_alloc(strlen(line_t) + 1);
        (void)sscanf(line_t, "%s", temp_t);
        int length_t;
        length_t = (int)strlen(temp_t);
        for (l = 0; l < length_t; l++) {
          temp_t[l] = toupper(temp_t[l]);
          if (!noconv && temp_t[l] == 'T')
            temp_t[l] = 'U';
        }
        string_t = (char *)vrna_alloc(length_t + 11);
        strcpy(string_t, "NNNNN");
        strcat(string_t + 5, temp_t);
        strcat(string_t + 5 + length_t, "NNNNN");
        free(temp_t);
        char *name_output;
        name_output = NULL;
        if (nice) {
          name_output = (char *)vrna_alloc(sizeof(char) * (length_t + length_s + 2));
          strcpy(name_output, name_t + 1);
          strcat(name_output, "_");
          strcat(name_output, name_s + 1);
          name_output[length_t + length_s + 1] = '\0';
        }

        if (delta >= 0) {
          snoopT  *subopt;
          snoopT  *sub;
          if (!fast && !plfold_up_flag) {
            subopt = snoop_subopt(string_t, string_s, delta, 5, penalty, threshloop,
                                  threshLE, threshRE, threshDE, threshTE, threshSE, threshD, distance,
                                  half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, fullStemEnergy);
            if (subopt == NULL) {
              printf("no target found under the given constraints\n");
              free(subopt);
              free(string_t);
              free(line_t);
              continue;
            }

            int count = 0;
            for (sub = subopt; sub->structure != NULL; sub++) {
              print_struc(sub, string_t, string_s, name_s, name_t, count++, nice);
              free(sub->structure);
            }
            free(subopt);
            free(string_t);
            free(line_t);
          } else if (!plfold_up_flag) {
            Lsnoop_subopt_list(string_t, string_s, delta, 5, penalty, threshloop,
                               threshLE, threshRE, threshDE, threshTE, threshSE, threshD, distance,
                               half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, alignment_length, name_output, fullStemEnergy);
            free(string_t);
            free(line_t);
          } else {
            if (plfold_up_flag == 1) {
              int   **access_s1;
              char  *file_s1;
              int   s1_len;/* k;*/ /*,j; */
              s1_len  = strlen(string_t);
              file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(name_t + 1) + strlen(access) + 9));
              strcpy(file_s1, access);
              strcat(file_s1, "/");
              strcat(file_s1, name_t + 1);
              strcat(file_s1, "_openen");
              access_s1 = read_plfold_i(file_s1, 1, s1_len);
              if (fast) {
                Lsnoop_subopt_list_XS(string_t, string_s, (const int **)access_s1, delta, 5, penalty, threshloop,
                                      threshLE, threshRE, threshDE, threshTE, threshSE, threshD, distance,
                                      half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, alignment_length, name_output, fullStemEnergy);
              } else {
                snoop_subopt_XS(string_t, string_s, (const int **)access_s1, delta, 5, penalty, threshloop,
                                threshLE, threshRE, threshDE, threshTE, threshSE, threshD, distance,
                                half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, alignment_length, name_output, fullStemEnergy);
              }

              /*                for(j=0; j<s1_len; j++){  */
              /*                  for(k=0; k<access_s1[0][0];k++){  */
              /*                    printf("%d \t",access_s1[k][j]);  */
              /*                  }  */
              /*                  printf("\n");  */
              /*                }  */
              free(file_s1);
              free(string_t);
              free(line_t);
              int i = access_s1[0][0];
              while (--i > -1)
                free(access_s1[i]);
              free(access_s1);
            } else if (plfold_up_flag == 2) {
              int   **access_s1;
              char  *file_s1;
              int   s1_len, k;/* ,j; */
              s1_len  = strlen(string_t);
              file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(name_t + 1) + strlen(suffix) + strlen(access) + 3));
              strcpy(file_s1, access);
              strcat(file_s1, "/");
              strcat(file_s1, name_t + 1);
              strcat(file_s1, "_");
              strcat(file_s1, suffix);
              access_s1 = read_rnaup(file_s1, 1, s1_len);
              if (fast) {
                Lsnoop_subopt_list_XS(string_t, string_s, (const int **)access_s1, delta, 5, penalty, threshloop,
                                      threshLE, threshRE, threshDE, threshTE, threshSE, threshD, distance,
                                      half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, alignment_length, name_output, fullStemEnergy);
              } else {
                snoop_subopt_XS(string_t, string_s, (const int **)access_s1, delta, 5, penalty, threshloop,
                                threshLE, threshRE, threshDE, threshTE, threshSE, threshD, distance,
                                half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, alignment_length, name_output, fullStemEnergy);
              }

              /*                for(j=0; j<s1_len; j++){  */
              /*                  for(k=0; k<access_s1[0][0];k++){  */
              /*                    printf("%d \t",access_s1[k][j]);  */
              /*                  }  */
              /*                  printf("\n");  */
              /*                }  */
              free(file_s1);
              free(string_t);
              free(line_t);
              k = access_s1[0][0];
              while (--k > -1)
                free(access_s1[k]);
              free(access_s1);
            }
          }
        } else {
          mfe = snoopfold(string_t, string_s, penalty, threshloop, threshLE, threshRE, threshDE, threshD,
                          half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, fullStemEnergy);
          if (mfe.energy < INF) {
            print_struc(&mfe, string_t, string_s, name_s, name_t, 0, 1);
            free(mfe.structure);
          }

          free(line_t);
          length_t = (int)strlen(string_t);
          free(string_t);
        }

        free(name_t);
        if (nice)
          free(name_output);
      } while (1);
      rewind(mrna);
      snofree_arrays(strlen(string_s));  /* free's base_pair */
      free(string_s);
      string_s = NULL;
      free(name_s);
      name_s = NULL;
    } while (1);
  } else {
    if (tname == NULL || sname == NULL)
      RNAsnoop_cmdline_parser_print_help();

    sno = fopen(sname, "r");
    if (sno == NULL) {
      printf("%s: Wrong snoRNA file name\n", sname);
      return 0;
    }

    mrna = fopen(tname, "r");
    if (mrna == NULL) {
      printf("%s: Wrong target file name\n", tname);
      return 0;
    }

    char  *temp1[MAX_NUM_NAMES], *temp2[MAX_NUM_NAMES], *AS1[MAX_NUM_NAMES], *AS2[MAX_NUM_NAMES], *names1[MAX_NUM_NAMES], *names2[MAX_NUM_NAMES];
    int   n_seq, n_seq2;
    n_seq   = read_clustal(mrna, temp1, names1);
    n_seq2  = read_clustal(sno, temp2, names2);
    if (n_seq != n_seq2) {
      for (i = 0; temp1[i]; i++) {
        free(temp1[i]);
        free(temp2[i]);
      }
      vrna_message_error("unequal number of seqs in alignments");
    }

    for (i = 0; temp1[i]; i++) {
      AS1[i]  = (char *)vrna_alloc((strlen(temp1[i]) + 11) * sizeof(char));
      AS2[i]  = (char *)vrna_alloc((strlen(temp2[i]) + 11) * sizeof(char));
      strcpy(AS1[i], "NNNNN");
      strcpy(AS2[i], "NNNNN");
      strcat(AS1[i], temp1[i]);
      strcat(AS2[i], temp2[i]);
      strcat(AS1[i], "NNNNN\0");
      strcat(AS2[i], "NNNNN\0");
    }
    for (i = 0; temp1[i]; i++) {
      free(temp1[i]);
      free(temp2[i]);
    }
    AS1[n_seq]  = NULL;
    AS2[n_seq]  = NULL;
    update_fold_params();
    alisnofold((const char **)AS2, max_asymm, threshloop, min_s2, max_s2, half_stem, max_half_stem);
    snoopT struc;
    struc = alisnoopfold((const char **)AS1, (const char **)AS2,
                         penalty, threshloop,
                         threshLE, threshRE, threshDE, threshD,
                         half_stem, max_half_stem,
                         min_s2, max_s2, min_s1,
                         max_s1, min_d1, min_d2);
    if (!(struc.structure == NULL)) {
      aliprint_struc(&struc, (const char **)AS1, (const char **)AS2, names1, names2, 0, nice);
      free(struc.structure);
    }

    snoopT *subopt = alisnoop_subopt((const char **)AS1, (const char **)AS2, delta, 5, penalty, threshloop,
                                     threshLE, threshRE, threshDE, threshD,
                                     threshTE, threshSE, distance,
                                     half_stem, max_half_stem, min_s2, max_s2,
                                     min_s1, max_s1, min_d1, min_d2);


    snoopT *sub;
    if (subopt == NULL) {
      printf("no target found under the given constraints\n");
    } else {
      int count = 1;
      for (sub = subopt; !(sub->structure == NULL); sub++) {
        aliprint_struc(sub, (const char **)AS1, (const char **)AS2, names1, names2, count, nice);
        free(sub->structure);
        count++;
      }
    }

    alisnofree_arrays(strlen(AS2[0]));
    free(subopt);
    for (i = 0; temp1[i]; i++) {
      free(AS1[i]);
      free(AS2[i]);
      free(names1[i]);
      free(names2[i]);
    }
  }

  fclose(sno);
  fclose(mrna);
  return 0;
}


static void
print_struc(snoopT      *dup,
            const char  *s1,
            const char  *s2,
            char        *name_s,
            char        *name_t,
            int         count,
            int         nice)
{
  int   l1;

  l1 = strchr(dup->structure, '&') - dup->structure;
  char  *target_struct;
  int   shift = 0, n2;
  char  *psoutput;
  psoutput = (char *)vrna_alloc(100 * sizeof(char));
  /*   if(dup->i > strlen(s1)-10){ */
  /*         dup->i--; */
  /*         l1--; */
  /*   } */
  /*   if(dup->i-l1< 0){ */
  /*         l1--; */
  /*         shift++; */
  /*   } */
  target_struct = (char *)vrna_alloc(sizeof(char) * (strlen(dup->structure) + 1));
  strncpy(target_struct, dup->structure + shift, l1);
  strncat(target_struct, dup->structure + (strchr(dup->structure, '&') - dup->structure), strlen(dup->structure) - (strchr(dup->structure, '&') - dup->structure));
  strcat(target_struct, "\0");
  char  *target;
  target = (char *)vrna_alloc(l1 + 1);
  strncpy(target, (s1 + dup->i - l1 + 5), l1);
  target[l1] = '\0';
  char  *s4;
  n2  = strlen(s2);
  s4  = (char *)vrna_alloc(sizeof(char) * (n2 - 9));
  strncpy(s4, s2 + 5, n2 - 10);
  s4[n2 - 10] = '\0';
  printf("%s %3d,%-3d;%3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f + %5.2f + 4.1 ) (%5.2f) \n%s&%s\n",
         target_struct, dup->i + 1 - l1,
         dup->i, dup->u, dup->j + 1, dup->j + (int)(strrchr(dup->structure, '>') - strchr(dup->structure, '>')) + 1,
         (dup->Loop_D + dup->Duplex_El + dup->Duplex_Er + dup->Loop_E) + 4.10,
         dup->Duplex_El, dup->Duplex_Er, dup->Loop_E, dup->Loop_D, dup->fullStemEnergy, target, s4);
  if (nice) {
    char  *temp_seq;
    char  *temp_struc;
    temp_seq    = (char *)vrna_alloc(sizeof(char) * (l1 + n2 - 9));
    temp_struc  = (char *)vrna_alloc(sizeof(char) * (l1 + n2 - 9));
    strcpy(temp_seq, target);
    strcat(temp_seq, s4);
    strncpy(temp_struc, target_struct, l1);
    strcat(temp_struc, target_struct + l1 + 1);
    temp_seq[n2 + l1 - 10]    = '\0';
    temp_struc[n2 + l1 - 10]  = '\0';
    cut_point                 = l1 + 1;

    psoutput = vrna_strdup_printf("sno_%d_u_%d_%s_%s.ps",
                                  count,
                                  dup->u,
                                  name_t + 1,
                                  name_s + 1);

    PS_rna_plot_snoop_a(temp_seq, temp_struc, psoutput, NULL, NULL);
    cut_point = -1;
    free(temp_seq);
    free(temp_struc);
    free(psoutput);
  }

  free(s4);
  free(target_struct);
  free(target);
}


static void
aliprint_struc(snoopT     *dup,
               const char **s1,
               const char **s2,
               char       **name_t,
               char       **name_s,
               int        count,
               int        nice)
{
  int   l1;

  l1 = strchr(dup->structure, '&') - dup->structure;
  char  *target_struct;
  int   shift = 0;
  /*   if(dup->i > strlen(s1[0])-10){ */
  /*     dup->i--; */
  /*     l1--; */
  /*   } */
  /*   if(dup->i-l1< 0){ */
  /*     l1--; */
  /*     shift++; */
  /*   } */
  int length_struct = strlen(dup->structure);
  target_struct = (char *)vrna_alloc(sizeof(char) * (length_struct + 1));
  strncpy(target_struct, dup->structure + shift, l1);
  strncat(target_struct, dup->structure + (strchr(dup->structure, '&') - dup->structure), length_struct - (strchr(dup->structure, '&') - dup->structure));
  strcat(target_struct, "\0");
  /* get the corresponding alignment slice */
  int   n_seq, s;
  for (s = 0; s1[s] != NULL; s++);
  n_seq = s;
  int   n1, n2;
  n1  = strlen(s1[0]);
  n2  = strlen(s2[0]);
  char  **target;
  target = (char **)vrna_alloc((n_seq + 1) * sizeof(char *));
  for (s = 0; s < n_seq; s++) {
    target[s] = (char *)vrna_alloc((l1 + n2 - 8) * sizeof(char));
    strncpy(target[s], (s1[s] + dup->i - l1 + 5), l1);
    strcat(target[s], "&");
    strncat(target[s], (s2[s] + 5), n2 - 10);
    target[s][l1 + n2 - 9] = '\0';
  }
  char *consens;
  consens     = consens_mis((const char **)target);
  consens[l1] = '&';

  printf("%s %3d,%-3d;%3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f + %5.2f + 4.1; duplex cov = %5.2f; stem cov = %5.2f )\n%s\n",
         dup->structure, dup->i + 1 - l1,
         dup->i, dup->u, dup->j + 1, dup->j + (int)(strrchr(dup->structure, '>') - strchr(dup->structure, '>')) + 1,
         (dup->Loop_D + dup->Duplex_El + dup->Duplex_Er + dup->Loop_E) / n_seq + 4.10,
         dup->Duplex_El / n_seq, dup->Duplex_Er / n_seq, dup->Loop_E / n_seq, dup->Loop_D / n_seq, dup->pscd / n_seq, dup->psct / n_seq, consens);
  if (nice) {
    char  *psoutput;
    psoutput = (char *)vrna_alloc(100 * sizeof(char));

    char  *temp_seq, *temp_struct, **temp_target;
    temp_seq    = (char *)vrna_alloc((l1 + n2 - 9) * sizeof(char));
    temp_struct = (char *)vrna_alloc((l1 + n2 - 9) * sizeof(char));
    temp_target = (char **)vrna_alloc((n_seq + 1) * sizeof(char *));
    for (s = 0; s < n_seq; s++) {
      temp_target[s] = (char *)vrna_alloc((l1 + n2 - 9) * sizeof(char));
      strncpy(temp_target[s], (s1[s] + dup->i - l1 + 5), l1);
      strncat(temp_target[s], (s2[s] + 5), n2 - 10);
      temp_target[s][l1 + n2 - 10] = '\0';
    }
    strncpy(temp_seq, consens, l1);
    strncpy(temp_struct, target_struct, l1);
    strcat(temp_seq, consens + l1 + 1);
    strcat(temp_struct, target_struct + l1 + 1);
    temp_seq[n2 - 10 + l1]    = '\0';
    temp_struct[n2 - 10 + l1] = '\0';
    char  str[16];
    char  upos[16];
    char  *temp;
    int   length_name = strlen(name_t[0]) + strlen(name_s[0]) + 1;
    temp = (char *)vrna_alloc(sizeof(char) * (length_name + 1));
    strcpy(temp, name_t[0]);
    strcat(temp, "_");
    strcat(temp, name_s[0]);
    temp[length_name] = '\0';
    strcpy(psoutput, "snoaln_");
    sprintf(str, "%d", count);
    strcat(psoutput, str);
    sprintf(upos, "%d", dup->u);
    strcat(psoutput, "_u_");
    strcat(psoutput, upos);
    strcat(psoutput, "_");
    for (s = 0; s < length_name; s++)
      if (temp[s] == '/')
        temp[s] = '-';

    strcat(psoutput, temp);
    strcat(psoutput, ".ps\0");
    /*   psoutput[strlen(temp)+4+strlen(str)+39]='\0'; */
    cut_point = l1 + 1;
    aliPS_color_aln(dup->structure, psoutput, (const char **)target, (const char **)name_t);
    psoutput[1] = 't';
    psoutput[2] = 'r';
    PS_rna_plot_snoop_a(temp_seq, temp_struct, psoutput, NULL, (const char **)temp_target);
    cut_point = -1;
    free(psoutput);
    for (s = 0; s < n_seq; s++)
      free(temp_target[s]);
    free(temp);
    free(temp_seq);
    free(temp_struct);
    free(temp_target);
  }

  for (s = 0; s < n_seq; s++)
    free(target[s]);
  free(consens);
  free(target_struct);
  free(target);
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
read_rnaup(char       *fname,
           const int  beg,
           const int  end)
{
  FILE  *in = fopen(fname, "r"); /*  open RNA_up file */

  int   i, j;
  int   **access;
  int   offset, temp;

  temp    = 0;
  offset  = 0;
  int   seq_pos;
  int   beg_r, end_r;
  beg_r = beg;
  end_r = end;

  if (in == NULL) {
    printf("%s", fname);
    perror("RNAup File open error here, Computing next target");

    exit(EXIT_FAILURE);
  }

  char tmp[2048] = {
    0x0
  };
  /* char *ptr; */


  if (fgets(tmp, sizeof(tmp), in) == 0)
    perror("Empty File");

  if (strchr(tmp, '>')) {
    vrna_message_error("file %s is not in RNAup format", fname);
    exit(EXIT_FAILURE);
  }

  while (!strstr(fgets(tmp, sizeof(tmp), in), "pos")) {
  }
  ;
  /*  (void) fgets(tmp,sizeof(tmp),in); //Datum  */
  /*   (void) fgets(tmp,sizeof(tmp),in); //white line */
  /*   (void) fgets(tmp,sizeof(tmp),in); //RNAup  */
  /*   (void) fgets(tmp,sizeof(tmp),in); //sequence length */
  /*   (void) fgets(tmp,sizeof(tmp),in); //sequence */

  int dim_x;
  dim_x   = get_max_u(tmp, 'S'); /*  get unpaired regions by conting tabs in second line */
  access  = (int **)vrna_alloc(sizeof(int *) * (dim_x + 2));
  for (i = 0; i < dim_x + 2; i++)
    access[i] = (int *)vrna_alloc(sizeof(int) * (end_r - beg_r + 7));
  for (i = 0; i < end_r - beg_r + 6; i++)
    for (j = 0; j < dim_x + 2; j++)
      access[j][i] = INF;
  access[0][0] = dim_x + 2;
  while (fgets(tmp, sizeof(tmp), in) != 0 && --end_r > -1) {
    /* read a record */
    float n;
    /* int i; */
    /* int u; */
    beg_r--;
    if (beg_r < 1) {
      if (sscanf(tmp, "%d%n", &seq_pos, &temp) == 1) {
        /*  read sequenz pos = 1. spalte */
        offset += temp;
        int count;
        count = 1;
        while (sscanf(tmp + offset, "%f%n", &n, &temp) == 1) {
          /* read columns */
          offset                        += temp;
          access[count][-beg_r + 5 + 1] = (int)rint(100 * n); /* seq_pos+5 */
          /*           printf("%d %d %f\n", count, -beg_r, access[count][-beg_r]); */
          count++;
        }
        offset = 0;
      }
    }
  }
  fclose(in);
  return access;
}


static int **
read_plfold_i(char      *fname,
              const int beg,
              const int end)
{
  FILE  *in = fopen(fname, "r");
  int   i, j;
  int   **access;
  int   offset, temp;

  temp    = 0;
  offset  = 0;
  int   seq_pos;
  int   beg_r, end_r;
  beg_r = beg;
  end_r = end;

  if (in == NULL) {
    perror(" open error");
    exit(EXIT_FAILURE);
  }

  char tmp[2048] = {
    0x0
  };
  /* char *ptr; */
  if (fgets(tmp, sizeof(tmp), in) == 0)
    perror("Empty File");

  if (strchr(tmp, '>')) {
    vrna_message_error("file %s is not in RNAplfold format", fname);
    exit(EXIT_FAILURE);
  }

  if (fgets(tmp, sizeof(tmp), in) == 0)
    perror("No accessibility data");

  int dim_x;
  dim_x   = get_max_u(tmp, '\t');
  access  = (int **)vrna_alloc(sizeof(int *) * (dim_x + 2));
  for (i = 0; i < dim_x + 2; i++)
    access[i] = (int *)vrna_alloc(sizeof(int) * (end_r - beg_r + 1));

  for (i = 0; i < end_r - beg_r + 1; i++)
    for (j = 0; j < dim_x + 2; j++)
      access[j][i] = INF;
  access[0][0] = dim_x + 2;
  while (fgets(tmp, sizeof(tmp), in) != 0 && --end_r > -1) {
    /* read a record */
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
          offset                      += temp;
          access[count][seq_pos + 5]  = (int)rint(100 * n); /* round the number */
          count++;
        }
        offset = 0;
      }
    }
  }
  fclose(in);
  return access;
}


static void
redraw_output(char  *fname,
              char  *output,
              int   plfold_up_flag,
              char  *suffix,
              char  *access,
              char  *sname,
              char  *tname)
{
  char  *line;

  line = NULL;
  int   count;
  short two_seq = 0;
  char  *results;
  char  *sequence;
  char  *query;
  char  *target;
  char  *structure;
  char  *pos;
  int   posi;
  int   begin = 0, end = 0, u = 0;
  if (output == NULL) {
    output = (char *)vrna_alloc(sizeof(char) * 2);
    strcpy(output, ".\0");
  }

  count = 0;
  if (sname == NULL && tname == NULL) {
    while ((line = vrna_read_line(stdin)) != NULL) {
      count++;
      if (two_seq == 0 && *line == '>') {
        query = (char *)vrna_alloc(strlen(line) + 1);
        (void)sscanf(line, "%s", query);
        free(line);
        line = NULL;
        memmove(query, query + 1, strlen(query));
        two_seq++;
      } else if (two_seq == 1 && *line == '>') {
        target = (char *)vrna_alloc(strlen(line) + 1);
        (void)sscanf(line, "%s", target);
        free(line);
        line = NULL;
        memmove(target, target + 1, strlen(target));
        two_seq++;
      } else if (two_seq == 2) {
        int *relative_access;
        relative_access = NULL;
        printf("%s %s\n", target, query);
        if (strchr(line, '(') && strchr(line, '&') && strchr(line, '(') && strchr(line, ',') && strchr(line, ':') && strchr(line, '-')) {
          results = line;
          int length;
          pos       = strchr(results, ' ');
          posi      = (int)(pos - results);
          length    = posi;
          structure = (char *)vrna_alloc((length + 1) * sizeof(char));
          sscanf(results, "%s", structure); /* parse structure */
          char *line2;
          if ((line2 = vrna_read_line(stdin)) != NULL) {
            sequence = (char *)vrna_alloc((length + 1) * sizeof(char));
            sscanf(line2, "%s", sequence);
            if (line2 != NULL)
              free(line2);
          } else {
            printf("your file is messed up");
            exit(0);
          } /* parse sequence */

          sscanf(pos, "%10d,%10d;%10d", &begin, &end, &u); /* parse coordinates */
          if (plfold_up_flag) {
            int   **access_s1;
            char  *file_s1;
            /* read_rnaup_file */
            file_s1 = (char *)vrna_alloc(sizeof(char) * (strlen(target) + strlen(suffix) + strlen(access) + 3));
            strcpy(file_s1, access);
            strcat(file_s1, "/");
            strcat(file_s1, target);
            strcat(file_s1, "_");
            strcat(file_s1, suffix);
            access_s1 = read_rnaup(file_s1, begin, end);
            free(file_s1);
            relative_access     = (int *)vrna_alloc(sizeof(int) * (end - begin + 2));
            relative_access[0]  = access_s1[1][1 + 5];
            int i;
            for (i = 2; i < (end - begin + 2); i++)
              relative_access[i - 1] = access_s1[i + 1][i + 5] - access_s1[i][i + 4];
            int l = access_s1[0][0];
            while (--l > -1)
              free(access_s1[l]);
            free(access_s1);
          }

          char  *catseq, *catstruct, *output_file;
          catseq    = (char *)vrna_alloc(strlen(sequence) * sizeof(char));
          catstruct = (char *)vrna_alloc(strlen(structure) * sizeof(char));
          int   l1 = strchr(structure, '&') - structure;
          strncpy(catseq, sequence, l1);
          strcat(catseq, sequence + l1 + 1);
          strncpy(catstruct, structure, l1);
          strcat(catstruct, structure + l1 + 1);
          strcat(catseq, "\0");
          strcat(catstruct, "\0");

          /* printf("%s\n%s\n%s\n%s", catseq,sequence,catstruct,structure); */
          cut_point   = l1 + 1;
          output_file = (char *)vrna_alloc((strlen(output) + strlen(query) + strlen(target) + 50) * sizeof(char));
          strcpy(output_file, output);
          strcat(output_file, "/");
          strcat(output_file, "sno_");
          strcat(output_file, query);
          strcat(output_file, "_");
          strcat(output_file, target);
          strcat(output_file, "_u_");
          char str[9];
          sprintf(str, "%d", u);
          strcat(output_file, str);
          strcat(output_file, "_");
          sprintf(str, "%d", count);
          strcat(output_file, str);
          strcat(output_file, ".ps\0");
          PS_rna_plot_snoop_a(catseq, catstruct, output_file, relative_access, NULL);
          if (relative_access)
            free(relative_access);

          free(output_file);
          output_file = NULL;
          free(catseq);
          catseq = NULL;
          free(catstruct);
          catstruct = NULL;
          free(structure);
          free(sequence);
          free(line);
          line = NULL;
        } else if (*line == '>') {
          free(query);
          free(target);
          two_seq = 1;
          query   = (char *)vrna_alloc(sizeof(char) * (strlen(line) + 1));
          (void)sscanf(line, "%s", query);
          free(line);
          memmove(query, query + 1, strlen(query));
        }
      }
    }
    free(target);
    free(query);
  } else if (!(tname == NULL && sname == NULL)) {
    FILE  *sno, *mrna;
    int   i;
    sno = fopen(sname, "r");
    if (sno == NULL)
      printf("%s: Wrong snoRNA file name\n", sname);

    mrna = fopen(tname, "r");
    if (mrna == NULL)
      printf("%s: Wrong target file name\n", tname);

    char  *temp1[MAX_NUM_NAMES], *temp2[MAX_NUM_NAMES], *AS1[MAX_NUM_NAMES], *AS2[MAX_NUM_NAMES], *names1[MAX_NUM_NAMES], *names2[MAX_NUM_NAMES];
    int   n_seq, n_seq2;
    n_seq   = read_clustal(mrna, temp1, names1);
    n_seq2  = read_clustal(sno, temp2, names2);
    if (n_seq != n_seq2) {
      for (i = 0; temp1[i]; i++) {
        free(temp1[i]);
        free(temp2[i]);
      }
      vrna_message_error("unequal number of seqs in alignments");
    }

    for (i = 0; temp1[i]; i++) {
      AS1[i]  = (char *)vrna_alloc((strlen(temp1[i]) + 11) * sizeof(char));
      AS2[i]  = (char *)vrna_alloc((strlen(temp2[i]) + 11) * sizeof(char));
      strcpy(AS1[i], "NNNNN");
      strcpy(AS2[i], "NNNNN");
      strcat(AS1[i], temp1[i]);
      strcat(AS2[i], temp2[i]);
      strcat(AS1[i], "NNNNN\0");
      strcat(AS2[i], "NNNNN\0");
    }
    for (i = 0; temp1[i]; i++) {
      free(temp1[i]);
      free(temp2[i]);
    }
    AS1[n_seq]  = NULL;
    AS2[n_seq]  = NULL;
    int count = 0;
    while ((line = vrna_read_line(stdin)) != NULL) {
      results = line;
      if (strchr(line, '(') && strchr(line, '&') && strchr(line, '(') && strchr(line, ',') && strchr(line, ':') && strchr(line, '-')) {
        count++;
        int     length;
        /* int sbegin, send; */
        /* int energy; */
        snoopT  dup;
        pos           = strchr(results, ' ');
        posi          = (int)(pos - results);
        length        = posi;
        dup.structure = (char *)vrna_alloc((length + 1) * sizeof(char));
        sscanf(results, "%s %10d,%10d;%10d : %10d,%10d (%10f = %10f + %10f + %10f + %10f + 4.1; duplex cov = %10f; stem cov = %10f",
               dup.structure,
               &begin,
               &(dup.i),
               &(dup.u),
               &(dup.j),
               &end,
               &(dup.energy),
               &(dup.Duplex_El),
               &(dup.Duplex_Er),
               &(dup.Loop_E),
               &(dup.Loop_D),
               &(dup.pscd),
               &(dup.psct)); /* parse duplex stuff; */
        dup.energy    *= n_seq;
        dup.Duplex_El *= n_seq;
        dup.Duplex_Er *= n_seq;
        dup.Loop_E    *= n_seq;
        dup.Loop_D    *= n_seq;
        aliprint_struc(&dup, (const char **)AS1, (const char **)AS2, names1, names2, count, 1);
        free(dup.structure);
        free(line);
      } else {
        free(line);
      }
    }

    for (i = 0; AS1[i]; i++) {
      free(AS1[i]);
      free(AS2[i]);
      free(names1[i]);
      free(names2[i]);
    }
    fclose(mrna);
    fclose(sno);
  }

  /*   free(output); */
}
