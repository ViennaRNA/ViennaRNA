/*
 *  ViennaRNA/utils/structure_tree.c
 *
 *  Various functions for tree representation of secondary strucures
 *
 *  c  Ivo L Hofacker, Walter Fontana, Ronny Lorenz
 *     ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/structures/tree.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE char *
annotate_enclosing_pairs(const char *structure);


PRIVATE char *
db2HIT(const char *structure);


PRIVATE char *
db2Shapiro(const char *structure,
           int        with_stems,
           int        with_weights,
           int        with_external);


PRIVATE char *
db2ExpandedTree(const char *structure);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_db_to_tree_string(const char   *structure,
                       unsigned int type)
{
  char *tree = NULL;

  if (structure) {
    switch (type) {
      case VRNA_STRUCTURE_TREE_HIT:
        tree = db2HIT(structure);
        break;

      case VRNA_STRUCTURE_TREE_SHAPIRO_SHORT:
        tree = db2Shapiro(structure, 0, 0, 0);
        break;

      case VRNA_STRUCTURE_TREE_SHAPIRO:
        tree = db2Shapiro(structure, 1, 0, 0);
        break;

      case VRNA_STRUCTURE_TREE_SHAPIRO_EXT:
        tree = db2Shapiro(structure, 1, 0, 1);
        break;

      case VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT:
        tree = db2Shapiro(structure, 1, 1, 1);
        break;

      case VRNA_STRUCTURE_TREE_EXPANDED:
        tree = db2ExpandedTree(structure);
        break;

      default:
        break;
    }
  }

  return tree;
}


PUBLIC char *
vrna_tree_string_unweight(const char *structure)
{
  unsigned int  i, l;
  char          *tree;

  tree = NULL;

  if (structure) {
    tree = (char *)vrna_alloc(strlen(structure) + 1);

    i = l = 0;
    while (structure[i]) {
      if (!isdigit((int)structure[i]))
        tree[l++] = structure[i];

      i++;
    }
    tree[l] = '\0';
    /* resize to actual length */
    tree = (char *)vrna_realloc(tree, sizeof(char) * (l + 1));
  }

  return tree;
}


PUBLIC char *
vrna_tree_string_to_db(const char *structure)
{
  unsigned int  i, j, k, l, n, w, *match_paren;
  char          id[10], *db;
  const char    *temp;
  int           o;
  vrna_cstr_t   tmp_struct;

  db = NULL;

  if (structure) {
    n           = strlen(structure);
    tmp_struct  = vrna_cstr(4 * n, NULL);
    match_paren = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n / 2 + 1));

    i     = n - 1;
    o     = 0;
    k     = 9;
    id[k] = '\0';

    while (1) {
      switch (structure[i]) {
        case '(':
          if (o < 0) {
            vrna_log_warning("vrna_tree_string_to_db(): "
                                 "Unbalanced parenthesis detected in tree string!"
                                 "Can't convert back to dot-bracket notation");
            goto tree_string_to_db_exit;
          }

          for (j = 0; j < match_paren[o]; j++)
            vrna_cstr_printf(tmp_struct, "(");

          match_paren[o--] = 0;
          break;

        case 'U':
          w = 1;
          sscanf(id + k, "%9u", &w);
          for (j = 0; j < w; j++)
            vrna_cstr_printf(tmp_struct, ".");

          k = 9;
          break;

        case 'P':
          w = 1;
          sscanf(id + k, "%9u", &w);
          for (j = 0; j < w; j++)
            vrna_cstr_printf(tmp_struct, ")");

          match_paren[o]  = w;
          k               = 9;
          break;

        case 'R':
          break;

        case ')':
          o++;
          break;

        default:
          if (!isdigit((int)structure[i])) {
            vrna_log_warning("vrna_tree_string_to_db(): "
                                 "Unsupported character \"%c\" detected in tree string! "
                                 "Can't convert back to dot-bracket notation",
                                 structure[i]);
            goto tree_string_to_db_exit;
          }

          if (k == 0) {
            vrna_log_warning("vrna_tree_string_unexpand(): "
                                 "Node weight too large! "
                                 "Can't convert back to dot-bracket notation");
            goto tree_string_to_db_exit;
          }

          id[--k] = structure[i];
          break;
      }

      if (i == 0)
        break;

      i--;
    }

    temp  = vrna_cstr_string(tmp_struct);
    l     = strlen(temp);
    db    = (char *)vrna_alloc(sizeof(char) * (l + 1));

    for (i = 0; i < l; i++)
      db[i] = temp[l - i - 1];

    db[l] = '\0';

tree_string_to_db_exit:

    vrna_cstr_discard(tmp_struct);
    vrna_cstr_free(tmp_struct);
    free(match_paren);
  }

  return db;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE char *
annotate_enclosing_pairs(const char *structure)
{
  int   *match;
  int   i, n, o, p;
  char  *annot;

  annot = NULL;

  if (structure) {
    n     = (int)strlen(structure);
    annot = strdup(structure);
    match = (int *)vrna_alloc(sizeof(int) * (n / 2 + 1));

    for (i = o = 0; i < n; i++)
      switch (annot[i]) {
        case '.':
          break;

        case '(':
          match[++o] = i;
          break;

        case ')':
          p = i;
          /* seek for outer-most pair in current helix */
          while ((annot[p + 1] == ')') && (match[o - 1] == match[o] - 1)) {
            p++;
            o--;
          }
          /* annotate enclosing pair */
          annot[p]        = ']';
          annot[match[o]] = '[';
          i               = p;
          o--;
          break;

        default:
          vrna_log_warning("annotate_enclosing_pairs: "
                               "Dot-braket string contains junk character \"%c\"",
                               annot[i]);
          free(annot);
          free(match);
          return NULL;
      }

    free(match);
  }

  return annot;
}


PRIVATE char *
db2HIT(const char *structure)
{
  unsigned int  i, p, u, n;
  char          *HIT, *annot_struct;
  vrna_cstr_t   tmp_struct;

  HIT           = NULL;
  annot_struct  = annotate_enclosing_pairs(structure);
  if (annot_struct) {
    n           = strlen(structure);
    tmp_struct  = vrna_cstr(4 * n, NULL);

    /* start of tree */
    vrna_cstr_printf(tmp_struct, "(");

    for (i = p = u = 0; i < n; i++) {
      switch (annot_struct[i]) {
        case '.':
          u++;
          break;

        case '[':
          if (u > 0) {
            vrna_cstr_printf(tmp_struct, "(U%d)", u);
            u = 0;
          }

          vrna_cstr_printf(tmp_struct, "(");
          break;

        case ')':
          if (u > 0) {
            vrna_cstr_printf(tmp_struct, "(U%d)", u);
            u = 0;
          }

          p++;
          break;

        case ']':
          if (u > 0) {
            vrna_cstr_printf(tmp_struct, "(U%d)", u);
            u = 0;
          }

          vrna_cstr_printf(tmp_struct, "P%d)", p + 1);
          p = 0;
          break;
      }
    }

    if (u > 0)
      vrna_cstr_printf(tmp_struct, "(U%d)", u);

    /* add root and end of tree */
    vrna_cstr_printf(tmp_struct, "R)");

    HIT = strdup(vrna_cstr_string(tmp_struct));

    vrna_cstr_discard(tmp_struct);
    vrna_cstr_free(tmp_struct);
    free(annot_struct);
  }

  return HIT;
}


PRIVATE char *
db2Shapiro(const char *structure,
           int        with_stems,
           int        with_weights,
           int        with_external)
{
  unsigned int  i, n, p, lp, loops, pairs, unpaired, *loop_size, *loop, *bulge,
                *loop_degree, *helix_size;
  char          *tree, *annot_struct;
  vrna_cstr_t   tmp_struct;

  tree          = NULL;
  annot_struct  = annotate_enclosing_pairs(structure);

  if (annot_struct) {
    n           = strlen(structure);
    tmp_struct  = vrna_cstr(4 * n, NULL);

    loop_size   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n / 2 + 1));
    helix_size  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n / 2 + 1));
    loop        = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n / 2 + 1));
    bulge       = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n / 2 + 1));
    loop_degree = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n / 2 + 1));

    p = lp = loops = pairs = unpaired = 0;

    for (i = 0; i < n; i++) {
      switch (annot_struct[i]) {
        case '.':
          unpaired++;
          loop_size[loop[lp]]++;
          break;

        case '[':
          vrna_cstr_printf(tmp_struct, "(");

          if (with_stems)
            vrna_cstr_printf(tmp_struct, "(");

          if ((i > 0) && (annot_struct[i - 1] == '(' || annot_struct[i - 1] == '['))
            bulge[lp] = 1;

          lp++;
          loop_degree[++loops]  = 1;
          loop[lp]              = loops;
          bulge[lp]             = 0;
          break;

        case ')':
          if (annot_struct[i - 1] == ']')
            bulge[lp] = 1;

          p++;
          break;

        case ']':
          if (annot_struct[i - 1] == ']')
            bulge[lp] = 1;

          switch (loop_degree[loop[lp]]) {
            case 1:
              vrna_cstr_printf(tmp_struct, "H");    /* hairpin */
              break;

            case 2:
              if (bulge[lp] == 1)
                vrna_cstr_printf(tmp_struct, "B");  /* bulge */
              else
                vrna_cstr_printf(tmp_struct, "I");  /* internal loop */

              break;

            default:
              vrna_cstr_printf(tmp_struct, "M");    /* multiloop */
          }

          helix_size[loop[lp]] = p + 1;

          if (with_weights)
            vrna_cstr_printf(tmp_struct,
                             "%d",
                             loop_size[loop[lp]]);

          vrna_cstr_printf(tmp_struct, ")");

          if (with_stems) {
            vrna_cstr_printf(tmp_struct, "S");
            if (with_weights)
              vrna_cstr_printf(tmp_struct,
                               "%d",
                               helix_size[loop[lp]]);

            vrna_cstr_printf(tmp_struct, ")");
          }

          pairs += p + 1;
          p     = 0;
          loop_degree[loop[--lp]]++;
          break;
      }
    }

    if ((with_external) && (loop_size[0])) {
      if (with_weights)
        tree = vrna_strdup_printf("((%sE%d)R)",
                                  vrna_cstr_string(tmp_struct),
                                  loop_size[0]);
      else
        tree = vrna_strdup_printf("((%sE)R)",
                                  vrna_cstr_string(tmp_struct));
    } else {
      /* add root and end of tree */
      tree = vrna_strdup_printf("(%sR)",
                                vrna_cstr_string(tmp_struct));
    }

    vrna_cstr_discard(tmp_struct);
    vrna_cstr_free(tmp_struct);
    free(loop_degree);
    free(loop_size);
    free(helix_size);
    free(loop);
    free(bulge);
    free(annot_struct);
  }

  return tree;
}


PRIVATE char *
db2ExpandedTree(const char *structure)
{
  unsigned int  i, n;
  char          *tree;
  vrna_cstr_t   tmp_struct;

  n           = strlen(structure);
  tmp_struct  = vrna_cstr(4 * n, NULL);

  for (i = 0; i < n; i++) {
    if (structure[i] == '(')
      vrna_cstr_printf(tmp_struct, "(");
    else if (structure[i] == ')')
      vrna_cstr_printf(tmp_struct, "P)");
    else
      vrna_cstr_printf(tmp_struct, "(U)");
  }

  tree = vrna_strdup_printf("(%sR)",
                            vrna_cstr_string(tmp_struct));

  vrna_cstr_discard(tmp_struct);
  vrna_cstr_free(tmp_struct);

  return tree;
}
