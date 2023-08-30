
#define   PS_MACRO_LAYOUT_BASE              1U
#define   PS_MACRO_LAYOUT_EXTRAS            2U
#define   PS_MACRO_DOTPLOT_BASE             4U
#define   PS_MACRO_DOTPLOT_SD               8U
#define   PS_MACRO_DOTPLOT_UD              16U
#define   PS_MACRO_DOTPLOT_SC_MOTIFS       32U
#define   PS_MACRO_DOTPLOT_LINEAR_DATA     64U
#define   PS_MACRO_DOTPLOT_ALL             (PS_MACRO_DOTPLOT_LINEAR_DATA | PS_MACRO_DOTPLOT_SD | PS_MACRO_DOTPLOT_UD | PS_MACRO_DOTPLOT_SC_MOTIFS)
#define   PS_MACRO_ALN_BASE               128U

static void
print_PS_header(FILE          *fh,
                const char    *title,
                int           bbox[4],
                vrna_md_t     *md,
                const char    *comments,
                const char    *dict_name,
                unsigned int  options)
{
  fprintf(fh,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Creator: ViennaRNA-%s\n"
          "%%%%CreationDate: %s"
          "%%%%Title: %s\n"
          "%%%%BoundingBox: %d %d %d %d\n"
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n\n",
          VRNA_VERSION,
          vrna_time_stamp(),
          title,
          bbox[0], bbox[1], bbox[2], bbox[3]);

  if (md)
    fprintf(fh, "%% Program options: %s\n\n", vrna_md_option_string(md));

  if (comments) {
    char      *c, *ptr, *save_ptr, *ptr2;
    unsigned  int i, l;
    c   = strdup(comments);

    ptr2  = c;
#ifdef _WIN32
# ifndef __MINGW32__
    ptr   = strtok_s(c, "\n", &save_ptr);
# else
    ptr   = strtok_r(c, "\n", &save_ptr);
# endif
#else
    ptr   = strtok_r(c, "\n", &save_ptr);
#endif

    l     = ptr - ptr2;

    while (ptr != NULL) {
      /* print empty newlines if they were present in input */
      for (i = 1; i < l; i++)
        fprintf(fh, "\n");

      /* print comment line */
      fprintf(fh, "%% %s\n", ptr);

      ptr2  = ptr;
      l     = strlen(ptr);
#ifdef _WIN32
# ifndef __MINGW32__
      ptr   = strtok_s(NULL, "\n", &save_ptr);
# else
      ptr   = strtok_r(NULL, "\n", &save_ptr);
# endif
#else
      ptr   = strtok_r(NULL, "\n", &save_ptr);
#endif

      /* compute number of delimiters between previous and current token */
      l   = ptr - ptr2 - l;
    }
    free(c);

    fprintf(fh, "\n");
  }

  /* create dictionary */
  fprintf(fh, "/%s 100 dict def\n\n", dict_name);

  /* put dictionary on stack top */
  fprintf(fh, "%s begin\n\n%%%%BeginProlog\n\n", dict_name);

  /* add macros */

  if (options & PS_MACRO_LAYOUT_BASE)
    fprintf(fh, "%s", PS_structure_plot_macro_base);

  if (options & PS_MACRO_LAYOUT_EXTRAS)
    fprintf(fh, "%s", PS_structure_plot_macro_extras);

  if (options & PS_MACRO_DOTPLOT_BASE)
    fprintf(fh,"%s", PS_dot_plot_macro_base);

  if (options & PS_MACRO_DOTPLOT_SD)  /* structured domains, a.k.a. gquads et al. */
    fprintf(fh,"%s", PS_dot_plot_macro_sd);

  if (options & PS_MACRO_DOTPLOT_UD)  /* unstructured domains */
    fprintf(fh,"%s", PS_dot_plot_macro_ud);

  if (options & PS_MACRO_DOTPLOT_SC_MOTIFS) /* soft constraint motifs */
    fprintf(fh,"%s", PS_dot_plot_macro_sc_motifs);

  if (options & PS_MACRO_DOTPLOT_LINEAR_DATA) /* linear data */
    fprintf(fh,"%s", PS_dot_plot_macro_linear_data);

  if (options & PS_MACRO_ALN_BASE)
    fprintf(fh, "%s", PS_aln_macro_base);

  /* remove dictionary from stack */
  fprintf(fh, "\n%%%%EndProlog\n\n");
}


PRIVATE void
print_PS_footer(FILE *fh)
{
  fprintf(fh, "showpage\n"
              "end\n"
              "%%%%EOF\n");
}


PRIVATE void
print_PS_sequence(FILE    *fh,
                  const char *sequence)
{
  unsigned int i, length;

  length = strlen(sequence);

  fprintf(fh, "/sequence { (\\\n");

  i = 0;
  while (i < length) {
    fprintf(fh, "%.255s\\\n", sequence + i);  /* no lines longer than 255 */
    i += 255;
  }
  fprintf(fh, ") } def\n");
  fprintf(fh, "/len { sequence length } bind def\n\n");
}


PRIVATE void
print_PS_coords(FILE          *fh,
                float         *X,
                float         *Y,
                unsigned int  n)
{
  unsigned int i;

  fprintf(fh, "/coor [\n");

  for (i = 0; i < n; i++)
    fprintf(fh, "[%3.8f %3.8f]\n", X[i], Y[i]);

  fprintf(fh, "] def\n");
}


