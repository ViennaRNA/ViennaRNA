


static void
print_PS_header(FILE        *fh,
                const char  *title,
                int         bbox_x0,
                int         bbox_y0,
                int         bbox_x1,
                int         bbox_y1)
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
          VERSION,
          vrna_time_stamp(),
          title,
          bbox_x0, bbox_y0, bbox_x1, bbox_y1);
}


static void
print_PS_options(FILE      *fh,
                 vrna_md_t *md)
{
  fprintf(fh,
          "%%Options: %s\n",
          (md) ?
            vrna_md_option_string(md) :
            option_string());
}


static void
print_PS_help_structure(FILE *fh)
{
  fprintf(fh, "%% to switch off outline pairs of sequence comment or\n"
          "%% delete the appropriate line near the end of the file\n\n");
}


PRIVATE void
print_PS_sequence(FILE    *fh,
                  const char *sequence)
{
  unsigned int i, length;

  length = strlen(sequence);

  fprintf(fh, "/sequence (\\\n");

  i = 0;
  while (i < length) {
    fprintf(fh, "%.255s\\\n", sequence + i);  /* no lines longer than 255 */
    i += 255;
  }
  fprintf(fh, ") def\n");
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


