
static void
print_SVG_header(FILE         *fh,
                const float   scale[2],
                const float   transform[2],
                unsigned int  options)
{
  fprintf(fh,
          "%s",
          SVG_structure_plot_header);

  fprintf(fh,
          "  <g transform=\"scale(%7f,%7f) translate(%7f,%7f)\">\n",
          scale[0],
          scale[1],
          transform[0],
          transform[1]);
}


PRIVATE void
print_SVG_footer(FILE *fh)
{
  fprintf(fh,
          "  </g>\n%s",
          SVG_structure_plot_footer);
}


PRIVATE void
print_SVG_bases(FILE          *fh,
                float         *X,
                float         *Y,
                const char    *string,
                unsigned int  n)
{
  unsigned int i;

  fprintf(fh,
          "    <g transform=\"translate(-4.6, 4)\" id=\"seq\">\n");

  for (i = 0; i < n; i++)
    fprintf(fh,
            "      <text class=\"nucleotide\" x=\"%.3f\" y=\"%.3f\">%c</text>\n",
            X[i], Y[i], string[i]);

  fprintf(fh,
          "    </g>\n");
}

PRIVATE void
print_SVG_backbone(FILE         *fh,
                   const float  *X,
                   const float  *Y,
                   unsigned int n)
{
  unsigned int i;

  fprintf(fh,
          "    <polyline class=\"backbone\" id=\"outline\" points=\"\n");

  for (i = 0; i < n; i++)
    fprintf(fh,
            "      %3.3f,%3.3f\n",
            X[i],
            Y[i]);

  fprintf(fh,
          "    \" />\n");
}


PRIVATE void
print_SVG_pairs(FILE          *fh,
                const short   *pt,
                const float   *X,
                const float   *Y,
                const float   *CX,
                const float   *CY,
                unsigned int  n,
                unsigned int  plot_type)
{
  unsigned int i, j;

  fprintf(fh,
          "    <g id=\"pairs\">\n");

  for (i = 1; i <= n; i++) {
    if ((j = (unsigned int)pt[i]) > i) {
      if (plot_type == VRNA_PLOT_TYPE_CIRCULAR) {
        fprintf(fh,
                "      <path class=\"basepairs\" id=\"%u,%u\" d=\"M %6.5f %6.5f C %6.5f,%6.5f %6.5f,%6.5f %6.5f %6.5f\" />\n",
                i,
                j,
                X[i - 1],
                Y[i - 1],
                CX[i - 1],
                CY[i - 1],
                CX[j - 1],
                CY[j - 1],
                X[j - 1],
                Y[j - 1]);
       } else {
        fprintf(fh,
                "      <line class=\"basepairs\" id=\"%u,%u\" x1=\"%6.5f\" y1=\"%6.5f\" x2=\"%6.5f\" y2=\"%6.5f\" />\n",
                i, j,
                X[i - 1], Y[i - 1],
                X[j - 1], Y[j - 1]);
      }
    }
  }
  fprintf(fh,
          "    </g>\n");
}
