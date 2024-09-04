/*#######################################*/
/* Get coordinates for xy plot           */
/*#######################################*/

%{
  COORDINATE *
  get_xy_coordinates(const char *structure)
  {
    int         i, ret;
    short       length  = (short) strlen(structure);
    COORDINATE  *coords = (COORDINATE *) vrna_alloc((length+1)*sizeof(COORDINATE));
    float       *X, *Y;

    ret = vrna_plot_coords(structure, &X, &Y, rna_plot_type);

    if (ret == (int)length) {
      for(i=0;i<=length;i++){
        coords[i].X = X[i];
        coords[i].Y = Y[i];
      }
    }
    free(X);
    free(Y);

    return(coords);
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") get_xy_coordinates;
%feature("kwargs") get_xy_coordinates;
#endif

COORDINATE *get_xy_coordinates(const char *structure);

%extend COORDINATE {
  COORDINATE *
  get(int i)
  {
    return self+i;
  }

}


%rename (simple_xy_coordinates) my_simple_xy_coordinates;
%rename (simple_circplot_coordinates) my_simple_circplot_coordinates;

#ifdef VRNA_WITH_NAVIEW_LAYOUT
%rename (naview_xy_coordinates) my_naview_xy_coordinates;
#endif

%{
#include <vector>
#include <string>

  std::vector<COORDINATE>
  my_simple_xy_coordinates(std::string structure)
  {
    std::vector<COORDINATE> ret;
    short *table  = vrna_ptable(structure.c_str());
    float *X      = (float *) vrna_alloc((table[0]+1)*sizeof(float));
    float *Y      = (float *) vrna_alloc((table[0]+1)*sizeof(float));
    simple_xy_coordinates(table, X, Y);

    for(int i = 0; i <= table[0]; i++){
      COORDINATE c;
      c.X = X[i];
      c.Y = Y[i];
      ret.push_back(c);
    }

    free(X);
    free(Y);
    free(table);
    return ret;
  }

  std::vector<COORDINATE>
  my_simple_circplot_coordinates(std::string structure)
  {
    std::vector<COORDINATE> ret;
    short *table  = vrna_ptable(structure.c_str());
    float *X      = (float *) vrna_alloc((table[0]+1)*sizeof(float));
    float *Y      = (float *) vrna_alloc((table[0]+1)*sizeof(float));
    simple_circplot_coordinates(table, X, Y);

    for(int i = 0; i <= table[0]; i++){
      COORDINATE c;
      c.X = X[i];
      c.Y = Y[i];
      ret.push_back(c);
    }

    free(X);
    free(Y);
    free(table);
    return ret;
  }

#ifdef VRNA_WITH_NAVIEW_LAYOUT
  std::vector<COORDINATE>
  my_naview_xy_coordinates(std::string structure)
  {
    std::vector<COORDINATE> ret;
    short *table  = vrna_ptable(structure.c_str());
    float *X      = (float *) vrna_alloc((table[0]+1)*sizeof(float));
    float *Y      = (float *) vrna_alloc((table[0]+1)*sizeof(float));
    naview_xy_coordinates(table, X, Y);

    for(int i = 0; i <= table[0]; i++){
      COORDINATE c;
      c.X = X[i];
      c.Y = Y[i];
      ret.push_back(c);
    }

    free(X);
    free(Y);
    free(table);
    return ret;
  }
#endif

%}

#ifdef SWIGPYTHON
%feature("autodoc") my_simple_circplot_coordinates;
%feature("kwargs") my_simple_circplot_coordinates;
# ifdef VRNA_WITH_NAVIEW_LAYOUT
%feature("autodoc") my_naview_xy_coordinates;
%feature("kwargs") my_naview_xy_coordinates;
# endif
#endif

std::vector<COORDINATE> my_simple_xy_coordinates(std::string);
std::vector<COORDINATE> my_simple_circplot_coordinates(std::string);

#ifdef VRNA_WITH_NAVIEW_LAYOUT
std::vector<COORDINATE> my_naview_xy_coordinates(std::string);
#endif

%ignore simple_xy_cordinates;
%ignore simple_circplot_coordinates;

#ifdef VRNA_WITH_NAVIEW_LAYOUT
%ignore naview_xy_coordinates;
#endif

/* class definitions for vrna_plot_layout_t */
%rename (plot_layout) vrna_plot_layout_t;

%nodefaultctor  vrna_plot_layout_t;
%nodefaultdtor  vrna_plot_layout_t;

typedef struct {} vrna_plot_layout_t;

%extend vrna_plot_layout_t {
  vrna_plot_layout_t(std::string  structure,
                     unsigned int plot_type = VRNA_PLOT_TYPE_DEFAULT)
  {
    return vrna_plot_layout(structure.c_str(),
                            plot_type);
  }

  ~vrna_plot_layout_t()
  {
    vrna_plot_layout_free($self);
  }
};

/* class definitions for vrna_plot_options_puzzler_t */
%rename (plot_options_puzzler)  vrna_plot_options_puzzler_t;

%nodefaultctor  vrna_plot_options_puzzler_t;
%nodefaultdtor  vrna_plot_options_puzzler_t;

typedef struct {
  short       checkAncestorIntersections;
  short       checkSiblingIntersections;
  short       checkExteriorIntersections;
  short       allowFlipping;
  short       optimize;
} vrna_plot_options_puzzler_t;

%extend vrna_plot_options_puzzler_t {
  vrna_plot_options_puzzler_t *
  vrna_plot_options_puzzler_t(void)
  {
    return vrna_plot_options_puzzler();
  }

  ~vrna_plot_options_puzzler_t()
  {
    vrna_plot_options_puzzler_free($self);
  }
};



%newobject  vrna_plot_layout_simple;
#ifdef VRNA_WITH_NAVIEW_LAYOUT
%newobject  vrna_plot_layout_naview;
#endif
%newobject  vrna_plot_layout_circular;
%newobject  vrna_plot_layout_turtle;
%newobject  vrna_plot_layout_puzzler;

%rename (plot_layout_simple)    vrna_plot_layout_simple;
#ifdef VRNA_WITH_NAVIEW_LAYOUT
%rename (plot_layout_naview)    vrna_plot_layout_naview;
#endif
%rename (plot_layout_circular)  vrna_plot_layout_circular;
%rename (plot_layout_turtle)    vrna_plot_layout_turtle;
%rename (plot_layout_puzzler)   vrna_plot_layout_puzzler;


%constant unsigned int PLOT_TYPE_SIMPLE   = VRNA_PLOT_TYPE_SIMPLE;
%constant unsigned int PLOT_TYPE_NAVIEW   = VRNA_PLOT_TYPE_NAVIEW;
%constant unsigned int PLOT_TYPE_CIRCULAR = VRNA_PLOT_TYPE_CIRCULAR;
%constant unsigned int PLOT_TYPE_TURTLE   = VRNA_PLOT_TYPE_TURTLE;
%constant unsigned int PLOT_TYPE_PUZZLER  = VRNA_PLOT_TYPE_PUZZLER;
%constant unsigned int PLOT_TYPE_DEFAULT  = VRNA_PLOT_TYPE_DEFAULT;


%include <ViennaRNA/plotting/layouts.h>




%rename (my_PS_rna_plot_snoop_a) PS_rna_plot_snoop_a;

%{
  int
  my_PS_rna_plot_snoop_a( std::string               sequence,
                          std::string               structure,
                          std::string               filename,
                          std::vector<int>          relative_access,
                          std::vector<std::string>  seqs)
  {
    std::vector<const char*> seqs_vec;
    std::transform(seqs.begin(), seqs.end(), std::back_inserter(seqs_vec), convert_vecstring2veccharcp);
    seqs_vec.push_back(NULL); /* mark end of sequences */

    return PS_rna_plot_snoop_a( sequence.c_str(),
                                structure.c_str(),
                                filename.c_str(),
                                &relative_access[0],
                                (const char **)&seqs_vec[0]);
  }

  int
  file_PS_rnaplot(std::string sequence,
                  std::string structure,
                  std::string filename,
                  vrna_md_t   *md_p = NULL)
  {
    return vrna_file_PS_rnaplot(sequence.c_str(), structure.c_str(), filename.c_str(), md_p);
  }

  int
  file_PS_rnaplot_a(std::string sequence,
                    std::string structure,
                    std::string filename,
                    std::string pre,
                    std::string post,
                    vrna_md_t   *md_p = NULL)
  {
    return vrna_file_PS_rnaplot_a(sequence.c_str(), structure.c_str(), filename.c_str(), pre.c_str(), post.c_str(), md_p);
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") my_PS_rna_plot_snoop_a;
%feature("kwargs") my_PS_rna_plot_snoop_a;
#endif

int my_PS_rna_plot_snoop_a( std::string sequence,
                            std::string structure,
                            std::string filename,
                            std::vector<int> relative_access,
                            std::vector<std::string> seqs);

int file_PS_rnaplot(std::string sequence,
                      std::string structure,
                      std::string filename,
                      vrna_md_t  *md_p);

int file_PS_rnaplot(std::string sequence,
                      std::string structure,
                      std::string filename);

int file_PS_rnaplot_a(std::string sequence,
                        std::string structure,
                        std::string filename,
                        std::string pre,
                        std::string post,
                        vrna_md_t  *md_p);

int file_PS_rnaplot_a(std::string sequence,
                        std::string structure,
                        std::string filename,
                        std::string pre,
                        std::string post);

%ignore PS_rna_plot_snoop_a;


%rename (plot_data) vrna_plot_data_t;

%nodefaultctor  vrna_plot_data_t;
%nodefaultdtor  vrna_plot_data_t;

typedef struct {
  char          *pre;
  char          *post;
  vrna_md_t     *md;
  unsigned int  options;
} vrna_plot_data_t;

%extend vrna_plot_data_t {
  vrna_plot_data_t(std::string  pre = "",
                   std::string  post = "",
                   vrna_md_t    *md = NULL,
                   unsigned int options = 0)
  {
    vrna_plot_data_t  *d = NULL;

    if ((pre.size() > 0) ||
        (post.size() > 0) ||
        (md)) {
      d           = (vrna_plot_data_t *)vrna_alloc(sizeof(vrna_plot_data_t));
      d->pre      = (pre.size() > 0) ? strdup(pre.c_str()) : NULL;
      d->post     = (post.size() > 0) ? strdup(post.c_str()) : NULL;
      d->md       = md;
      d->options  = options;
    }

    return d;
  }

  ~vrna_plot_data_t()
  {
    free($self->pre);
    free($self->post);
    free($self);
  }
};

%rename (plot_structure)      my_plot_structure;
%rename (plot_structure_svg)  my_plot_structure_svg;
%rename (plot_structure_eps)  my_plot_structure_eps;
%rename (plot_structure_gml)  my_plot_structure_gml;
%rename (plot_structure_ssv)  my_plot_structure_ssv;
%rename (plot_structure_xrna) my_plot_structure_xrna;


#ifdef SWIGPYTHON
%feature("autodoc") my_plot_structure;
%feature("kwargs")  my_plot_structure;
%feature("autodoc") my_plot_structure_svg;
%feature("kwargs")  my_plot_structure_svg;
%feature("autodoc") my_plot_structure_eps;
%feature("kwargs")  my_plot_structure_eps;
%feature("autodoc") my_plot_structure_gml;
%feature("kwargs")  my_plot_structure_gml;
%feature("autodoc") my_plot_structure_ssv;
%feature("kwargs")  my_plot_structure_ssv;
%feature("autodoc") my_plot_structure_xrna;
%feature("kwargs")  my_plot_structure_xrna;
#endif

%{
  int
  my_plot_structure(std::string         filename,
                    std::string         sequence,
                    std::string         structure,
                    unsigned int        file_format = VRNA_FILE_FORMAT_PLOT_DEFAULT,
                    vrna_plot_layout_t  *layout     = NULL,
                    vrna_plot_data_t    *data       = NULL)
  {
    return vrna_plot_structure(filename.c_str(),
                               sequence.c_str(),
                               structure.c_str(),
                               file_format,
                               layout,
                               data);
  }


  int
  my_plot_structure_eps(std::string         filename,
                        std::string         sequence,
                        std::string         structure,
                        vrna_plot_layout_t  *layout = NULL,
                        vrna_plot_data_t    *data   = NULL)
  {
    return vrna_plot_structure_eps(filename.c_str(),
                                   sequence.c_str(),
                                   structure.c_str(),
                                   layout,
                                   data);
  }


  int
  my_plot_structure_svg(std::string         filename,
                        std::string         sequence,
                        std::string         structure,
                        vrna_plot_layout_t  *layout = NULL,
                        vrna_plot_data_t    *data   = NULL)
  {
    return vrna_plot_structure_svg(filename.c_str(),
                                 sequence.c_str(),
                                 structure.c_str(),
                                 layout,
                                 data);
  }


  int
  my_plot_structure_gml(std::string         filename,
                        std::string         sequence,
                        std::string         structure,
                        vrna_plot_layout_t  *layout = NULL,
                        vrna_plot_data_t    *data   = NULL,
                        char                option  = 'x')
  {
    return vrna_plot_structure_gml(filename.c_str(),
                                   sequence.c_str(),
                                   structure.c_str(),
                                   layout,
                                   data,
                                   option);
  }


  int
  my_plot_structure_ssv(std::string         filename,
                        std::string         sequence,
                        std::string         structure,
                        vrna_plot_layout_t  *layout = NULL,
                        vrna_plot_data_t    *data   = NULL)
  {
    return vrna_plot_structure_ssv(filename.c_str(),
                                   sequence.c_str(),
                                   structure.c_str(),
                                   layout,
                                   data);
  }


  int
  my_plot_structure_xrna(std::string        filename,
                         std::string        sequence,
                         std::string        structure,
                         vrna_plot_layout_t *layout = NULL,
                         vrna_plot_data_t   *data   = NULL)
  {
    return vrna_plot_structure_xrna(filename.c_str(),
                                    sequence.c_str(),
                                    structure.c_str(),
                                    layout,
                                    data);
  }


%}


int
my_plot_structure(std::string         filename,
                  std::string         sequence,
                  std::string         structure,
                  unsigned int        file_format = VRNA_FILE_FORMAT_PLOT_DEFAULT,
                  vrna_plot_layout_t  *layout     = NULL,
                  vrna_plot_data_t    *data       = NULL);


int
my_plot_structure_eps(std::string         filename,
                      std::string         sequence,
                      std::string         structure,
                      vrna_plot_layout_t  *layout = NULL,
                      vrna_plot_data_t    *data   = NULL);


int
my_plot_structure_svg(std::string         filename,
                      std::string         sequence,
                      std::string         structure,
                      vrna_plot_layout_t  *layout = NULL,
                      vrna_plot_data_t    *data   = NULL);


int
my_plot_structure_gml(std::string         filename,
                      std::string         sequence,
                      std::string         structure,
                      vrna_plot_layout_t  *layout = NULL,
                      vrna_plot_data_t    *data   = NULL,
                      char                option  = 'x');


int
my_plot_structure_ssv(std::string         filename,
                      std::string         sequence,
                      std::string         structure,
                      vrna_plot_layout_t  *layout = NULL,
                      vrna_plot_data_t    *data   = NULL);


int
my_plot_structure_xrna(std::string        filename,
                       std::string        sequence,
                       std::string        structure,
                       vrna_plot_layout_t *layout = NULL,
                       vrna_plot_data_t   *data   = NULL);


%constant unsigned int FILE_FORMAT_EPS          = VRNA_FILE_FORMAT_EPS;
%constant unsigned int FILE_FORMAT_SVG          = VRNA_FILE_FORMAT_SVG;        
%constant unsigned int FILE_FORMAT_GML          = VRNA_FILE_FORMAT_GML;        
%constant unsigned int FILE_FORMAT_SSV          = VRNA_FILE_FORMAT_SSV;
%constant unsigned int FILE_FORMAT_XRNA         = VRNA_FILE_FORMAT_XRNA;        
%constant unsigned int FILE_FORMAT_PLOT_DEFAULT = VRNA_FILE_FORMAT_PLOT_DEFAULT;


%include <ViennaRNA/plotting/structures.h>


/*#######################################*/
/* Create colored alignment plots        */
/*#######################################*/

#ifdef SWIGPYTHON
%feature("autodoc") file_PS_aln;
%feature("kwargs") file_PS_aln;
#endif

%{
  int
  file_PS_aln(std::string               filename,
              std::vector<std::string>  alignment,
              std::vector<std::string>  identifiers,
              std::string               structure,
              unsigned int              start       = 0,
              unsigned int              end         = 0,
              int                       offset      = 0,
              unsigned int              columns     = 60)
  {
    std::vector<const char*> aln_vec;
    std::vector<const char*> id_vec;

    std::transform(alignment.begin(),
                   alignment.end(),
                   std::back_inserter(aln_vec),
                   convert_vecstring2veccharcp);

    std::transform(alignment.begin(),
                   alignment.end(),
                   std::back_inserter(id_vec),
                   convert_vecstring2veccharcp);

    aln_vec.push_back(NULL); /* mark end of sequences */
    id_vec.push_back(NULL); /* mark end of sequences */

    return vrna_file_PS_aln_slice(filename.c_str(),
                                  (const char **)&aln_vec[0],
                                  (const char **)&id_vec[0],
                                  structure.c_str(),
                                  start,
                                  end,
                                  offset,
                                  columns);
  }

%}

int file_PS_aln(std::string               filename,
                std::vector<std::string>  alignment,
                std::vector<std::string>  identifiers,
                std::string               structure,
                unsigned int              start   = 0,
                unsigned int              end     = 0,
                int                       offset  = 0,
                unsigned int              columns = 60);


%ignore PS_color_aln;
%ignore aliPS_color_aln;
%ignore vrna_file_PS_aln;
%ignore vrna_file_PS_aln_slice;

%include <ViennaRNA/plotting/alignments.h>


%{
  int
  plot_dp_EPS(std::string             filename,
              std::string             sequence,
              std::vector<vrna_ep_t>  upper_triangle,
              std::vector<vrna_ep_t>  lower_triangle = std::vector<vrna_ep_t>(),
              vrna_dotplot_auxdata_t  *auxdata = NULL,
              unsigned int            options = VRNA_PLOT_PROBABILITIES_DEFAULT)
  {
    vrna_ep_t               pp;

    if ((sequence.length() > 0) &&
        (filename.length() > 0)) {
      pp.i = pp.j = 0;
      pp.p = 0.;
      pp.type = VRNA_PLIST_TYPE_BASEPAIR;
      upper_triangle.push_back(pp);
      lower_triangle.push_back(pp);

      return vrna_plot_dp_EPS(filename.c_str(),
                              sequence.c_str(),
                              &upper_triangle[0],
                              &lower_triangle[0],
                              auxdata,
                              options);
    }

    return 0;
  }
%}

int
plot_dp_EPS(std::string             filename,
            std::string             sequence,
            std::vector<vrna_ep_t>  upper_triangle,
            std::vector<vrna_ep_t>  lower_triangle = std::vector<vrna_ep_t>(),
            vrna_dotplot_auxdata_t  *auxdata = NULL,
            unsigned int            options = VRNA_PLOT_PROBABILITIES_DEFAULT);


%include <ViennaRNA/plotting/probabilities.h>

/*#################################*/
/* END get coordinates for xy plot */
/*#################################*/

