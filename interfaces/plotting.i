/*#######################################*/
/* Get coordinates for xy plot           */
/*#######################################*/

%{
  COORDINATE *get_xy_coordinates(const char *structure){
    int i;
    short *table = vrna_ptable(structure);
    short length = (short) strlen(structure);

    COORDINATE *coords = (COORDINATE *) vrna_alloc((length+1)*sizeof(COORDINATE));
    float *X = (float *) vrna_alloc((length+1)*sizeof(float));
    float *Y = (float *) vrna_alloc((length+1)*sizeof(float));

    switch(rna_plot_type){
      case VRNA_PLOT_TYPE_SIMPLE:   simple_xy_coordinates(table, X, Y);
                                    break;
      case VRNA_PLOT_TYPE_CIRCULAR: simple_circplot_coordinates(table, X, Y);
                                    break;
      default:                      naview_xy_coordinates(table, X, Y);
                                    break;
    }

    for(i=0;i<=length;i++){
      coords[i].X = X[i];
      coords[i].Y = Y[i];
    }
    free(table);
    free(X);
    free(Y);
    return(coords);
  }
%}

COORDINATE *get_xy_coordinates(const char *structure);

%extend COORDINATE {
  COORDINATE *get(int i) {
    return self+i;
  }

}


%rename (simple_xy_coordinates) my_simple_xy_coordinates;
%rename (simple_circplot_coordinates) my_simple_circplot_coordinates;
%rename (naview_xy_coordinates) my_naview_xy_coordinates;

%{
#include <vector>
#include <string>

  std::vector<COORDINATE> my_simple_xy_coordinates(std::string structure){
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

  std::vector<COORDINATE> my_simple_circplot_coordinates(std::string structure){
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

  std::vector<COORDINATE> my_naview_xy_coordinates(std::string structure){
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

%}

std::vector<COORDINATE> my_simple_xy_coordinates(std::string);
std::vector<COORDINATE> my_simple_circplot_coordinates(std::string);
std::vector<COORDINATE> my_naview_xy_coordinates(std::string);

%ignore simple_xy_cordinates;
%ignore simple_circplot_coordinates;
%ignore naview_xy_coordinates;

%include <ViennaRNA/plot_layouts.h>




%rename (my_PS_rna_plot_snoop_a) PS_rna_plot_snoop_a;

%{
  int my_PS_rna_plot_snoop_a( std::string sequence,
                              std::string structure,
                              std::string filename,
                              std::vector<int> relative_access,
                              std::vector<std::string> seqs)
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

  int file_PS_rnaplot(std::string sequence,
                      std::string structure,
                      std::string filename,
                      vrna_md_t  *md_p = NULL)
  {
    return vrna_file_PS_rnaplot(sequence.c_str(), structure.c_str(), filename.c_str(), md_p);
  }

  int file_PS_rnaplot_a(std::string sequence,
                        std::string structure,
                        std::string filename,
                        std::string pre,
                        std::string post,
                        vrna_md_t  *md_p = NULL)
  {
    return vrna_file_PS_rnaplot_a(sequence.c_str(), structure.c_str(), filename.c_str(), pre.c_str(), post.c_str(), md_p);
  }

%}

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

%include <ViennaRNA/plot_structure.h>


/*#######################################*/
/* Create colored alignment plots        */
/*#######################################*/

%rename (PS_color_aln)  my_PS_color_aln;
%rename (aliPS_color_aln) my_aliPS_color_aln;

%{
  int my_PS_color_aln(std::string structure,
                      std::string filename,
                      std::vector<std::string> alignment,
                      std::vector<std::string> identifiers)
  {
    std::vector<const char*> aln_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(aln_vec), convert_vecstring2veccharcp);
    aln_vec.push_back(NULL); /* mark end of sequences */
    std::vector<const char*> id_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(id_vec), convert_vecstring2veccharcp);
    id_vec.push_back(NULL); /* mark end of sequences */

    return vrna_file_PS_aln(filename.c_str(), (const char **)&aln_vec[0], (const char **)&id_vec[0], structure.c_str(), 60);
  }

  int my_aliPS_color_aln( std::string structure,
                          std::string filename,
                          std::vector<std::string> alignment,
                          std::vector<std::string> identifiers)
  {
    std::vector<const char*> aln_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(aln_vec), convert_vecstring2veccharcp);
    aln_vec.push_back(NULL); /* mark end of sequences */
    std::vector<const char*> id_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(id_vec), convert_vecstring2veccharcp);
    id_vec.push_back(NULL); /* mark end of sequences */

    return aliPS_color_aln(structure.c_str(), filename.c_str(), (const char **)&aln_vec[0], (const char **)&id_vec[0]);
  }

  int file_PS_aln(std::string filename,
                  std::vector<std::string> alignment,
                  std::vector<std::string> identifiers,
                  std::string structure,
                  int columns = 60)
  {
    std::vector<const char*> aln_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(aln_vec), convert_vecstring2veccharcp);
    aln_vec.push_back(NULL); /* mark end of sequences */
    std::vector<const char*> id_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(id_vec), convert_vecstring2veccharcp);
    id_vec.push_back(NULL); /* mark end of sequences */

    return vrna_file_PS_aln(filename.c_str(), (const char **)&aln_vec[0], (const char **)&id_vec[0], structure.c_str(), columns);
  }

  int file_PS_aln_sub(std::string filename,
                      std::vector<std::string> alignment,
                      std::vector<std::string> identifiers,
                      std::string structure,
                      int start,
                      int end,
                      int columns = 60)
  {
    std::vector<const char*> aln_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(aln_vec), convert_vecstring2veccharcp);
    aln_vec.push_back(NULL); /* mark end of sequences */
    std::vector<const char*> id_vec;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(id_vec), convert_vecstring2veccharcp);
    id_vec.push_back(NULL); /* mark end of sequences */

    return vrna_file_PS_aln_sub(filename.c_str(), (const char **)&aln_vec[0], (const char **)&id_vec[0], structure.c_str(), start, end, columns);
  }

%}

int my_PS_color_aln(std::string structure,
                    std::string filename,
                    std::vector<std::string> alignment,
                    std::vector<std::string> identifiers);

int my_aliPS_color_aln( std::string structure,
                        std::string filename,
                        std::vector<std::string> alignment,
                        std::vector<std::string> identifiers);

int file_PS_aln(std::string filename,
                std::vector<std::string> alignment,
                std::vector<std::string> identifiers,
                std::string structure,
                int columns);

int file_PS_aln(std::string filename,
                std::vector<std::string> alignment,
                std::vector<std::string> identifiers,
                std::string structure);

int file_PS_aln_sub(std::string filename,
                    std::vector<std::string> alignment,
                    std::vector<std::string> identifiers,
                    std::string structure,
                    int start,
                    int end,
                    int columns);

int file_PS_aln_sub(std::string filename,
                    std::vector<std::string> alignment,
                    std::vector<std::string> identifiers,
                    std::string structure,
                    int start,
                    int end);

%ignore PS_color_aln;
%ignore aliPS_color_aln;

%include <ViennaRNA/plot_aln.h>


%include <ViennaRNA/PS_dot.h>

/*#################################*/
/* END get coordinates for xy plot */
/*#################################*/

