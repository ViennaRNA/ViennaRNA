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

%include <ViennaRNA/plot_layouts.h>

%include <ViennaRNA/plot_structure.h>

%include <ViennaRNA/plot_aln.h>


%include <ViennaRNA/PS_dot.h>

/*#################################*/
/* END get coordinates for xy plot */
/*#################################*/

