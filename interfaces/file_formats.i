/**********************************************/
/* BEGIN interface for reading/writing files  */
/**********************************************/


%apply int    *OUTPUT { int *status };
%apply std::string *OUTPUT { std::string *shape_sequence };

%rename (file_SHAPE_read) my_file_SHAPE_read;

%{

  std::vector<double> my_file_SHAPE_read( std::string file_name,
                                          int length,
                                          double default_value,
                                          std::string *shape_sequence,
                                          int *status){

    std::vector<double> values (length+1, -999);
    char *seq = (char *)vrna_alloc(sizeof(char) * (length + 1));

    *status = vrna_file_SHAPE_read(file_name.c_str(), length, default_value, seq, (double *)&values[0]);

    *shape_sequence = std::string(seq);

    free(seq);
    return values;
  }

%}

std::vector<double> my_file_SHAPE_read( const char *file_name,
                                          int length,
                                          double default_value,
                                          std::string *shape_sequence,
                                          int *status);

%clear int *status;
%clear std::string *shape_sequence;

%include <ViennaRNA/file_formats.h>
