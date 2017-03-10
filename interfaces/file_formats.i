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

/**********************************************/
/* BEGIN interface for reading/writing MSA    */
/**********************************************/

%rename (file_msa_detect_format) vrna_file_msa_detect_format;

%constant unsigned int FILE_FORMAT_MSA_CLUSTAL   = VRNA_FILE_FORMAT_MSA_CLUSTAL;
%constant unsigned int FILE_FORMAT_MSA_DEFAULT   = VRNA_FILE_FORMAT_MSA_DEFAULT;
%constant unsigned int FILE_FORMAT_MSA_FASTA     = VRNA_FILE_FORMAT_MSA_FASTA;
%constant unsigned int FILE_FORMAT_MSA_MAF       = VRNA_FILE_FORMAT_MSA_MAF;
%constant unsigned int FILE_FORMAT_MSA_NOCHECK   = VRNA_FILE_FORMAT_MSA_NOCHECK;
%constant unsigned int FILE_FORMAT_MSA_STOCKHOLM = VRNA_FILE_FORMAT_MSA_STOCKHOLM;
%constant unsigned int FILE_FORMAT_MSA_MIS       = VRNA_FILE_FORMAT_MSA_MIS;
%constant unsigned int FILE_FORMAT_MSA_UNKNOWN   = VRNA_FILE_FORMAT_MSA_UNKNOWN;
%constant unsigned int FILE_FORMAT_MSA_QUIET     = VRNA_FILE_FORMAT_MSA_QUIET;
%constant unsigned int FILE_FORMAT_MSA_SILENT    = VRNA_FILE_FORMAT_MSA_SILENT;

%include <ViennaRNA/file_formats_msa.h>
