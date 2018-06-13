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

%include <ViennaRNA/io/file_formats.h>

/**********************************************/
/* BEGIN interface for reading/writing MSA    */
/**********************************************/

%rename (file_msa_detect_format)  my_file_msa_detect_format;
%rename (file_msa_read)           my_file_msa_read;
%rename (file_msa_read_record)    my_file_msa_read_record;
%rename (file_msa_write)          my_file_msa_write;

%{

  unsigned int
  my_file_msa_detect_format(std::string   filename,
                            unsigned int  options = VRNA_FILE_FORMAT_MSA_DEFAULT)
  {
    return vrna_file_msa_detect_format(filename.c_str(), options);
  }


  int my_file_msa_read( std::string               filename,
                        std::vector<std::string>  *names,
                        std::vector<std::string>  *alignment,
                        std::string               *id,
                        std::string               *structure,
                        unsigned int              options = VRNA_FILE_FORMAT_MSA_STOCKHOLM)
  {
    char **c_names, **c_aln, *c_id, *c_structure;

    int ret = vrna_file_msa_read(filename.c_str(), &c_names, &c_aln, &c_id, &c_structure, options);

    if (ret != -1) {
      names->clear();
      alignment->clear();

      names->reserve(ret);
      alignment->reserve(ret);

      for (int i = 0; i < ret; i++) {
        std::string id(c_names[i]);
        std::string seq(c_aln[i]);
        names->push_back(id);
        alignment->push_back(seq);
        free(c_names[i]);
        free(c_aln[i]);
      }
      *id        = (c_id) ? c_id : "";
      *structure = (c_structure) ? c_structure : "";

      free(c_names);
      free(c_aln);
      free(c_id);
      free(c_structure);
    }

    return ret;
  }


  int my_file_msa_read_record(FILE                      *filehandle,
                              std::vector<std::string>  *names,
                              std::vector<std::string>  *alignment,
                              std::string               *id,
                              std::string               *structure,
                              unsigned int              options = VRNA_FILE_FORMAT_MSA_STOCKHOLM)
  {
    char **c_names, **c_aln, *c_id, *c_structure;

    int ret = vrna_file_msa_read_record(filehandle, &c_names, &c_aln, &c_id, &c_structure, options);

    if (ret != -1) {
      names->clear();
      alignment->clear();

      names->reserve(ret);
      alignment->reserve(ret);

      for (int i = 0; i < ret; i++) {
        std::string id(c_names[i]);
        std::string seq(c_aln[i]);
        names->push_back(id);
        alignment->push_back(seq);
        free(c_names[i]);
        free(c_aln[i]);
      }
      *id        = (c_id) ? c_id : "";
      *structure = (c_structure) ? c_structure : "";

      free(c_names);
      free(c_aln);
      free(c_id);
      free(c_structure);
    }

    return ret;
  }


  int my_file_msa_write(std::string               filename,
                        std::vector<std::string>  names,
                        std::vector<std::string>  alignment,
                        std::string               id = "",
                        std::string               structure = "",
                        std::string               source = "",
                        unsigned int              options = VRNA_FILE_FORMAT_MSA_STOCKHOLM | VRNA_FILE_FORMAT_MSA_APPEND)
  {
    std::vector<const char*>  v_names;
    std::vector<const char*>  v_aln;

    transform(names.begin(), names.end(), back_inserter(v_names), convert_vecstring2veccharcp);
    v_names.push_back(NULL); /* mark end of sequences */
    transform(alignment.begin(), alignment.end(), back_inserter(v_aln), convert_vecstring2veccharcp);
    v_aln.push_back(NULL); /* mark end of sequences */

    return vrna_file_msa_write(filename.c_str(),
                              (const char **)&v_names[0],
                              (const char **)&v_aln[0],
                              (id != "") ? id.c_str() : NULL,
                              (structure != "") ? structure.c_str() : NULL,
                              (source != "") ? source.c_str() : NULL,
                              options);
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") my_file_msa_detect_format;
%feature("kwargs") my_file_msa_detect_format;
%feature("autodoc") my_file_msa_write;
%feature("kwargs") my_file_msa_write;
%feature("autodoc") my_file_msa_read;
%feature("kwargs") my_file_msa_read;
%feature("autodoc") my_file_msa_read_record;
%feature("kwargs") my_file_msa_read_record;
#endif

unsigned int my_file_msa_detect_format( std::string   filename,
                                        unsigned int  options = VRNA_FILE_FORMAT_MSA_DEFAULT);

int my_file_msa_write(std::string               filename,
                      std::vector<std::string>  names,
                      std::vector<std::string>   alignment,
                      std::string               id = "",
                      std::string               structure = "",
                      std::string               source = "",
                      unsigned int              options = VRNA_FILE_FORMAT_MSA_STOCKHOLM | VRNA_FILE_FORMAT_MSA_APPEND);

%apply std::string              *OUTPUT { std::string *id, std::string *structure };
%apply std::vector<std::string> *OUTPUT { std::vector<std::string> *names, std::vector<std::string>  *alignment };

int my_file_msa_read( std::string               filename,
                      std::vector<std::string>  *names,
                      std::vector<std::string>  *alignment,
                      std::string               *id,
                      std::string               *structure,
                      unsigned int              options = VRNA_FILE_FORMAT_MSA_STOCKHOLM);

int my_file_msa_read_record(FILE                      *filehandle,
                            std::vector<std::string>  *names,
                            std::vector<std::string>  *alignment,
                            std::string               *id,
                            std::string               *structure,
                            unsigned int              options = VRNA_FILE_FORMAT_MSA_STOCKHOLM);

%clear std::string *id, std::string *structure;
%clear std::vector<std::string>& names, std::vector<std::string>& alignment;

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
%constant unsigned int FILE_FORMAT_MSA_APPEND    = VRNA_FILE_FORMAT_MSA_APPEND;

%include <ViennaRNA/io/file_formats_msa.h>
