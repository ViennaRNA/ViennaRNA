/**********************************************/
/* BEGIN interface for structure SHAPE constraints */
/**********************************************/

  
%extend vrna_fold_compound_t {

  int sc_add_SHAPE_deigan(std::vector<double> reactivities,
                          double m,
                          double b,
                          unsigned int options=VRNA_OPTION_DEFAULT){

    return vrna_sc_add_SHAPE_deigan($self,(const double *)&reactivities[0],m,b,options);
  }

  int sc_add_SHAPE_deigan_ali(std::vector<string> shape_files,
                              std::vector<int> shape_file_association,
                              double m,
                              double b,
                              unsigned int options=VRNA_OPTION_DEFAULT){

    std::vector<const char*>  vc;
    transform(shape_files.begin(), shape_files.end(), back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of vector */
    return vrna_sc_add_SHAPE_deigan_ali($self,(const char **) &vc[0], (const int *) &shape_file_association[0],m,b,options);
  }

  int sc_add_SHAPE_zarringhalam(std::vector<double> reactivities,
                                double b,
                                double default_value,
                                const char * shape_conversion,
                                unsigned int options=VRNA_OPTION_DEFAULT){

    return vrna_sc_add_SHAPE_zarringhalam($self,(const double *) &reactivities[0],b,default_value,shape_conversion,options);
  }
}

%include  <ViennaRNA/constraints_SHAPE.h>
