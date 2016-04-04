/**********************************************/
/* BEGIN interface for structure soft constraints */
/**********************************************/



%extend vrna_fold_compound_t {


  void sc_remove()
  {
	  vrna_sc_remove($self);
  }
  void sc_init()
  {
	  vrna_sc_init($self);
  }
  
  
  /*not abel to use double vector multidimensional, not recognized
  void sc_add_bp(std::vector<std::vector<int>> constraints,unsigned int options=VRNA_OPTION_MFE)
  {
	  std::cout <<"geht";
  }*/
  


 void sc_add_up(std::vector<FLT_OR_DBL> constraint,unsigned int options=VRNA_OPTION_MFE)
 {
	vrna_sc_add_up($self, (const FLT_OR_DBL *)&constraint[0], options);
 }

  
  #ifdef SWIGPYTHON
  void add_auxdata(PyObject *data, PyObject *free_data){
    fc_add_pydata($self, data, free_data);
  }

  void add_callback(PyObject *PyFunc){
    fc_add_pycallback($self, PyFunc);
  }

  void sc_add_data(PyObject *data, PyObject *free_data){
    sc_add_pydata($self, data, free_data);
  }
  
  void sc_add_f(PyObject *PyFunc){
    sc_add_f_pycallback($self, PyFunc);
  }

  void sc_add_exp_f(PyObject *PyFunc){
    sc_add_exp_f_pycallback($self, PyFunc);
  }

#endif

#ifdef SWIGPERL5
  void add_auxdata(SV *data, SV *free_data){
    fc_add_perl_data($self, data, free_data);
  }

  void add_callback(SV *PerlFunc){
    fc_add_perl_callback($self, PerlFunc);
  }

  void sc_add_data(SV *data, SV *free_data){
    sc_add_perl_data($self, data, free_data);
  }
  
  void sc_add_f(SV *PerlFunc){
    sc_add_f_perl_callback($self, PerlFunc);
  }

  void sc_add_exp_f(SV *PerlFunc){
    sc_add_exp_f_perl_callback($self, PerlFunc);
  }

#endif
  
}



%include  <ViennaRNA/constraints_soft.h>
