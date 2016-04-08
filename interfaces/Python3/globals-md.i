/*
 *  get / set global model detail variables
 */

// Python typemaps
#ifdef SWIGPYTHON
%typemap(varin) double temperature {
  vrna_md_defaults_temperature(PyFloat_AsDouble($input));
}

%typemap(varout) double temperature {
  $result = PyFloat_FromDouble(vrna_md_defaults_temperature_get());
}

%typemap(varin) double betaScale {
  vrna_md_defaults_betaScale(PyFloat_AsDouble($input));
}

%typemap(varout) double betaScale {
  $result = PyFloat_FromDouble(vrna_md_defaults_betaScale_get());
}

%typemap(varin) int dangles {
  vrna_md_defaults_dangles(PyInt_AsLong($input));
}

%typemap(varout) int dangles {
  $result = PyInt_FromLong(vrna_md_defaults_dangles_get());
}

%typemap(varin) int tetra_loop {
  vrna_md_defaults_special_hp(PyInt_AsLong($input));
}

%typemap(varout) int tetra_loop {
  $result = PyInt_FromLong(vrna_md_defaults_special_hp_get());
}

%typemap(varin) int special_hp {
  vrna_md_defaults_special_hp(PyInt_AsLong($input));
}

%typemap(varout) int special_hp {
  $result = PyInt_FromLong(vrna_md_defaults_special_hp_get());
}

%typemap(varin) int noLonelyPairs {
  vrna_md_defaults_noLP(PyInt_AsLong($input));
}

%typemap(varout) int noLonelyPairs {
  $result = PyInt_FromLong(vrna_md_defaults_noLP_get());
}

%typemap(varin) int noLP {
  vrna_md_defaults_noLP(PyInt_AsLong($input));
}

%typemap(varout) int noLP {
  $result = PyInt_FromLong(vrna_md_defaults_noLP_get());
}

%typemap(varin) int noGU {
  vrna_md_defaults_noGU(PyInt_AsLong($input));
}

%typemap(varout) int noGU {
  $result = PyInt_FromLong(vrna_md_defaults_noGU_get());
}

%typemap(varin) int no_closingGU {
  vrna_md_defaults_noGUclosure(PyInt_AsLong($input));
}

%typemap(varout) int no_closingGU {
  $result = PyInt_FromLong(vrna_md_defaults_noGUclosure_get());
}

%typemap(varin) int noGUclosure {
  vrna_md_defaults_noGUclosure(PyInt_AsLong($input));
}

%typemap(varout) int noGUclosure {
  $result = PyInt_FromLong(vrna_md_defaults_noGUclosure_get());
}

%typemap(varin) int logML {
  vrna_md_defaults_logML(PyInt_AsLong($input));
}

%typemap(varout) int logML {
  $result = PyInt_FromLong(vrna_md_defaults_logML_get());
}

%typemap(varin) int circ {
  vrna_md_defaults_circ(PyInt_AsLong($input));
}

%typemap(varout) int circ {
  $result = PyInt_FromLong(vrna_md_defaults_circ_get());
}

%typemap(varin) int gquad {
  vrna_md_defaults_gquad(PyInt_AsLong($input));
}

%typemap(varout) int gquad {
  $result = PyInt_FromLong(vrna_md_defaults_gquad_get());
}

%typemap(varin) int uniq_ML {
  vrna_md_defaults_uniq_ML(PyInt_AsLong($input));
}

%typemap(varout) int uniq_ML {
  $result = PyInt_FromLong(vrna_md_defaults_uniq_ML_get());
}

%typemap(varin) int energy_set {
  vrna_md_defaults_energy_set(PyInt_AsLong($input));
}

%typemap(varout) int energy_set {
  $result = PyInt_FromLong(vrna_md_defaults_energy_set_get());
}

%typemap(varin) int backtrack {
  vrna_md_defaults_backtrack(PyInt_AsLong($input));
}

%typemap(varout) int backtrack {
  $result = PyInt_FromLong(vrna_md_defaults_backtrack_get());
}

/* backtrack_type implementation still missing */

%typemap(varin) int do_backtrack {
  vrna_md_defaults_compute_bpp(PyInt_AsLong($input));
}

%typemap(varout) int do_backtrack {
  $result = PyInt_FromLong(vrna_md_defaults_compute_bpp_get());
}

%typemap(varin) int compute_bpp {
  vrna_md_defaults_compute_bpp(PyInt_AsLong($input));
}

%typemap(varout) int compute_bpp {
  $result = PyInt_FromLong(vrna_md_defaults_compute_bpp_get());
}

%typemap(varin) int max_bp_span {
  vrna_md_defaults_max_bp_span(PyInt_AsLong($input));
}

%typemap(varout) int max_bp_span {
  $result = PyInt_FromLong(vrna_md_defaults_max_bp_span_get());
}

%typemap(varin) int min_loop_size {
  vrna_md_defaults_min_loop_size(PyInt_AsLong($input));
}

%typemap(varout) int min_loop_size {
  $result = PyInt_FromLong(vrna_md_defaults_min_loop_size_get());
}

%typemap(varin) int window_size {
  vrna_md_defaults_window_size(PyInt_AsLong($input));
}

%typemap(varout) int window_size {
  $result = PyInt_FromLong(vrna_md_defaults_window_size_get());
}

%typemap(varin) int olAliEn {
  vrna_md_defaults_olAliEn(PyInt_AsLong($input));
}

%typemap(varout) int olAliEn {
  $result = PyInt_FromLong(vrna_md_defaults_olAliEn_get());
}

%typemap(varin) int ribo {
  vrna_md_defaults_ribo(PyInt_AsLong($input));
}

%typemap(varout) int ribo {
  $result = PyInt_FromLong(vrna_md_defaults_ribo_get());
}

%typemap(varin) double cv_fact {
  vrna_md_defaults_cv_fact(PyFloat_AsDouble($input));
}

%typemap(varout) double cv_fact {
  $result = PyFloat_FromDouble(vrna_md_defaults_cv_fact_get());
}

%typemap(varin) double nc_fact {
  vrna_md_defaults_nc_fact(PyFloat_AsDouble($input));
}

%typemap(varout) double nc_fact {
  $result = PyFloat_FromDouble(vrna_md_defaults_nc_fact_get());
}

%typemap(varin) double sfact {
  vrna_md_defaults_sfact(PyFloat_AsDouble($input));
}

%typemap(varout) double sfact {
  $result = PyFloat_FromDouble(vrna_md_defaults_sfact_get());
}

#endif

