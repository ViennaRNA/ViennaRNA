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
  vrna_md_defaults_dangles((int)PyLong_AsLong($input));
}

%typemap(varout) int dangles {
  $result = PyLong_FromLong((long)vrna_md_defaults_dangles_get());
}

%typemap(varin) int tetra_loop {
  vrna_md_defaults_special_hp((int)PyLong_AsLong($input));
}

%typemap(varout) int tetra_loop {
  $result = PyLong_FromLong((long)vrna_md_defaults_special_hp_get());
}

%typemap(varin) int special_hp {
  vrna_md_defaults_special_hp((int)PyLong_AsLong($input));
}

%typemap(varout) int special_hp {
  $result = PyLong_FromLong((long)vrna_md_defaults_special_hp_get());
}

%typemap(varin) int noLonelyPairs {
  vrna_md_defaults_noLP((int)PyLong_AsLong($input));
}

%typemap(varout) int noLonelyPairs {
  $result = PyLong_FromLong((long)vrna_md_defaults_noLP_get());
}

%typemap(varin) int noLP {
  vrna_md_defaults_noLP((int)PyLong_AsLong($input));
}

%typemap(varout) int noLP {
  $result = PyLong_FromLong((long)vrna_md_defaults_noLP_get());
}

%typemap(varin) int noGU {
  vrna_md_defaults_noGU((int)PyLong_AsLong($input));
}

%typemap(varout) int noGU {
  $result = PyLong_FromLong((long)vrna_md_defaults_noGU_get());
}

%typemap(varin) int no_closingGU {
  vrna_md_defaults_noGUclosure((int)PyLong_AsLong($input));
}

%typemap(varout) int no_closingGU {
  $result = PyLong_FromLong((long)vrna_md_defaults_noGUclosure_get());
}

%typemap(varin) int noGUclosure {
  vrna_md_defaults_noGUclosure((int)PyLong_AsLong($input));
}

%typemap(varout) int noGUclosure {
  $result = PyLong_FromLong((long)vrna_md_defaults_noGUclosure_get());
}

%typemap(varin) int logML {
  vrna_md_defaults_logML((int)PyLong_AsLong($input));
}

%typemap(varout) int logML {
  $result = PyLong_FromLong((long)vrna_md_defaults_logML_get());
}

%typemap(varin) int circ {
  vrna_md_defaults_circ((int)PyLong_AsLong($input));
}

%typemap(varout) int circ {
  $result = PyLong_FromLong((long)vrna_md_defaults_circ_get());
}

%typemap(varin) int gquad {
  vrna_md_defaults_gquad((int)PyLong_AsLong($input));
}

%typemap(varout) int gquad {
  $result = PyLong_FromLong((long)vrna_md_defaults_gquad_get());
}

%typemap(varin) int uniq_ML {
  vrna_md_defaults_uniq_ML((int)PyLong_AsLong($input));
}

%typemap(varout) int uniq_ML {
  $result = PyLong_FromLong((long)vrna_md_defaults_uniq_ML_get());
}

%typemap(varin) int energy_set {
  vrna_md_defaults_energy_set((int)PyLong_AsLong($input));
}

%typemap(varout) int energy_set {
  $result = PyLong_FromLong((long)vrna_md_defaults_energy_set_get());
}

%typemap(varin) int backtrack {
  vrna_md_defaults_backtrack((int)PyLong_AsLong($input));
}

%typemap(varout) int backtrack {
  $result = PyLong_FromLong((long)vrna_md_defaults_backtrack_get());
}

/* backtrack_type implementation still missing */

%typemap(varin) int do_backtrack {
  vrna_md_defaults_compute_bpp((int)PyLong_AsLong($input));
}

%typemap(varout) int do_backtrack {
  $result = PyLong_FromLong((long)vrna_md_defaults_compute_bpp_get());
}

%typemap(varin) int compute_bpp {
  vrna_md_defaults_compute_bpp((int)PyLong_AsLong($input));
}

%typemap(varout) int compute_bpp {
  $result = PyLong_FromLong((long)vrna_md_defaults_compute_bpp_get());
}

%typemap(varin) int max_bp_span {
  vrna_md_defaults_max_bp_span((int)PyLong_AsLong($input));
}

%typemap(varout) int max_bp_span {
  $result = PyLong_FromLong((long)vrna_md_defaults_max_bp_span_get());
}

%typemap(varin) int min_loop_size {
  vrna_md_defaults_min_loop_size((int)PyLong_AsLong($input));
}

%typemap(varout) int min_loop_size {
  $result = PyLong_FromLong((long)vrna_md_defaults_min_loop_size_get());
}

%typemap(varin) int window_size {
  vrna_md_defaults_window_size((int)PyLong_AsLong($input));
}

%typemap(varout) int window_size {
  $result = PyLong_FromLong((long)vrna_md_defaults_window_size_get());
}

%typemap(varin) int olAliEn {
  vrna_md_defaults_olAliEn((int)PyLong_AsLong($input));
}

%typemap(varout) int olAliEn {
  $result = PyLong_FromLong((long)vrna_md_defaults_olAliEn_get());
}

%typemap(varin) int ribo {
  vrna_md_defaults_ribo((int)PyLong_AsLong($input));
}

%typemap(varout) int ribo {
  $result = PyLong_FromLong((long)vrna_md_defaults_ribo_get());
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

