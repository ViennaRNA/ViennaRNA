/*
 *  get / set global model detail variables
 */

#ifdef SWIGPERL5
%typemap(varin) double temperature {
  vrna_md_defaults_temperature((double)SvNV($input));
}

%typemap(varout) double temperature {
  sv_setnv($result, (double) vrna_md_defaults_temperature_get());
}

%typemap(varin) double betaScale {
  vrna_md_defaults_betaScale((double)SvNV($input));
}

%typemap(varout) double betaScale {
  sv_setnv($result, (double) vrna_md_defaults_betaScale_get());
}

%typemap(varin) int dangles {
  vrna_md_defaults_dangles(SvIV($input));
}

%typemap(varout) int dangles {
  sv_setiv($result, (IV) vrna_md_defaults_dangles_get());
}

%typemap(varin) int tetra_loop {
  vrna_md_defaults_special_hp(SvIV($input));
}

%typemap(varout) int tetra_loop {
  sv_setiv($result, (IV) vrna_md_defaults_special_hp_get());
}

%typemap(varin) int special_hp {
  vrna_md_defaults_special_hp(SvIV($input));
}

%typemap(varout) int special_hp {
  sv_setiv($result, (IV) vrna_md_defaults_special_hp_get());
}

%typemap(varin) int noLonelyPairs {
  vrna_md_defaults_noLP(SvIV($input));
}

%typemap(varout) int noLonelyPairs {
  sv_setiv($result, (IV) vrna_md_defaults_noLP_get());
}

%typemap(varin) int noLP {
  vrna_md_defaults_noLP(SvIV($input));
}

%typemap(varout) int noLP {
  sv_setiv($result, (IV) vrna_md_defaults_noLP_get());
}

%typemap(varin) int noGU {
  vrna_md_defaults_noGU(SvIV($input));
}

%typemap(varout) int noGU {
  sv_setiv($result, (IV) vrna_md_defaults_noGU_get());
}

%typemap(varin) int no_closingGU {
  vrna_md_defaults_noGUclosure(SvIV($input));
}

%typemap(varout) int no_closingGU {
  sv_setiv($result, (IV) vrna_md_defaults_noGUclosure_get());
}

%typemap(varin) int noGUclosure {
  vrna_md_defaults_noGUclosure(SvIV($input));
}

%typemap(varout) int noGUclosure {
  sv_setiv($result, (IV) vrna_md_defaults_noGUclosure_get());
}

%typemap(varin) int logML {
  vrna_md_defaults_logML(SvIV($input));
}

%typemap(varout) int logML {
  sv_setiv($result, (IV) vrna_md_defaults_logML_get());
}

%typemap(varin) int circ {
  vrna_md_defaults_circ(SvIV($input));
}

%typemap(varout) int circ {
  sv_setiv($result, (IV) vrna_md_defaults_circ_get());
}

%typemap(varin) int gquad {
  vrna_md_defaults_gquad(SvIV($input));
}

%typemap(varout) int gquad {
  sv_setiv($result, (IV) vrna_md_defaults_gquad_get());
}

%typemap(varin) int uniq_ML {
  vrna_md_defaults_uniq_ML(SvIV($input));
}

%typemap(varout) int uniq_ML {
  sv_setiv($result, (IV) vrna_md_defaults_uniq_ML_get());
}

%typemap(varin) int energy_set {
  vrna_md_defaults_energy_set(SvIV($input));
}

%typemap(varout) int energy_set {
  sv_setiv($result, (IV) vrna_md_defaults_energy_set_get());
}

%typemap(varin) int backtrack {
  vrna_md_defaults_backtrack(SvIV($input));
}

%typemap(varout) int backtrack {
  sv_setiv($result, (IV) vrna_md_defaults_backtrack_get());
}

/* backtrack_type implementation still missing */

%typemap(varin) int do_backtrack {
  vrna_md_defaults_compute_bpp(SvIV($input));
}

%typemap(varout) int do_backtrack {
  sv_setiv($result, (IV) vrna_md_defaults_compute_bpp_get());
}

%typemap(varin) int compute_bpp {
  vrna_md_defaults_compute_bpp(SvIV($input));
}

%typemap(varout) int compute_bpp {
  sv_setiv($result, (IV) vrna_md_defaults_compute_bpp_get());
}

%typemap(varin) int max_bp_span {
  vrna_md_defaults_max_bp_span(SvIV($input));
}

%typemap(varout) int max_bp_span {
  sv_setiv($result, (IV) vrna_md_defaults_max_bp_span_get());
}

%typemap(varin) int min_loop_size {
  vrna_md_defaults_min_loop_size(SvIV($input));
}

%typemap(varout) int min_loop_size {
  sv_setiv($result, (IV) vrna_md_defaults_min_loop_size_get());
}

%typemap(varin) int window_size {
  vrna_md_defaults_window_size(SvIV($input));
}

%typemap(varout) int window_size {
  sv_setiv($result, (IV) vrna_md_defaults_window_size_get());
}

%typemap(varin) int oldAliEn {
  vrna_md_defaults_oldAliEn(SvIV($input));
}

%typemap(varout) int oldAliEn {
  sv_setiv($result, (IV) vrna_md_defaults_oldAliEn_get());
}

%typemap(varin) int ribo {
  vrna_md_defaults_ribo(SvIV($input));
}

%typemap(varout) int ribo {
  sv_setiv($result, (IV) vrna_md_defaults_ribo_get());
}

%typemap(varin) double cv_fact {
  vrna_md_defaults_cv_fact((double)SvNV($input));
}

%typemap(varout) double cv_fact {
  sv_setnv($result, (double) vrna_md_defaults_cv_fact_get());
}

%typemap(varin) double nc_fact {
  vrna_md_defaults_nc_fact((double)SvNV($input));
}

%typemap(varout) double nc_fact {
  sv_setnv($result, (double) vrna_md_defaults_nc_fact_get());
}

%typemap(varin) double sfact {
  vrna_md_defaults_sfact((double)SvNV($input));
}

%typemap(varout) double sfact {
  sv_setnv($result, (double) vrna_md_defaults_sfact_get());
}

#endif
