/**********************************************/
/* BEGIN interface for MFE window callbacks   */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5

%{

typedef struct {
  SV  *cb;
  SV  *data;
} perl_mfe_window_callback_t;

static perl_mfe_window_callback_t * bind_mfe_window_callback(SV *PerlFunc, SV *PerlData);

static void perl_wrap_mfe_window_cb(int start, int end, const char *stucture, float energy, void *data);

#ifdef VRNA_WITH_SVM
static void perl_wrap_mfe_window_zscore_cb(int start, int end, const char *stucture, float energy, float zscore, void *data);
#endif

static perl_mfe_window_callback_t *
bind_mfe_window_callback(SV *PerlFunc, SV *PerlData){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if((!SvOK(PerlFunc)) || (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV)){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::mfe_window_cb, must be code reference\n");
    return NULL;
  }

  perl_mfe_window_callback_t *cb = (perl_mfe_window_callback_t *)vrna_alloc(sizeof(perl_mfe_window_callback_t));

  cb->cb = PerlFunc;      /* store callback */
  cb->data = PerlData;    /* bind data */

  return cb;
}

static void
perl_wrap_mfe_window_cb(int start, int end, const char *stucture, float energy, void *data){

  SV *func;
  perl_mfe_window_callback_t *cb = (perl_mfe_window_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV *startSV     = sv_newmortal();
    SV *stopSV      = sv_newmortal();
    SV *structureSV = sv_newmortal();
    SV *energySV    = sv_newmortal();
    /* add start, end, structure, and free energy to perl stack */
    sv_setiv(startSV, (IV)start);
    sv_setiv(stopSV, (IV)end);
    sv_setpv(structureSV, stucture);
    sv_setnv(energySV, (double)energy);
    XPUSHs(startSV);
    XPUSHs(stopSV);
    XPUSHs(structureSV);
    XPUSHs(energySV);
    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_VOID);
    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}

#ifdef VRNA_WITH_SVM
static void
perl_wrap_mfe_window_zscore_cb(int start, int end, const char *stucture, float energy, float zscore, void *data){

  SV *func;
  perl_mfe_window_callback_t *cb = (perl_mfe_window_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV *startSV     = sv_newmortal();
    SV *stopSV      = sv_newmortal();
    SV *structureSV = sv_newmortal();
    SV *energySV    = sv_newmortal();
    SV *zscoreSV    = sv_newmortal();
    /* add start, end, structure, and free energy to perl stack */
    sv_setiv(startSV, (IV)start);
    sv_setiv(stopSV, (IV)end);
    sv_setpv(structureSV, stucture);
    sv_setnv(energySV, (double)energy);
    sv_setnv(zscoreSV, (double)zscore);
    XPUSHs(startSV);
    XPUSHs(stopSV);
    XPUSHs(structureSV);
    XPUSHs(energySV);
    XPUSHs(zscoreSV);
    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_VOID);
    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}
#endif

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  float mfe_window_cb(SV *PerlFunc, SV *PerlData = NULL) {
    float en;
    perl_mfe_window_callback_t *cb = bind_mfe_window_callback(PerlFunc, PerlData);
    en = vrna_mfe_window_cb($self, &perl_wrap_mfe_window_cb, (void *)cb);
    free(cb);
    return en;
  }

#ifdef VRNA_WITH_SVM
  float mfe_window_zscore_cb(double min_z, SV *PerlFunc, SV *PerlData = NULL) {
    float en;
    perl_mfe_window_callback_t *cb = bind_mfe_window_callback(PerlFunc, PerlData);
    en = vrna_mfe_window_zscore_cb($self, min_z, &perl_wrap_mfe_window_zscore_cb, (void *)cb);
    free(cb);
    return en;
  }
#endif

}

%rename (Lfold_cb) my_Lfold_cb;
%rename (Lfoldz_cb) my_Lfoldz_cb;
%rename (aliLfold_cb) my_aliLfold_cb;

%{
  float my_Lfold_cb(char *string, int window_size, SV *PerlFunc, SV *PerlData = NULL) {
    float en;
    perl_mfe_window_callback_t *cb = bind_mfe_window_callback(PerlFunc, PerlData);
    en = vrna_Lfold_cb(string, window_size, &perl_wrap_mfe_window_cb, (void *)cb);
    free(cb);
    return en;
  }

#ifdef VRNA_WITH_SVM
  float my_Lfoldz_cb(char *string, int window_size, double min_z, SV *PerlFunc, SV *PerlData = NULL) {
    float en;
    perl_mfe_window_callback_t *cb = bind_mfe_window_callback(PerlFunc, PerlData);
    en = vrna_Lfoldz_cb(string, window_size, min_z, &perl_wrap_mfe_window_zscore_cb, (void *)cb);
    free(cb);
    return en;
  }
#endif

  float my_aliLfold_cb(std::vector<std::string> alignment, int window_size, SV *PerlFunc, SV *PerlData = NULL) {
    float en;
    perl_mfe_window_callback_t *cb = bind_mfe_window_callback(PerlFunc, PerlData);
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */
    en = vrna_aliLfold_cb((const char **)&vc[0], window_size, &perl_wrap_mfe_window_cb, (void *)cb);
    free(cb);
    return en;
  }

%}

float my_Lfold_cb(char *string, int window_size, SV *PerlFunc, SV *PerlData = NULL);
%ignore vrna_Lfold_cb;

#ifdef VRNA_WITH_SVM
float my_Lfoldz_cb(char *string, int window_size, double min_z, SV *PerlFunc, SV *PerlData = NULL);
%ignore vrna_Lfoldz_cb;
#endif

float my_aliLfold_cb(std::vector<std::string> alignment, int window_size, SV *PerlFunc, SV *PerlData = NULL);
%ignore vrna_aliLfold_cb;
%ignore aliLfold_cb;


#endif
