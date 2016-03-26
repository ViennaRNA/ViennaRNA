/**********************************************/
/* BEGIN interface for fold compound status   */
/* callback                                   */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5
%{

typedef struct {
  SV  *cb;
  SV  *data;
  SV  *delete_data;
} perl_callback_t;

typedef struct {
  SV  *cb_f;
  SV  *cb_bt;
  SV  *cb_exp_f;
  SV  *data;
  SV  *delete_data;
} perl_sc_callback_t;

static void       perl_wrap_fc_status_callback(unsigned char status, void *data);
static int        perl_wrap_sc_f_callback(int i, int j, int k, int l, char d, void *data);
static FLT_OR_DBL perl_wrap_sc_exp_f_callback(int i, int j, int k, int l, char d, void *data);

static void
delete_perl_callback(void * data){

  perl_callback_t *cb = (perl_callback_t *)data;
  /* first delete user data */
  if(cb->data && SvOK(cb->data)){
    if(cb->delete_data && SvOK(cb->delete_data)){
      /* call Perl subroutine */
      dSP;
      ENTER;
      SAVETMPS;
      PUSHMARK(sp);
      XPUSHs(cb->data);
      PUTBACK;
      perl_call_sv(cb->delete_data, G_VOID);
      FREETMPS;
      LEAVE;
      SvREFCNT_dec(cb->delete_data);
    }
    SvREFCNT_dec(cb->data);
  }

  /* now dispose of the callback */
  if(cb->cb && SvOK(cb->cb))
    SvREFCNT_dec(cb->cb);

  /* finally free perl callback */
  free(cb);
}

static void
delete_perl_sc_callback(void * data){

  perl_sc_callback_t *cb = (perl_sc_callback_t *)data;
  /* first delete user data */
  if(cb->data && SvOK(cb->data)){
    if(cb->delete_data && SvOK(cb->delete_data)){
      dSP;
      ENTER;
      SAVETMPS;
      PUSHMARK(sp);
      XPUSHs(cb->data);
      PUTBACK;
      perl_call_sv(cb->delete_data, G_VOID);
      FREETMPS;
      LEAVE;
      SvREFCNT_dec(cb->delete_data);
    }
    SvREFCNT_dec(cb->data);
  }

  /* now dispose of the registered callbacks */
  if(cb->cb_f && SvOK(cb->cb_f))
    SvREFCNT_dec(cb->cb_f);
  if(cb->cb_bt && SvOK(cb->cb_bt))
    SvREFCNT_dec(cb->cb_bt);
  if(cb->cb_exp_f && SvOK(cb->cb_exp_f))
    SvREFCNT_dec(cb->cb_exp_f);

  /* finally free pycallback */
  free(cb);
}


static void
fc_add_perl_data(vrna_fold_compound_t *vc,
                SV *data,
                SV *PerlFunc){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(!SvOK(data)){
    fprintf(stderr, "Warning: argument 1 for fold_compound::add_auxdata, must be defined\n");
    return;
  }

  if(PerlFunc && SvOK(PerlFunc)){
    if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
      fprintf(stderr, "Warning: argument 2 for fold_compound::add_auxdata, must be undef or code reference\n");
      return;
    }
  }

  perl_callback_t * cb;
  /* try to dispose of previous data */
  if(vc->auxdata){
    cb = (perl_callback_t *)vc->auxdata;
    if(cb->data && SvOK(cb->data)){
      if(cb->delete_data && SvOK(cb->delete_data)){
        /* call Perl subroutine */
        dSP;
        ENTER;
        SAVETMPS;
        PUSHMARK(sp);
        XPUSHs(cb->data);
        PUTBACK;
        perl_call_sv(cb->delete_data, G_VOID);
        FREETMPS;
        LEAVE;
        SvREFCNT_dec(cb->delete_data);
      }
      SvREFCNT_dec(cb->data);
    }
  } else {
    cb              = (perl_callback_t *)vrna_alloc(sizeof(perl_callback_t));
    cb->cb          = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->data        = data;     /* remember data */
  cb->delete_data = PerlFunc; /* remember delete data function */
  /* increase reference counter */
  if(data && SvOK(data))
    SvREFCNT_inc(data);
  if(PerlFunc && SvOK(PerlFunc))
    SvREFCNT_inc(PerlFunc);

  vc->auxdata = (void *)cb;
  if(!vc->free_auxdata){
    vc->free_auxdata = &delete_perl_callback;
  }
}

static void
fc_add_perl_callback( vrna_fold_compound_t *vc,
                      SV *PerlFunc){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument for fold_compound::add_callback, must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl_callback_t * cb;
  if(vc->auxdata){
    cb = (perl_callback_t *)vc->auxdata;
    /* release previous callback */
    SvREFCNT_dec(cb->cb);
  } else {
    cb = (perl_callback_t *)vrna_alloc(sizeof(perl_callback_t));
    cb->data = NULL;
    cb->delete_data = NULL;
  }
  cb->cb = PerlFunc;      /* remember callback */
  SvREFCNT_inc(PerlFunc); /* Increase reference counter */

  /* finaly bind callback wrapper to fold compound */
  vc->auxdata = (void *)cb;
  if(!vc->free_auxdata){
    vc->free_auxdata = &delete_perl_callback;
  }
  vrna_fold_compound_add_callback(vc, &perl_wrap_fc_status_callback);
}

static void
sc_add_f_perl_callback( vrna_fold_compound_t *vc,
                        SV *PerlFunc){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument for fold_compound::sc_add_f, must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl_sc_callback_t * cb;
  vrna_sc_add_f(vc, &perl_wrap_sc_f_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (perl_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    if(cb->cb_f && SvOK(cb->cb_f))
      SvREFCNT_dec(cb->cb_f);
  } else {
    cb = (perl_sc_callback_t *)vrna_alloc(sizeof(perl_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->cb_f = PerlFunc;      /* remember callback */
  SvREFCNT_inc(PerlFunc);   /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_perl_sc_callback;
  }
}

static void
sc_add_exp_f_perl_callback( vrna_fold_compound_t *vc,
                            SV *PerlFunc){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument for fold_compound::sc_add_exp_f, must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl_sc_callback_t * cb;
  vrna_sc_add_exp_f(vc, &perl_wrap_sc_exp_f_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (perl_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    if(cb->cb_exp_f && SvOK(cb->cb_exp_f))
      SvREFCNT_dec(cb->cb_exp_f);
  } else {
    cb = (perl_sc_callback_t *)vrna_alloc(sizeof(perl_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->cb_exp_f = PerlFunc;  /* remember callback */
  SvREFCNT_inc(PerlFunc);   /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_perl_sc_callback;
  }
}

#if 0
static vrna_basepair_t *
sc_add_bt_perl_callback(vrna_fold_compound_t *vc,
                        SV *PerlFunc){

  /* try to dispose of previous callback */
  perl_sc_callback_t * cb;
  vrna_sc_add_bt(vc, &perl_wrap_sc_bt_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (perl_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    if(SvOK(cb->cb_bt))
      SvREFCNT_dec(cb->cb_bt);
  } else {
    cb = (perl_sc_callback_t *)vrna_alloc(sizeof(perl_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->cb_bt = PerlFunc;   /* remember callback */
  SvREFCNT_inc(PerlFunc); /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_perl_sc_callback;
  }
}
#endif

static void
sc_add_perl_data( vrna_fold_compound_t *vc,
                  SV *data,
                  SV *PerlFunc){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(!SvOK(data)){
    fprintf(stderr, "Warning: argument 1 for fold_compound::sc_add_data, must be defined\n");
    return;
  }

  if(PerlFunc && SvOK(PerlFunc)){
    if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
      fprintf(stderr, "Warning: argument 2 for fold_compound::sc_add_data, must be undef or code reference\n");
      return;
    }
  }

  perl_sc_callback_t * cb;

  /* create soft constraints data structure */
  if(!vc->sc)
    vrna_sc_init(vc);

  /* try to dispose of previous data */
  if(vc->sc->data){
    cb = (perl_sc_callback_t *)vc->auxdata;
    if(cb->data && SvOK(cb->data)){
      if(cb->delete_data && SvOK(cb->delete_data)){
        dSP;
        ENTER;
        SAVETMPS;
        PUSHMARK(sp);
        XPUSHs(cb->data);
        PUTBACK;
        perl_call_sv(cb->delete_data, G_VOID);
        FREETMPS;
        LEAVE;
        SvREFCNT_dec(cb->delete_data);
      }
      SvREFCNT_dec(cb->data);
    }
  } else {
    cb              = (perl_sc_callback_t *)vrna_alloc(sizeof(perl_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->data        = data;     /* remember data */
  cb->delete_data = PerlFunc; /* remember delete data function */
  /* increase reference counter */
  if(data && SvOK(data))
    SvREFCNT_inc(data);
  if(PerlFunc && SvOK(PerlFunc))
    SvREFCNT_inc(PerlFunc);

  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_perl_sc_callback;
  }
}

static void
perl_wrap_fc_status_callback( unsigned char status,
                              void *data){

  SV *func;
  perl_callback_t *cb = (perl_callback_t *)data;
  func = cb->cb;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV* pSV = sv_newmortal();
    sv_setiv(pSV, (IV)status);  /* add status value to perl stack */
    XPUSHs(pSV);
    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_VOID);
    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}

static int
perl_wrap_sc_f_callback(int i,
                        int j,
                        int k,
                        int l,
                        char d,
                        void *data){

  int ret, count;
  SV *func;
  perl_sc_callback_t *cb = (perl_sc_callback_t *)data;

  func = cb->cb_f;

  ret = 0;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV* pSVi = sv_newmortal();
    sv_setiv(pSVi, (IV)i);
    XPUSHs(pSVi);
    SV* pSVj = sv_newmortal();
    sv_setiv(pSVj, (IV)j);
    XPUSHs(pSVj);
    SV* pSVk = sv_newmortal();
    sv_setiv(pSVk, (IV)k);
    XPUSHs(pSVk);
    SV* pSVl = sv_newmortal();
    sv_setiv(pSVl, (IV)l);
    XPUSHs(pSVl);
    SV* pSVd = sv_newmortal();
    sv_setiv(pSVd, (IV)d);
    XPUSHs(pSVd);

    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    count = perl_call_sv(func, G_SCALAR);

    SPAGAIN; /* refresh local copy of the perl stack pointer */

    if (count != 1)
      croak("sc_f_callback function did not return a single value\n");

    ret = POPi;

    FREETMPS;
    LEAVE;
  }

  return ret;
}

static FLT_OR_DBL
perl_wrap_sc_exp_f_callback(int i,
                            int j,
                            int k,
                            int l,
                            char d,
                            void *data){

  int count;
  FLT_OR_DBL ret;
  SV *func;

  perl_sc_callback_t *cb = (perl_sc_callback_t *)data;

  func = cb->cb_f;
  ret = 1.0;

  if(SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV* pSVi = sv_newmortal();
    sv_setiv(pSVi, (IV)i);
    XPUSHs(pSVi);
    SV* pSVj = sv_newmortal();
    sv_setiv(pSVj, (IV)j);
    XPUSHs(pSVj);
    SV* pSVk = sv_newmortal();
    sv_setiv(pSVk, (IV)k);
    XPUSHs(pSVk);
    SV* pSVl = sv_newmortal();
    sv_setiv(pSVl, (IV)l);
    XPUSHs(pSVl);
    SV* pSVd = sv_newmortal();
    sv_setiv(pSVd, (IV)d);
    XPUSHs(pSVd);

    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    count = perl_call_sv(func, G_SCALAR);

    SPAGAIN; /* refresh local copy of the perl stack pointer */

    if (count != 1)
      croak("sc_exp_f_callback function did not return a single value\n");

    ret = (FLT_OR_DBL)POPn;

    FREETMPS;
    LEAVE;
  }

  return ret;
}

%}

static void fc_add_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
static void fc_add_perl_data(vrna_fold_compound_t *vc, SV *data, SV *PerlFunc);
static void sc_add_f_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
//static void sc_add_bt_perl_callback(vrna_fold_compound_t *vc, SV *PyFunc);
static void sc_add_exp_f_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
static void sc_add_perl_data(vrna_fold_compound_t *vc, SV *data, SV *PerlFunc);

%typemap(in) SV* {
  $1 = $input;
}
#endif
