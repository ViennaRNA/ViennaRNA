/**********************************************/
/* BEGIN interface for fold compound status   */
/* callback                                   */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5
static void fc_add_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
static void fc_add_perl_data(vrna_fold_compound_t *vc, SV *data, SV *PerlFunc);

%{

typedef struct {
  SV  *cb;
  SV  *data;
  SV  *delete_data;
} perl_callback_t;

static void perl_wrap_fc_status_callback(unsigned char status, void *data);

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

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {
  void add_auxdata(SV *data, SV *PerlFunc){
    fc_add_perl_data($self, data, PerlFunc);
  }

  void add_callback(SV *PerlFunc){
    fc_add_perl_callback($self, PerlFunc);
  }
}


#endif
