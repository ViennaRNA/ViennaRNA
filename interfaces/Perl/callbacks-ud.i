/**********************************************/
/* BEGIN interface for unstructured domains   */
/* feature callbacks                          */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5
%{

typedef struct {
  SV  *prod_rule;
  SV  *exp_prod_rule;
  SV  *energy;
  SV  *exp_energy;
  SV  *data;
  SV  *delete_data;
  SV  *prob_add;
  SV  *prob_get;
} perl5_ud_callback_t;

static vrna_callback_ud_production      perl5_wrap_ud_prod_rule;
static vrna_callback_ud_exp_production  perl5_wrap_ud_exp_prod_rule;
static vrna_callback_ud_energy          perl5_wrap_ud_energy;
static vrna_callback_ud_exp_energy      perl5_wrap_ud_exp_energy;
static vrna_callback_ud_probs_add       perl5_wrap_ud_prob_add;
static vrna_callback_ud_probs_get       perl5_wrap_ud_prob_get;

static perl5_ud_callback_t *
new_perl_ud_cb(void){

  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)vrna_alloc(sizeof(perl5_ud_callback_t));

  cb->prod_rule     = NULL;
  cb->exp_prod_rule = NULL;
  cb->energy        = NULL;
  cb->exp_energy    = NULL;
  cb->data          = NULL;
  cb->delete_data   = NULL;
  cb->prob_add      = NULL;
  cb->prob_get      = NULL;

  return cb;
}



static void
delete_perl_ud_callback(void * data){

  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)data;
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
  if(cb->prod_rule && SvOK(cb->prod_rule))
    SvREFCNT_dec(cb->prod_rule);
  if(cb->exp_prod_rule && SvOK(cb->exp_prod_rule))
    SvREFCNT_dec(cb->exp_prod_rule);
  if(cb->energy && SvOK(cb->energy))
    SvREFCNT_dec(cb->energy);
  if(cb->exp_energy && SvOK(cb->exp_energy))
    SvREFCNT_dec(cb->exp_energy);
  if(cb->prob_add && SvOK(cb->prob_add))
    SvREFCNT_dec(cb->prob_add);
  if(cb->prob_get && SvOK(cb->prob_get))
    SvREFCNT_dec(cb->prob_get);

  /* finally free callback */
  free(cb);
}


static void
ud_set_data(vrna_fold_compound_t *vc,
            SV *data,
            SV *PerlFunc){

  if(!SvOK(data)){
    fprintf(stderr, "Warning: argument 1 for fold_compound::ud_set_data must be defined\n");
    return;
  }

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(PerlFunc && SvOK(PerlFunc)){
    if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
      fprintf(stderr, "Warning: argument 2 for fold_compound::ud_set_data must be undef or code reference\n");
      return;
    }
  }

  perl5_ud_callback_t * cb;
  /* try to dispose of previous data */
  if(vc->domains_up && vc->domains_up->data){
    cb = (perl5_ud_callback_t *)vc->domains_up->data;
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
    cb  = new_perl_ud_cb();
  }
  cb->data        = data;     /* remember data */
  cb->delete_data = PerlFunc; /* remember delete data function */
  /* increase reference counter */
  if(data && SvOK(data))
    SvREFCNT_inc(data);
  if(PerlFunc && SvOK(PerlFunc))
    SvREFCNT_inc(PerlFunc);

  /* bind callback wrapper to fold compound */
  vrna_ud_set_data(vc, (void *)cb, &delete_perl_ud_callback);
}

static void
ud_set_prod_rule_cb(vrna_fold_compound_t *vc,
                    SV *prod_cb,
                    SV *eval_cb){

  /* check whether prod_cb and eval_cb are references to Perl subroutines */
  if(SvTYPE(SvRV(prod_cb)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::ud_set_prod_rule_cb must be code reference\n");
    return;
  }
  if(SvTYPE(SvRV(eval_cb)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument 2 for fold_compound::ud_set_prod_rule_cb must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl5_ud_callback_t * cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (perl5_ud_callback_t *)vc->domains_up->data;
    /* release previous callbacks */
    if(cb->prod_rule && SvOK(cb->prod_rule))
      SvREFCNT_dec(cb->prod_rule);
    if(cb->energy && SvOK(cb->energy))
      SvREFCNT_dec(cb->energy);
  } else {
    cb = new_perl_ud_cb();
    vrna_ud_set_data(vc, (void *)cb, &delete_perl_ud_callback);
  }
  cb->prod_rule = prod_cb;      /* remember callback */
  cb->energy    = eval_cb;      /* remember callback */
  SvREFCNT_inc(prod_cb);   /* Increase reference counter */
  SvREFCNT_inc(eval_cb);   /* Increase reference counter */

  vrna_ud_set_prod_rule_cb(vc, &perl5_wrap_ud_prod_rule, &perl5_wrap_ud_energy);
}

static void
ud_set_exp_prod_rule_cb(vrna_fold_compound_t *vc,
                        SV *prod_cb,
                        SV *eval_cb){

  /* check whether prod_cb and eval_cb are references to Perl subroutines */
  if(SvTYPE(SvRV(prod_cb)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::ud_set_exp_prod_rule_cb must be code reference\n");
    return;
  }
  if(SvTYPE(SvRV(eval_cb)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument 2 for fold_compound::ud_set_exp_prod_rule_cb must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl5_ud_callback_t * cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (perl5_ud_callback_t *)vc->domains_up->data;
    /* release previous callbacks */
    if(cb->exp_prod_rule && SvOK(cb->exp_prod_rule))
      SvREFCNT_dec(cb->exp_prod_rule);
    if(cb->exp_energy && SvOK(cb->exp_energy))
      SvREFCNT_dec(cb->exp_energy);
  } else {
    cb = new_perl_ud_cb();
    vrna_ud_set_data(vc, (void *)cb, &delete_perl_ud_callback);
  }
  cb->exp_prod_rule = prod_cb;      /* remember callback */
  cb->exp_energy    = eval_cb;      /* remember callback */
  SvREFCNT_inc(prod_cb);   /* Increase reference counter */
  SvREFCNT_inc(eval_cb);   /* Increase reference counter */

  vrna_ud_set_exp_prod_rule_cb(vc, &perl5_wrap_ud_exp_prod_rule, &perl5_wrap_ud_exp_energy);
}

static void
ud_set_prob_cb( vrna_fold_compound_t *vc,
                SV *setter,
                SV *getter){

  /* check whether prod_cb and eval_cb are references to Perl subroutines */
  if(SvTYPE(SvRV(setter)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::ud_set_prob_cb must be code reference\n");
    return;
  }
  if(SvTYPE(SvRV(getter)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument 2 for fold_compound::ud_set_prob_cb must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl5_ud_callback_t * cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (perl5_ud_callback_t *)vc->domains_up->data;
    /* release previous callbacks */
    if(cb->prob_add && SvOK(cb->prob_add))
      SvREFCNT_dec(cb->prob_add);
    if(cb->prob_get && SvOK(cb->prob_get))
      SvREFCNT_dec(cb->prob_get);
  } else {
    cb = new_perl_ud_cb();
    vrna_ud_set_data(vc, (void *)cb, &delete_perl_ud_callback);
  }
  cb->prob_add = setter;  /* remember callback */
  cb->prob_get = getter;  /* remember callback */
  SvREFCNT_inc(setter);   /* Increase reference counter */
  SvREFCNT_inc(getter);   /* Increase reference counter */

  vrna_ud_set_prob_cb(vc, &perl5_wrap_ud_prob_add, &perl5_wrap_ud_prob_get);
}


static void
perl5_wrap_ud_prod_rule(vrna_fold_compound_t *vc,
                        void *data){

  SV *func;
  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)data;
  func = cb->prod_rule;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    if(cb->data && SvOK(cb->data))  /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_VOID);
    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}


static void
perl5_wrap_ud_exp_prod_rule(vrna_fold_compound_t *vc,
                            void *data){

  SV *func;
  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)data;
  func = cb->exp_prod_rule;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    if(cb->data && SvOK(cb->data))  /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_VOID);
    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}


static int
perl5_wrap_ud_energy( vrna_fold_compound_t *vc,
                      int i,
                      int j,
                      unsigned int looptype,
                      void *data){

  int ret, count;
  SV *func;
  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)data;

  func = cb->energy;

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
    SV* pSVt = sv_newmortal();
    sv_setiv(pSVt, (IV)looptype);
    XPUSHs(pSVt);

    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    count = perl_call_sv(func, G_SCALAR);

    SPAGAIN; /* refresh local copy of the perl stack pointer */

    if (count != 1)
      croak("ud_energy_callback function did not return a single value\n");

    ret = POPi;

    FREETMPS;
    LEAVE;
  }

  return ret;
}


static FLT_OR_DBL
perl5_wrap_ud_exp_energy( vrna_fold_compound_t *vc,
                          int i,
                          int j,
                          unsigned int looptype,
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
    SV* pSVt = sv_newmortal();
    sv_setiv(pSVt, (IV)looptype);
    XPUSHs(pSVt);

    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    count = perl_call_sv(func, G_SCALAR);

    SPAGAIN; /* refresh local copy of the perl stack pointer */

    if (count != 1)
      croak("ud_exp_energy_callback function did not return a single value\n");

    ret = (FLT_OR_DBL)POPn;

    FREETMPS;
    LEAVE;
  }

  return ret;
}


static void
perl5_wrap_ud_prob_add( vrna_fold_compound_t *vc,
                        int i,
                        int j,
                        unsigned int looptype,
                        FLT_OR_DBL prob,
                        void *data){

  SV *func;
  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)data;
  func = cb->prob_add;

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
    SV* pSVt = sv_newmortal();
    sv_setiv(pSVt, (IV)looptype);
    XPUSHs(pSVt);
    SV* pSVp = sv_newmortal();
    sv_setnv(pSVp, (double)prob);
    XPUSHs(pSVp);

    if(cb->data && SvOK(cb->data))  /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_VOID);
    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}


static FLT_OR_DBL
perl5_wrap_ud_prob_get( vrna_fold_compound_t *vc,
                        int i,
                        int j,
                        unsigned int looptype,
                        int motif,
                        void *data){

  int count;
  FLT_OR_DBL ret;
  SV *func;
  perl5_ud_callback_t *cb = (perl5_ud_callback_t *)data;

  func  = cb->prob_add;
  ret   = 0.;

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
    SV* pSVt = sv_newmortal();
    sv_setiv(pSVt, (IV)looptype);
    XPUSHs(pSVt);
    SV* pSVm = sv_newmortal();
    sv_setiv(pSVm, (IV)motif);
    XPUSHs(pSVt);

    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    count = perl_call_sv(func, G_SCALAR);

    SPAGAIN; /* refresh local copy of the perl stack pointer */

    if (count != 1)
      croak("ud_probs_get_callback function did not return a single value\n");

    ret = (FLT_OR_DBL)POPn;

    FREETMPS;
    LEAVE;
  }

  return ret;
}

%}

static void ud_set_data(vrna_fold_compound_t *vc, SV *data, SV *PerlFunc);
static void ud_set_prod_rule_cb(vrna_fold_compound_t *vc, SV *prod_cb, SV *eval_cb);
static void ud_set_exp_prod_rule_cb(vrna_fold_compound_t *vc, SV *prod_cb, SV *eval_cb);
static void ud_set_prob_cb(vrna_fold_compound_t *vc, SV *setter, SV *getter);

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  void ud_set_data(SV *data, SV *PerlFunc){
    ud_set_data($self, data, PerlFunc);
  }

  void ud_set_prod_rule_cb(SV *prod_cb, SV *eval_cb){

    ud_set_prod_rule_cb($self, prod_cb, eval_cb);
  }

  void ud_set_exp_prod_rule_cb(SV *prod_cb, SV *eval_cb){

    ud_set_exp_prod_rule_cb($self, prod_cb, eval_cb);
  }

  void ud_set_prob_cb(SV *setter, SV *getter){

    ud_set_prob_cb($self, setter, getter);
  }
}


#endif
