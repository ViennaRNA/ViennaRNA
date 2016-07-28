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
  SV  *cb_f;
  SV  *cb_bt;
  SV  *cb_exp_f;
  SV  *data;
  SV  *delete_data;
} perl_sc_callback_t;

static int              perl_wrap_sc_f_callback(int i, int j, int k, int l, char d, void *data);
static vrna_basepair_t  *perl_wrap_sc_bt_callback(int i, int j, int k, int l, char d, void *data);
static FLT_OR_DBL       perl_wrap_sc_exp_f_callback(int i, int j, int k, int l, char d, void *data);

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

static void
sc_add_bt_perl_callback(vrna_fold_compound_t *vc,
                        SV *PerlFunc){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if(SvTYPE(SvRV(PerlFunc)) != SVt_PVCV){
    fprintf(stderr, "Warning: invalid argument for fold_compound::sc_add_bt, must be code reference\n");
    return;
  }

  /* try to dispose of previous callback */
  perl_sc_callback_t * cb;
  vrna_sc_add_bt(vc, &perl_wrap_sc_bt_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (perl_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    if(cb->cb_bt && SvOK(cb->cb_bt))
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
    cb = (perl_sc_callback_t *)vc->sc->data;
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

static vrna_basepair_t *
perl_wrap_sc_bt_callback( int i,
                          int j,
                          int k,
                          int l,
                          char d,
                          void *data){

  int c, count, len, num_pairs;
  SV *func, *bp;
  perl_sc_callback_t *cb = (perl_sc_callback_t *)data;
  vrna_basepair_t *ptr, *pairs = NULL;
  func = cb->cb_bt;
  /* compose function call */
  if(SvOK(func)){
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
    count = perl_call_sv(func, G_ARRAY);

    SPAGAIN;

    if(count == 0){ /* no additional base pairs */
      PUTBACK;
      FREETMPS;
      LEAVE;
      return NULL;
    }

    /* when we get here, we should have got an array of something */
    len       = 10;
    num_pairs = 0;
    pairs     = (vrna_basepair_t *)vrna_alloc(sizeof(vrna_basepair_t) * len);
    /* result should be list of pairs */
    for(c=0; c < count; c++){
      /* pop SV off the stack */
      bp = POPs;
      /* if list entry is defined and a reference */
      if(SvOK(bp) && SvROK(bp)){
        /* maybe the user was so kind to create a list of vrna_basepair_t? */
        if(SWIG_ConvertPtr(bp, (void **) &ptr, SWIGTYPE_p_vrna_basepair_t, SWIG_POINTER_EXCEPTION) == 0){
          pairs[num_pairs] = *ptr;
          num_pairs++;
        } /* check whether we've got a reference to a hash */
        else if(SvTYPE(SvRV(bp)) == SVt_PVHV){
          HV *pair = (HV*)SvRV(bp);
          /* check whether the hash has 'i' and 'j' keys */
          if(hv_exists(pair, "i", 1) && hv_exists(pair, "j", 1)){
            pairs[num_pairs].i = (int)SvIV(* hv_fetch(pair, "i", 1, 0));
            pairs[num_pairs].j = (int)SvIV(* hv_fetch(pair, "j", 1, 0));
            num_pairs++;
          }
        }
        /* check whether we've got a refrence to an array */
        else if(SvTYPE(SvRV(bp)) == SVt_PVAV){
          AV *pair = (AV*)SvRV(bp);
          if(av_len(pair) == 1){ /* size of array must be 2, av_len() returns highest index */
            SV **pair_i = av_fetch(pair, 0, 0);
            SV **pair_j = av_fetch(pair, 1, 0);
            if(pair_i && pair_j){
              pairs[num_pairs].i = (int)SvIV(* pair_i);
              pairs[num_pairs].j = (int)SvIV(* pair_j);
              num_pairs++;
            }
          }
        } else {
          continue;
        }
      }
      /* increase length if necessary */
      if(num_pairs == len){
        len = (int)(1.2 * len);
        pairs = (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t)*len);
      }
    }
    /* put end marker in list */
    pairs[num_pairs].i = pairs[num_pairs].j = 0;
    pairs = (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t)*(num_pairs+1));

    PUTBACK;
    FREETMPS;
    LEAVE;
  }
  return pairs;
}

%}

static void sc_add_f_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
static void sc_add_bt_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
static void sc_add_exp_f_perl_callback(vrna_fold_compound_t *vc, SV *PerlFunc);
static void sc_add_perl_data(vrna_fold_compound_t *vc, SV *data, SV *PerlFunc);

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  void sc_add_data(SV *data, SV *PerlFunc){
    sc_add_perl_data($self, data, PerlFunc);
  }
  
  void sc_add_f(SV *PerlFunc){
    sc_add_f_perl_callback($self, PerlFunc);
  }

  void sc_add_bt(SV *PerlFunc){
    sc_add_bt_perl_callback($self, PerlFunc);
  }

  void sc_add_exp_f(SV *PerlFunc){
    sc_add_exp_f_perl_callback($self, PerlFunc);
  }
}


#endif
