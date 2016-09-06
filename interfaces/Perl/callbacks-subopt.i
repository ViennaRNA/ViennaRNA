/**********************************************/
/* BEGIN interface for subopt callbacks       */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5

%{

typedef struct {
  SV  *cb;
  SV  *data;
} perl_subopt_callback_t;

static perl_subopt_callback_t * bind_subopt_callback(SV *PerlFunc, SV *PerlData);

static void perl_wrap_subopt_cb(const char *stucture, float energy, void *data);

static perl_subopt_callback_t *
bind_subopt_callback(SV *PerlFunc, SV *PerlData){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if((!SvOK(PerlFunc)) || (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV)){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::subopt_cb, must be code reference\n");
    return NULL;
  }

  perl_subopt_callback_t *cb = (perl_subopt_callback_t *)vrna_alloc(sizeof(perl_subopt_callback_t));

  cb->cb = PerlFunc;      /* store callback */
  cb->data = PerlData;    /* bind data */

  return cb;
}

static void
perl_wrap_subopt_cb(const char *stucture, float energy, void *data){

  SV *func;
  perl_subopt_callback_t *cb = (perl_subopt_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV *structureSV = sv_newmortal();
    SV *energySV    = sv_newmortal();
    /* add structure and free energy to perl stack */
    sv_setpv(structureSV, stucture);
    sv_setnv(energySV, (double)energy);
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

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  void subopt_cb(int delta, SV *PerlFunc, SV *PerlData = NULL){

    perl_subopt_callback_t *cb = bind_subopt_callback(PerlFunc, PerlData);
    vrna_subopt_cb($self, delta, &perl_wrap_subopt_cb, (void *)cb);
    free(cb);
  }

}


#endif
