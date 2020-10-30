/**********************************************/
/* BEGIN interface for heat_capacity callbacks       */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5

%{

typedef struct {
  SV  *cb;
  SV  *data;
} perl_heat_capacity_callback_t;

static perl_heat_capacity_callback_t * bind_heat_capacity_callback(SV *PerlFunc, SV *PerlData);

static void perl_wrap_heat_capacity_cb(float temperature, float hc, void *data);

static perl_heat_capacity_callback_t *
bind_heat_capacity_callback(SV *PerlFunc, SV *PerlData){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if((!SvOK(PerlFunc)) || (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV)){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::heat_capacity_cb, must be code reference\n");
    return NULL;
  }

  perl_heat_capacity_callback_t *cb = (perl_heat_capacity_callback_t *)vrna_alloc(sizeof(perl_heat_capacity_callback_t));

  cb->cb = PerlFunc;      /* store callback */
  cb->data = PerlData;    /* bind data */

  return cb;
}

static void
perl_wrap_heat_capacity_cb(float temp, float hc, void *data){

  SV *func;
  perl_heat_capacity_callback_t *cb = (perl_heat_capacity_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    dSP;

    SV *err_tmp;

    /* call Perl subroutine */
    ENTER;
    SAVETMPS;
    PUSHMARK(SP);

    SV *tempSV = sv_newmortal();
    SV *hcSV   = sv_newmortal();
    /* add structure and free energy to perl stack */
    sv_setnv(tempSV, (double)temp);
    sv_setnv(hcSV, (double)hc);
    XPUSHs(tempSV);
    XPUSHs(hcSV);
    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;

    perl_call_sv(func, G_EVAL | G_DISCARD);

    SPAGAIN;

    err_tmp = ERRSV;
    if (SvTRUE(err_tmp)) {
      croak ("Some error occurred while executing heat_capacity callback - %s\n", SvPV_nolen(err_tmp));
    }

    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  void
  heat_capacity_cb(float        T_min,
                   float        T_max,
                   float        T_increment,
                   unsigned int mpoints,
                   SV           *PerlFunc,
                   SV           *PerlData = NULL)
  {
    perl_heat_capacity_callback_t *cb = bind_heat_capacity_callback(PerlFunc, PerlData);
    vrna_heat_capacity_cb($self, T_min, T_max, T_increment, mpoints, &perl_wrap_heat_capacity_cb, (void *)cb);
    free(cb);
  }

}


#endif
