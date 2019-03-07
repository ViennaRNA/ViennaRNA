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
} perl_bs_callback_t;

static perl_bs_callback_t * bind_bs_callback(SV *PerlFunc, SV *PerlData);

static void perl_wrap_bs_cb(const char *stucture, void *data);

static perl_bs_callback_t *
bind_bs_callback(SV *PerlFunc, SV *PerlData){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if((!SvOK(PerlFunc)) || (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV)){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::pbacktrack*_cb, must be code reference\n");
    return NULL;
  }

  perl_bs_callback_t *cb = (perl_bs_callback_t *)vrna_alloc(sizeof(perl_bs_callback_t));

  cb->cb = PerlFunc;      /* store callback */
  cb->data = PerlData;    /* bind data */

  return cb;
}

static void
perl_wrap_bs_cb(const char *stucture, void *data){

  SV *func;
  perl_bs_callback_t *cb = (perl_bs_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    SV *err_tmp;

    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);
    SV *structureSV = sv_newmortal();
    /* add structure and free energy to perl stack */
    sv_setpv(structureSV, stucture);
    XPUSHs(structureSV);
    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);
    PUTBACK;
    perl_call_sv(func, G_EVAL | G_DISCARD);

    err_tmp = ERRSV;
    if (SvTRUE(err_tmp)) {
      croak ("Some error occurred while executing Boltzmann sampling callback - %s\n", SvPV_nolen(err_tmp));
    }

    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  unsigned int
  pbacktrack_cb(unsigned int length,
                unsigned int num_samples,
                SV           *PerlFunc,
                SV           *PerlData = NULL)
  {
    unsigned int i;
    perl_bs_callback_t *cb = bind_bs_callback(PerlFunc, PerlData);

    if (length == 0)
      i = vrna_pbacktrack_num_cb($self,
                                 num_samples,
                                 &perl_wrap_bs_cb,
                                 (void *)cb);
    else
      i = vrna_pbacktrack5_num_cb($self,
                                  length,
                                  num_samples,
                                  &perl_wrap_bs_cb,
                                  (void *)cb);
    free(cb);

    return i;
  }

  unsigned int
  pbacktrack_nr_cb(unsigned int num_samples,
                   SV           *PerlFunc,
                   SV           *PerlData = NULL)
  {
    unsigned int i;
    perl_bs_callback_t *cb = bind_bs_callback(PerlFunc, PerlData);

    i = vrna_pbacktrack_nr_cb($self,
                              num_samples,
                              &perl_wrap_bs_cb,
                              (void *)cb);

    free(cb);

    return i;
  }

  %apply vrna_nr_memory_t *INOUT { vrna_nr_memory_t *nr_memory };

  unsigned int
  pbacktrack_nr_resume_cb(unsigned int     num_samples,
                          SV               *PerlFunc,
                          SV               *PerlData,
                          vrna_nr_memory_t *nr_memory)
  {
    unsigned int i;
    perl_bs_callback_t *cb = bind_bs_callback(PerlFunc, PerlData);

    i = vrna_pbacktrack_nr_resume_cb($self,
                                     num_samples,
                                     &perl_wrap_bs_cb,
                                     (void *)cb,
                                     nr_memory);

    free(cb);

    return i;
  }

  %clear vrna_nr_memory_t *nr_memory;
}

#endif
