/**********************************************/
/* BEGIN interface for PF window callback     */
/**********************************************/

/* see also
  http://perldoc.perl.org/perlguts.html#Working-with-SVs
*/

#ifdef SWIGPERL5

%{

typedef struct {
  SV  *cb;
  SV  *data;
} perl_pf_window_callback_t;

static perl_pf_window_callback_t * bind_pf_window_callback(SV *PerlFunc, SV *PerlData);

static void perl_wrap_pf_window_cb(FLT_OR_DBL *pr, int pr_size, int i, int max, unsigned int type, void *data);

static perl_pf_window_callback_t *
bind_pf_window_callback(SV *PerlFunc, SV *PerlData){

  /* check whether PerlFunc actually is a reference to a Perl subroutine */
  if((!SvOK(PerlFunc)) || (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV)){
    fprintf(stderr, "Warning: invalid argument 1 for fold_compound::probs_window, must be code reference\n");
    return NULL;
  }

  perl_pf_window_callback_t *cb = (perl_pf_window_callback_t *)vrna_alloc(sizeof(perl_pf_window_callback_t));

  cb->cb = PerlFunc;      /* store callback */
  cb->data = PerlData;    /* bind data */

  return cb;
}

static void
perl_wrap_pf_window_cb(FLT_OR_DBL *pr, int pr_size, int i, int max, unsigned int type, void *data){

  SV  *func, *pr_sizeSV, *iSV, *maxSV, *typeSV;
  AV  *prAV;

  perl_pf_window_callback_t *cb = (perl_pf_window_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    dSP;

    SV *err_tmp;

    /* call Perl subroutine */
    ENTER;
    SAVETMPS;
    PUSHMARK(SP);

    pr_sizeSV = sv_newmortal();
    iSV       = sv_newmortal();
    maxSV     = sv_newmortal();
    typeSV    = sv_newmortal();

    prAV = newAV();

    if (type & VRNA_PROBS_WINDOW_UP) { /* We distinguish output for unpaired probabilities */

      /* 0th element */
      av_push(prAV, newSV(0));

      /* actual values in range [1, MIN(i, max)] */
      for (int cnt = 1; cnt <= pr_size; cnt++)
        av_push(prAV, newSVnv((double)pr[cnt]));

      /* empty values in range [0,i - 1] */
      for (int cnt = pr_size + 1; cnt <= max; cnt++)
        av_push(prAV, newSV(0));

    } else { /* and pairing/stacking probabilities for pair (i, j) or ensemble free energies for subsegment [i, j] */

      /* empty values in range [0, i] */
      for (int cnt = 0; cnt <= i; cnt++)
        av_push(prAV, newSV(0));

      /* actual values in range [i + 1, pr_size] */
      for (int cnt = i + 1; cnt <= pr_size; cnt++) {
        av_push(prAV, newSVnv((double)pr[cnt]));
      }
    }

    /* add pr_size, i, max, and type to perl stack */
    sv_setiv(pr_sizeSV, (IV)pr_size);
    sv_setiv(iSV, (IV)i);
    sv_setiv(maxSV, (IV)max);
    sv_setuv(typeSV, (UV)type);

    XPUSHs(sv_2mortal(newRV_noinc((SV*) prAV)));
    XPUSHs(pr_sizeSV);
    XPUSHs(iSV);
    XPUSHs(maxSV);
    XPUSHs(typeSV);

    if(cb->data && SvOK(cb->data))          /* add data object to perl stack (if any) */
      XPUSHs(cb->data);

    PUTBACK;

    perl_call_sv(func, G_EVAL | G_DISCARD);

    SPAGAIN;

    err_tmp = ERRSV;
    if (SvTRUE(err_tmp)) {
      croak ("Some error occurred while executing sliding window partition function callback - %s\n", SvPV_nolen(err_tmp));
    }

    FREETMPS;
    LEAVE;
  }

  return /*void*/;
}

%}

/* now we bind vrna_probs_window() as method to the fold_compound object using the above callback wrapper */
%extend vrna_fold_compound_t {

  int probs_window(int ulength, unsigned int options, SV *PerlFunc, SV *PerlData = NULL) {
    perl_pf_window_callback_t *cb = bind_pf_window_callback(PerlFunc, PerlData);
    int r = vrna_probs_window($self, ulength, options, &perl_wrap_pf_window_cb, (void *)cb);
    free(cb);
    return r;
  }
}


%{

  int pfl_fold_cb(std::string sequence, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL) {
    perl_pf_window_callback_t *cb = bind_pf_window_callback(PerlFunc, PerlData);
    int r = vrna_pfl_fold_cb(sequence.c_str(), window_size, max_bp_span, &perl_wrap_pf_window_cb, (void *)cb);
    free(cb);
    return r;
  }

  int pfl_fold_up_cb(std::string sequence, int ulength, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL) {
    perl_pf_window_callback_t *cb = bind_pf_window_callback(PerlFunc, PerlData);
    int r = vrna_pfl_fold_up_cb(sequence.c_str(), ulength, window_size, max_bp_span, &perl_wrap_pf_window_cb, (void *)cb);
    free(cb);
    return r;
  }

%}

int pfl_fold_cb(std::string sequence, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL);
int pfl_fold_up_cb(std::string sequence, int ulength, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL);

#endif
