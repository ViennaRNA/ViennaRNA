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

  int array_len;
  SV  *func, **pr_listSVS, *pr_sizeSV, *iSV, *maxSV, *typeSV, *prRV;
  AV  *prAV;

  perl_pf_window_callback_t *cb = (perl_pf_window_callback_t *)data;

  func  = cb->cb;

  if(func && SvOK(func)){
    /* call Perl subroutine */
    dSP;
    ENTER;
    SAVETMPS;
    PUSHMARK(sp);

    pr_sizeSV = sv_newmortal();
    iSV       = sv_newmortal();
    maxSV     = sv_newmortal();
    typeSV    = sv_newmortal();

    if (type & VRNA_PROBS_WINDOW_UP) { /* We distinguish output for unpaired probabilities */

      array_len = max + 1;

      /* create Perl array for unpaired probabilities */
      pr_listSVS = (SV **)vrna_alloc(sizeof(SV*) * array_len);

      /* 0th element */
      pr_listSVS[0] = newSV(0);

      /* actual values in range [1, MIN(i, max)] */
      for (int cnt = 1; cnt <= pr_size; cnt++) {
        pr_listSVS[cnt] = sv_newmortal();
        sv_setnv(pr_listSVS[cnt], (double)pr[cnt]);
      }

      /* empty values in range [0,i - 1] */
      for (int cnt = pr_size + 1; cnt <= max; cnt++)
        pr_listSVS[cnt] = newSV(0);

    } else { /* and pairing/stacking probabilities for pair (i, j) or ensemble free energies for subsegment [i, j] */

      array_len = pr_size + 1;

      /* create Perl array for pr values */
      pr_listSVS = (SV **)vrna_alloc(sizeof(SV*) * array_len);

      /* empty values in range [0, i] */
      for (int cnt = 0; cnt <= i; cnt++)
        pr_listSVS[cnt] = newSV(0);

      /* actual values in range [i + 1, pr_size] */
      for (int cnt = i + 1; cnt <= pr_size; cnt++) {
        pr_listSVS[cnt] = sv_newmortal();
        sv_setnv(pr_listSVS[cnt], (double)pr[cnt]);
      }
    }

    prAV = av_make(array_len, pr_listSVS);

    free(pr_listSVS);

    prRV = newRV_inc((SV*) prAV);

    /* add pr_size, i, max, and type to perl stack */
    sv_setiv(pr_sizeSV, (IV)pr_size);
    sv_setiv(iSV, (IV)i);
    sv_setiv(maxSV, (IV)max);
    sv_setuv(typeSV, (UV)type);

    XPUSHs(prRV);
    XPUSHs(pr_sizeSV);
    XPUSHs(iSV);
    XPUSHs(maxSV);
    XPUSHs(typeSV);

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

/* now we bind vrna_probs_window() as method to the fold_compound object using the above callback wrapper */
%extend vrna_fold_compound_t {

  void probs_window(int ulength, unsigned int options, SV *PerlFunc, SV *PerlData = NULL) {
    perl_pf_window_callback_t *cb = bind_pf_window_callback(PerlFunc, PerlData);
    vrna_probs_window($self, ulength, options, &perl_wrap_pf_window_cb, (void *)cb);
    free(cb);
  }
}


%{

  void pfl_fold_cb(std::string sequence, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL) {
    perl_pf_window_callback_t *cb = bind_pf_window_callback(PerlFunc, PerlData);
    vrna_pfl_fold_cb(sequence.c_str(), window_size, max_bp_span, &perl_wrap_pf_window_cb, (void *)cb);
    free(cb);
  }

  void pfl_fold_up_cb(std::string sequence, int ulength, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL) {
    perl_pf_window_callback_t *cb = bind_pf_window_callback(PerlFunc, PerlData);
    vrna_pfl_fold_up_cb(sequence.c_str(), ulength, window_size, max_bp_span, &perl_wrap_pf_window_cb, (void *)cb);
    free(cb);
  }

%}

void pfl_fold_cb(std::string sequence, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL);
void pfl_fold_up_cb(std::string sequence, int ulength, int window_size, int max_bp_span, SV *PerlFunc, SV *PerlData = NULL);

#endif
