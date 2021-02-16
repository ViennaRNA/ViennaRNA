/*
 * ********************************************
 * BEGIN interface for fold compound status
 * callback
 *********************************************
 */

/* see also
 * http://perldoc.perl.org/perlguts.html#Working-with-SVs
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

  static int
  perl_wrap_sc_f_callback(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          void          *data);


  static vrna_basepair_t *
  perl_wrap_sc_bt_callback(int            i,
                           int            j,
                           int            k,
                           int            l,
                           unsigned char  d,
                           void           *data);


  static FLT_OR_DBL
  perl_wrap_sc_exp_f_callback(int           i,
                              int           j,
                              int           k,
                              int           l,
                              unsigned char d,
                              void          *data);


  static perl_sc_callback_t *
  get_new_cb(void)
  {
    perl_sc_callback_t *cb;

    cb              = (perl_sc_callback_t *)vrna_alloc(sizeof(perl_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;

    return cb;
  }


  static void
  clean_cb_data(perl_sc_callback_t *cb)
  {
    if ((cb->data) &&
        (SvOK(cb->data))) {
      if ((cb->delete_data) &&
          (SvOK(cb->delete_data))) {
        dSP;

        SV *err_tmp;

        ENTER;
        SAVETMPS;
        PUSHMARK(SP);
        XPUSHs(cb->data);
        PUTBACK;

        perl_call_sv(cb->delete_data, G_EVAL | G_DISCARD);

        SPAGAIN;

        /* Check the eval first */
        err_tmp = ERRSV;

        if (SvTRUE(err_tmp))
          croak(
            "Some error occurred while executing generic soft constraint delete_data() callback - %s\n",
            SvPV_nolen(err_tmp));

        FREETMPS;
        LEAVE;
        SvREFCNT_dec(cb->delete_data);
      }

      SvREFCNT_dec(cb->data);

      cb->data        = NULL;
      cb->delete_data = NULL;
    }
  }


  static perl_sc_callback_t *
  reuse_or_new_cb_f(vrna_sc_t *sc)
  {
    perl_sc_callback_t *cb;

    cb = (sc->data) ? (perl_sc_callback_t *)sc->data : get_new_cb();

    /* release previous callback */
    if ((cb->cb_f) &&
        (SvOK(cb->cb_f)))
      SvREFCNT_dec(cb->cb_f);

    return cb;
  }


  static perl_sc_callback_t *
  reuse_or_new_cb_exp_f(vrna_sc_t *sc)
  {
    perl_sc_callback_t *cb;

    cb = (sc->data) ? (perl_sc_callback_t *)sc->data : get_new_cb();

    /* release previous callback */
    if ((cb->cb_exp_f) &&
        (SvOK(cb->cb_exp_f)))
      SvREFCNT_dec(cb->cb_exp_f);

    return cb;
  }


  static perl_sc_callback_t *
  reuse_or_new_cb_data(vrna_sc_t *sc)
  {
    perl_sc_callback_t *cb;

    cb = (sc->data) ? (perl_sc_callback_t *)sc->data : get_new_cb();

    clean_cb_data(cb);

    return cb;
  }


  static void
  delete_perl_sc_callback(void *data)
  {
    perl_sc_callback_t *cb = (perl_sc_callback_t *)data;

    /* first delete user data */
    clean_cb_data(cb);

    /* now dispose of the registered callbacks */
    if ((cb->cb_f) &&
        (SvOK(cb->cb_f)))
      SvREFCNT_dec(cb->cb_f);

    if ((cb->cb_bt) &&
        (SvOK(cb->cb_bt)))
      SvREFCNT_dec(cb->cb_bt);

    if ((cb->cb_exp_f) &&
        (SvOK(cb->cb_exp_f)))
      SvREFCNT_dec(cb->cb_exp_f);

    /* finally free pycallback */
    free(cb);
  }


  static int
  sc_add_f_perl_callback(vrna_fold_compound_t *vc,
                         SV                   *callback)
  {
    unsigned int        s;

    perl_sc_callback_t  *cb;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* check whether PerlFunc actually is a reference to a Perl subroutine */
        if (SvTYPE(SvRV(callback)) != SVt_PVCV) {
          fprintf(stderr,
                  "Warning: invalid argument for fold_compound::sc_add_f, must be code reference\n");
          return 0;
        } else if (vrna_sc_add_f(vc, &perl_wrap_sc_f_callback)) {
          /*
           *  The above call returns 0 on any error.
           *  Otherwise it binds the wrapper function and
           *  prepares the soft constraint data structure
           *  inside vc
           */

          /* now bind the Perl function to the wrapper structure */
          cb = reuse_or_new_cb_f(vc->sc);

          SvREFCNT_inc(callback); /* Increase referenc counter */
          cb->cb_f = callback;    /* remember callback */

          /* finaly bind callback wrapper to fold compound */
          vc->sc->data = (void *)cb;

          vc->sc->free_data = &delete_perl_sc_callback;

          return 1;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        /* Find out if the function argument actually is defined and an array reference */
        if ((SvOK(callback)) &&
            (SvROK(callback)) &&
            (SvTYPE(SvRV(callback)) == SVt_PVAV)) {
          if ((av_top_index((AV *)SvRV(callback)) + 1) < (int)vc->n_seq) {
            fprintf(stderr,
                    "Warning: invalid argument for fold_compound::sc_add_f: Too few code references in array\n");
            return 0;
          }

          if (!vc->scs)
            vrna_sc_init(vc);

          for (s = 0; s < vc->n_seq; s++) {
            /* first, dispose of previous callback if necessary */
            cb = reuse_or_new_cb_f(vc->scs[s]);

            /* try to retrieve one subroutine at a time */
            SV **f = av_fetch((AV *)SvRV(callback),
                              (SSize_t)s,
                              (I32)0);

            /* check if refrence retrieved is defined */
            if ((f) &&
                (SvOK(*f))) {
              /* warn if not a code refernce */
              if ((!SvROK(*f)) ||
                  (SvTYPE(SvRV(*f)) != SVt_PVCV)) {
                fprintf(stderr,
                        "Warning: invalid argument for fold_compound::sc_add_f, Not a code reference at %d\n",
                        s);
              } else {
                /* actually prepare and bind callback wrapper */
                SvREFCNT_inc(*f); /* Increase reference counter */
                cb->cb_f = *f;    /* remember callback */
              }
            }

            /* finaly bind callback wrapper to fold compound */
            vc->scs[s]->f         = &perl_wrap_sc_f_callback;
            vc->scs[s]->data      = (void *)cb;
            vc->scs[s]->free_data = &delete_perl_sc_callback;
          }

          return 1;
        } else {
          fprintf(stderr,
                  "Warning: invalid argument for fold_compound::sc_add_f, must be reference to array of code references\n");
        }

        break;
    }

    return 0;
  }


  static int
  sc_add_exp_f_perl_callback(vrna_fold_compound_t *vc,
                             SV                   *callback)
  {
    unsigned int        s;
    perl_sc_callback_t  *cb;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* check whether PerlFunc actually is a reference to a Perl subroutine */
        if (SvTYPE(SvRV(callback)) != SVt_PVCV) {
          fprintf(stderr,
                  "Warning: invalid argument for fold_compound::sc_add_exp_f, must be code reference\n");
          return 0;
        } else if (vrna_sc_add_exp_f(vc, &perl_wrap_sc_exp_f_callback)) {
          /*
           *  The above call returns 0 on any error.
           *  Otherwise it binds the wrapper function and
           *  prepares the soft constraint data structure
           *  inside vc
           */

          /* now bind the python function to the wrapper structure */
          cb = reuse_or_new_cb_exp_f(vc->sc);

          SvREFCNT_inc(callback);   /* Increase referenc counter */
          cb->cb_exp_f = callback;  /* remember callback */

          /* finaly bind callback wrapper to fold compound */
          vc->sc->data = (void *)cb;

          vc->sc->free_data = &delete_perl_sc_callback;

          return 1;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        /* Find out if the function argument actually is defined and an array reference */
        if ((SvOK(callback)) &&
            (SvROK(callback)) &&
            (SvTYPE(SvRV(callback)) == SVt_PVAV)) {
          if ((av_top_index((AV *)SvRV(callback)) + 1) < (int)vc->n_seq) {
            fprintf(stderr,
                    "Warning: invalid argument for fold_compound::sc_add_exp_f: Too few code references in array\n");
            return 0;
          }

          if (!vc->scs)
            vrna_sc_init(vc);

          for (s = 0; s < vc->n_seq; s++) {
            /* first, dispose of previous callback if necessary */
            cb = reuse_or_new_cb_exp_f(vc->scs[s]);

            /* try to retrieve one subroutine at a time */
            SV **f = av_fetch((AV *)SvRV(callback),
                              (SSize_t)s,
                              (I32)0);

            /* check if refrence retrieved is defined */
            if ((f) &&
                (SvOK(*f))) {
              /* warn if not a code refernce */
              if ((!SvROK(*f)) ||
                  (SvTYPE(SvRV(*f)) != SVt_PVCV)) {
                fprintf(stderr,
                        "Warning: invalid argument for fold_compound::sc_add_exp_f, Not a code reference at %d\n",
                        s);
              } else {
                /* actually prepare and bind callback wrapper */
                SvREFCNT_inc(*f);   /* Increase reference counter */
                cb->cb_exp_f = *f;  /* remember callback */
              }
            }

            /* finaly bind callback wrapper to fold compound */
            vc->scs[s]->exp_f     = &perl_wrap_sc_exp_f_callback;
            vc->scs[s]->data      = (void *)cb;
            vc->scs[s]->free_data = &delete_perl_sc_callback;
          }

          return 1;
        } else {
          fprintf(stderr,
                  "Warning: invalid argument for fold_compound::sc_add_exp_f, must be reference to array of code references\n");
        }

        break;
    }

    return 0;
  }


  static int
  sc_add_bt_perl_callback(vrna_fold_compound_t  *vc,
                          SV                    *PerlFunc)
  {
    /* check whether PerlFunc actually is a reference to a Perl subroutine */
    if (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV) {
      fprintf(stderr,
              "Warning: invalid argument for fold_compound::sc_add_bt, must be code reference\n");
      return 0;
    }

    /* try to dispose of previous callback */
    perl_sc_callback_t *cb;

    if (vrna_sc_add_bt(vc, &perl_wrap_sc_bt_callback)) {
      /* now bind the Perl function to the wrapper structure */
      if (vc->sc->data) {
        cb = (perl_sc_callback_t *)vc->sc->data;
        /* release previous callback */
        if ((cb->cb_bt) &&
            (SvOK(cb->cb_bt)))
          SvREFCNT_dec(cb->cb_bt);
      } else {
        cb = get_new_cb();
      }

      cb->cb_bt = PerlFunc;   /* remember callback */
      SvREFCNT_inc(PerlFunc); /* Increase referenc counter */

      /* finaly bind callback wrapper to fold compound */
      vc->sc->data      = (void *)cb;
      vc->sc->free_data = &delete_perl_sc_callback;

      return 1;
    }

    return 0;
  }


  static int
  sc_add_perl_data(vrna_fold_compound_t *vc,
                   SV                   *data,
                   SV                   *PerlFunc)
  {
    unsigned int        s;
    perl_sc_callback_t  *cb;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* check whether data is defined */
        if (!SvOK(data)) {
          fprintf(stderr, "Warning: argument 1 for fold_compound::sc_add_data, must be defined\n");
          return 0;
        }

        /* check whether free-callback is defined and callable */
        if ((PerlFunc) &&
            (SvOK(PerlFunc))) {
          if (SvTYPE(SvRV(PerlFunc)) != SVt_PVCV) {
            fprintf(stderr,
                    "Warning: argument 2 for fold_compound::sc_add_data, must be undef or code reference\n");
            return 0;
          }
        }

        /* create soft constraints data structure */
        if (!vc->sc)
          vrna_sc_init(vc);

        /* try to dispose of previous data */
        cb = reuse_or_new_cb_data(vc->sc);

        cb->data        = data;     /* remember data */
        cb->delete_data = PerlFunc; /* remember delete data function */
        /* increase reference counter */
        if ((data) &&
            (SvOK(data)))
          SvREFCNT_inc(data);

        if ((PerlFunc) &&
            (SvOK(PerlFunc)))
          SvREFCNT_inc(PerlFunc);

        vc->sc->data      = (void *)cb;
        vc->sc->free_data = &delete_perl_sc_callback;

        return 1;

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        int data_good, free_data_good;

        data_good = free_data_good = 0;

        /* Find out if the data argument actually is defined and an array reference */
        if ((SvOK(data)) &&
            (SvROK(data)) &&
            (SvTYPE(SvRV(data)) == SVt_PVAV)) {
          if ((av_top_index((AV *)SvRV(data)) + 1) < (int)vc->n_seq) {
            fprintf(stderr,
                    "Warning: invalid argument for fold_compound::sc_add_data: Too few code references in array\n");
            return 0;
          }

          data_good = 1;
        }

        if ((SvOK(PerlFunc)) &&
            (SvROK(PerlFunc)) &&
            (SvTYPE(SvRV(PerlFunc)) == SVt_PVAV))
          if ((av_top_index((AV *)SvRV(PerlFunc)) + 1) >= (int)vc->n_seq)
            free_data_good = 1;

        if (data_good) {
          if (!vc->scs)
            vrna_sc_init(vc);

          for (s = 0; s < vc->n_seq; s++) {
            /* first, dispose of previous callback if necessary */
            cb = reuse_or_new_cb_data(vc->scs[s]);

            /* try to retrieve one subroutine at a time */
            SV **d = av_fetch((AV *)SvRV(data),
                              (SSize_t)s,
                              (I32)0);

            /* check if refrence retrieved is defined */
            if ((d) &&
                (SvOK(*d))) {
              /* actually prepare and bind callback wrapper */
              SvREFCNT_inc(*d); /* Increase reference counter */
              cb->data = *d;    /* remember data */
            }

            if (free_data_good) {
              SV **f = av_fetch((AV *)SvRV(PerlFunc),
                                (SSize_t)s,
                                (I32)0);

              if ((f) &&
                  (SvOK(*f))) {
                /* warn if not a code refernce */
                if ((!SvROK(*f)) ||
                    (SvTYPE(SvRV(*f)) != SVt_PVCV)) {
                  fprintf(stderr,
                          "Warning: invalid argument for fold_compound::sc_add_data, Not a code reference at %d\n",
                          s);
                } else {
                  /* actually prepare and bind callback wrapper */
                  SvREFCNT_inc(*f);     /* Increase reference counter */
                  cb->delete_data = *f; /* remember callback */
                }
              }
            }

            /* finaly bind callback wrapper to fold compound */
            vc->scs[s]->data      = (void *)cb;
            vc->scs[s]->free_data = &delete_perl_sc_callback;
          }

          return 1;
        }

        break;
    }

    return 0;
  }


  static int
  perl_wrap_sc_f_callback(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          void          *data)
  {
    int                 ret, count;
    SV                  *func;

    perl_sc_callback_t  *cb = (perl_sc_callback_t *)data;

    func = cb->cb_f;

    ret = 0;

    if (func && SvOK(func)) {
      dSP;

      SV  *err_tmp;
      I32 ax; /* expected by ST(x) macro */

      /* call Perl subroutine */
      ENTER;
      SAVETMPS;
      PUSHMARK(SP);

      SV  *pSVi = sv_newmortal();
      sv_setiv(pSVi, (IV)i);
      XPUSHs(pSVi);
      SV  *pSVj = sv_newmortal();
      sv_setiv(pSVj, (IV)j);
      XPUSHs(pSVj);
      SV  *pSVk = sv_newmortal();
      sv_setiv(pSVk, (IV)k);
      XPUSHs(pSVk);
      SV  *pSVl = sv_newmortal();
      sv_setiv(pSVl, (IV)l);
      XPUSHs(pSVl);
      SV  *pSVd = sv_newmortal();
      sv_setiv(pSVd, (IV)d);
      XPUSHs(pSVd);

      if (cb->data && SvOK(cb->data))        /* add data object to perl stack (if any) */
        XPUSHs(cb->data);

      PUTBACK;

      count = perl_call_sv(func, G_EVAL | G_ARRAY);

      SPAGAIN; /* refresh local copy of the perl stack pointer */

      /* prepare the stack such that we can use the ST(x) macro */
      SP  -= count;
      ax  = (SP - PL_stack_base) + 1;

      /* Check the eval first */
      err_tmp = ERRSV;
      if (SvTRUE(err_tmp)) {
        croak("Some error occurred while executing generic soft constraint callback - %s\n",
              SvPV_nolen(err_tmp));
        POPs;
      } else {
        /* we expect a single value to be returned */
        if (count != 1)         /* more or less than 1 return value */
          croak("Generic soft constraint callback must return exactly 1 item\n");
        else if (!SvOK(ST(0)))  /* return value is undefined */
          croak("Generic soft constraint callback must not return undef\n");
        else if (!SvIOK(ST(0))) /* return value is not an integer scalar */
          croak("Generic soft constraint callback must return pseudo energy value\n");
        else
          ret = SvIV(ST(0));
      }

      PUTBACK;
      FREETMPS;
      LEAVE;
    }

    return ret;
  }


  static FLT_OR_DBL
  perl_wrap_sc_exp_f_callback(int           i,
                              int           j,
                              int           k,
                              int           l,
                              unsigned char d,
                              void          *data)
  {
    int                 count;
    FLT_OR_DBL          ret;
    SV                  *func;

    perl_sc_callback_t  *cb = (perl_sc_callback_t *)data;

    func  = cb->cb_f;
    ret   = 1.0;

    if (SvOK(func)) {
      dSP;

      SV  *err_tmp;
      I32 ax; /* expected by ST(x) macro */

      /* call Perl subroutine */
      ENTER;
      SAVETMPS;
      PUSHMARK(SP);

      SV  *pSVi = sv_newmortal();
      sv_setiv(pSVi, (IV)i);
      XPUSHs(pSVi);
      SV  *pSVj = sv_newmortal();
      sv_setiv(pSVj, (IV)j);
      XPUSHs(pSVj);
      SV  *pSVk = sv_newmortal();
      sv_setiv(pSVk, (IV)k);
      XPUSHs(pSVk);
      SV  *pSVl = sv_newmortal();
      sv_setiv(pSVl, (IV)l);
      XPUSHs(pSVl);
      SV  *pSVd = sv_newmortal();
      sv_setiv(pSVd, (IV)d);
      XPUSHs(pSVd);

      if (cb->data && SvOK(cb->data))        /* add data object to perl stack (if any) */
        XPUSHs(cb->data);

      PUTBACK;

      count = perl_call_sv(func, G_EVAL | G_ARRAY);

      SPAGAIN; /* refresh local copy of the perl stack pointer */

      /* prepare the stack such that we can use the ST(x) macro */
      SP  -= count;
      ax  = (SP - PL_stack_base) + 1;

      /* Check the eval first */
      err_tmp = ERRSV;
      if (SvTRUE(err_tmp)) {
        croak(
          "Some error occurred while executing generic soft constraint callback (partition function) - %s\n",
          SvPV_nolen(err_tmp));
        POPs;
      } else {
        /* we expect a single value to be returned */
        if (count != 1)         /* more or less than 1 return value */
          croak("Generic soft constraint callback (partition function) must return exactly 1 item\n");
        else if (!SvOK(ST(0)))  /* return value is undefined */
          croak("Generic soft constraint callback (partition function) must not return undef\n");
        else if (!SvIOK(ST(0))) /* return value is not an integer scalar */
          croak(
            "Generic soft constraint callback (partition function) must return Boltzmann weighted pseudo energy value\n");
        else
          ret = SvIV(ST(0));
      }

      PUTBACK;
      FREETMPS;
      LEAVE;
    }

    return ret;
  }


  static vrna_basepair_t *
  perl_wrap_sc_bt_callback(int            i,
                           int            j,
                           int            k,
                           int            l,
                           unsigned char  d,
                           void           *data)
  {
    int                 c, count, len, num_pairs;
    SV                  *func, *bp;
    perl_sc_callback_t  *cb = (perl_sc_callback_t *)data;
    vrna_basepair_t     *ptr, *pairs = NULL;

    func = cb->cb_bt;

    /* compose function call */
    if (SvOK(func)) {
      dSP;

      SV *err_tmp;

      ENTER;
      SAVETMPS;
      PUSHMARK(SP);

      SV  *pSVi = sv_newmortal();
      sv_setiv(pSVi, (IV)i);
      XPUSHs(pSVi);
      SV  *pSVj = sv_newmortal();
      sv_setiv(pSVj, (IV)j);
      XPUSHs(pSVj);
      SV  *pSVk = sv_newmortal();
      sv_setiv(pSVk, (IV)k);
      XPUSHs(pSVk);
      SV  *pSVl = sv_newmortal();
      sv_setiv(pSVl, (IV)l);
      XPUSHs(pSVl);
      SV  *pSVd = sv_newmortal();
      sv_setiv(pSVd, (IV)d);
      XPUSHs(pSVd);

      if (cb->data && SvOK(cb->data))        /* add data object to perl stack (if any) */
        XPUSHs(cb->data);

      PUTBACK;

      count = perl_call_sv(func, G_EVAL | G_ARRAY);

      SPAGAIN;

      /* Check the eval first */
      err_tmp = ERRSV;
      if (SvTRUE(err_tmp)) {
        croak(
          "Some error occurred while executing generic soft constraint backtrack callback - %s\n",
          SvPV_nolen(err_tmp));
        POPs;
        PUTBACK;
        FREETMPS;
        LEAVE;
        return NULL;
      } else if (count == 0) {
        /* no additional base pairs */
        PUTBACK;
        FREETMPS;
        LEAVE;
        return NULL;
      } else {
        /* when we get here, we should have got an array of something */
        len       = 10;
        num_pairs = 0;
        pairs     = (vrna_basepair_t *)vrna_alloc(sizeof(vrna_basepair_t) * len);
        /* result should be list of pairs */
        for (c = 0; c < count; c++) {
          /* pop SV off the stack */
          bp = POPs;
          /* if list entry is defined and a reference */
          if (SvOK(bp) && SvROK(bp)) {
            /* maybe the user was so kind to create a list of vrna_basepair_t? */
            if (SWIG_ConvertPtr(bp, (void **)&ptr, SWIGTYPE_p_vrna_basepair_t,
                                SWIG_POINTER_EXCEPTION) == 0) {
              pairs[num_pairs] = *ptr;
              num_pairs++;
            } else if (SvTYPE(SvRV(bp)) == SVt_PVHV) {
              /* check whether we've got a reference to a hash */
              HV *pair = (HV *)SvRV(bp);
              /* check whether the hash has 'i' and 'j' keys */
              if (hv_exists(pair, "i", 1) && hv_exists(pair, "j", 1)) {
                pairs[num_pairs].i  = (int)SvIV(*hv_fetch(pair, "i", 1, 0));
                pairs[num_pairs].j  = (int)SvIV(*hv_fetch(pair, "j", 1, 0));
                num_pairs++;
              }
            } else if (SvTYPE(SvRV(bp)) == SVt_PVAV) {
              /* check whether we've got a refrence to an array */
              AV *pair = (AV *)SvRV(bp);
              if (av_len(pair) == 1) {
                /* size of array must be 2, av_len() returns highest index */
                SV  **pair_i  = av_fetch(pair, 0, 0);
                SV  **pair_j  = av_fetch(pair, 1, 0);
                if (pair_i && pair_j) {
                  pairs[num_pairs].i  = (int)SvIV(*pair_i);
                  pairs[num_pairs].j  = (int)SvIV(*pair_j);
                  num_pairs++;
                }
              }
            } else {
              continue;
            }
          }

          /* increase length if necessary */
          if (num_pairs == len) {
            len   = (int)(1.2 * len);
            pairs = (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t) * len);
          }
        }
      }

      /* put end marker in list */
      pairs[num_pairs].i  = pairs[num_pairs].j = 0;
      pairs               =
        (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t) * (num_pairs + 1));

      PUTBACK;
      FREETMPS;
      LEAVE;
    }

    return pairs;
  }


%}


static int
sc_add_f_perl_callback(vrna_fold_compound_t *vc,
                       SV                   *PerlFunc);


static int
sc_add_bt_perl_callback(vrna_fold_compound_t  *vc,
                        SV                    *PerlFunc);


static int
sc_add_exp_f_perl_callback(vrna_fold_compound_t *vc,
                           SV                   *PerlFunc);


static int
sc_add_perl_data(vrna_fold_compound_t *vc,
                 SV                   *data,
                 SV                   *PerlFunc);


/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {
  int
  sc_add_data(SV  *data,
              SV  *PerlFunc)
  {
    return sc_add_perl_data($self, data, PerlFunc);
  }


  int
  sc_add_f(SV *PerlFunc)
  {
    return sc_add_f_perl_callback($self, PerlFunc);
  }


  int
  sc_add_bt(SV *PerlFunc)
  {
    return sc_add_bt_perl_callback($self, PerlFunc);
  }


  int
  sc_add_exp_f(SV *PerlFunc)
  {
    return sc_add_exp_f_perl_callback($self, PerlFunc);
  }
}


#endif
