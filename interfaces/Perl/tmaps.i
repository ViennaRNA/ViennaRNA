// convert between perl and C file handle
%typemap(in) FILE * {
  if (SvOK($input)) /* check for undef */
    $1 = PerlIO_findFILE(IoIFP(sv_2io($input)));
  else  $1 = NULL;
}

%typemap(typecheck) FILE * {
  if (SvOK($input))
    $1 = (sv_2io($input)) ? 1 : 0;
}


// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
        AV *tempav;
        I32 len;
        int i;
        SV  **tv;
        if (!SvROK($input))
            croak("Argument $argnum is not a reference.");
        if (SvTYPE(SvRV($input)) != SVt_PVAV)
            croak("Argument $argnum is not an array.");
        tempav = (AV*)SvRV($input);
        len = av_len(tempav);
        $1 = (char **) malloc((len+2)*sizeof(char *));
        for (i = 0; i <= len; i++) {
            tv = av_fetch(tempav, i, 0);
            $1[i] = (char *) SvPV(*tv,PL_na);
        }
        $1[i] = NULL;
};

// This tells SWIG to treat char *[], const char **, and const char *[] the same as char **
%apply char ** { char *[], const char **, const char *[] };

// Creates a new Perl array and places a NULL-terminated char ** into it
%typemap(out) char ** {
        SV **svs;
        int i = 0,len = 0;
        /* Figure out how many elements we have */
        while ($1[len])
           len++;
        svs = (SV **) malloc(len*sizeof(SV *));
        for (i = 0; i < len ; i++) {
            svs[i] = sv_newmortal();
            sv_setpv((SV*)svs[i], $1[i]);
        }
        $result = newRV_noinc((SV*)av_make(len,svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

/**
 *  Handle return parameters in argument list
 */
%typemap(argout)  std::vector<std::string> *OUTPUT {
        SV **svs;
        int i = 0;
        svs = (SV **) malloc($1->size() * sizeof(SV *));
        for (std::vector<std::string>::iterator it = $1->begin(); it != $1->end(); it++, i++) {
            svs[i] = sv_newmortal();
            sv_setpv((SV*)svs[i], it->c_str());
        }
        $result = newRV_noinc((SV*)av_make($1->size(), svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

// We simply ignore input here
%typemap(in,numinputs = 0) std::vector<std::string> *OUTPUT(std::vector<std::string> junk) {
  $1 = &junk;
}


%typemap(argout)  std::vector<unsigned int> *OUTPUT {
        SV **svs;
        int i = 0;
        svs = (SV **) malloc($1->size() * sizeof(SV *));
        for (std::vector<unsigned int>::iterator it = $1->begin(); it != $1->end(); it++, i++) {
            svs[i] = sv_newmortal();
            sv_setiv((SV*)svs[i], *it);
        }
        $result = newRV_noinc((SV*)av_make($1->size(), svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

// We simply ignore input here
%typemap(in,numinputs = 0) std::vector<unsigned int> *OUTPUT(std::vector<unsigned int> junk) {
  $1 = &junk;
}


%typemap(argout)  std::vector<int> *OUTPUT {
        SV **svs;
        int i = 0;
        svs = (SV **) malloc($1->size() * sizeof(SV *));
        for (std::vector<int>::iterator it = $1->begin(); it != $1->end(); it++, i++) {
            svs[i] = sv_newmortal();
            sv_setiv((SV*)svs[i], *it);
        }
        $result = newRV_noinc((SV*)av_make($1->size(), svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

// We simply ignore input here
%typemap(in,numinputs = 0) std::vector<int> *OUTPUT(std::vector<int> junk) {
  $1 = &junk;
}


%typemap(argout)  std::vector<float> *OUTPUT {
        SV **svs;
        int i = 0;
        svs = (SV **) malloc($1->size() * sizeof(SV *));
        for (std::vector<float>::iterator it = $1->begin(); it != $1->end(); it++, i++) {
            svs[i] = sv_newmortal();
            sv_setnv((SV*)svs[i], (double)*it);
        }
        $result = newRV_noinc((SV*)av_make($1->size(), svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

// We simply ignore input here
%typemap(in,numinputs = 0) std::vector<float> *OUTPUT(std::vector<float> junk) {
  $1 = &junk;
}


%typemap(argout)  std::vector<double> *OUTPUT {
        SV **svs;
        int i = 0;
        svs = (SV **) malloc($1->size() * sizeof(SV *));
        for (std::vector<double>::iterator it = $1->begin(); it != $1->end(); it++, i++) {
            svs[i] = sv_newmortal();
            sv_setnv((SV*)svs[i], *it);
        }
        $result = newRV_noinc((SV*)av_make($1->size(), svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

// We simply ignore input here
%typemap(in,numinputs = 0) std::vector<double> *OUTPUT(std::vector<double> junk) {
  $1 = &junk;
}


/**
 *  Handle nested vectors in function return
 */

namespace std {
  class vector;

  %typemap(out) vector<vector<double> > {
    AV *arr = newAV();
    for(unsigned int i = 0; i < $1.size(); i++) {
      AV *vec = newAV();
      for(unsigned int j = 0; j < $1[i].size(); j++) {
        SV *v = newSVnv($1[i][j]);
        if (!av_store(vec, j, v))
          SvREFCNT_dec(v);
      }
      /* store reference to array */
      av_store(arr, i, newRV_noinc((SV*) vec));
    }

    $result = newRV_noinc((SV*) arr);
    sv_2mortal( $result );
    argvi++;
  }

  %typemap(out) vector<vector<int> > {
    AV *arr = newAV();
    for(unsigned int i = 0; i < $1.size(); i++) {
      AV *vec = newAV();
      for(unsigned int j = 0; j < $1[i].size(); j++) {
        SV *v = newSViv($1[i][j]);
        if (!av_store(vec, j, v))
          SvREFCNT_dec(v);
      }
      /* store reference to array */
      av_store(arr, i, newRV_noinc((SV*) vec));
    }

    $result = newRV_noinc((SV*) arr);
    sv_2mortal( $result );
    argvi++;
  }
}


%typemap(out) int [ANY] 
{ 
  AV* av = newAV();
  int i = 0,len = 0;
  len = $1_dim0;

  for (i = 0; i < len ; i++) {
      SV* perlval = newSV(0);
      sv_setiv(perlval, (IV)$1[i]);
      av_push(av, perlval);
  }
  $result = newRV_noinc((SV*) av );
  sv_2mortal( $result );
  argvi++;
}


/*
  we need this crazy piece of argout typemap only because we don't
  want the vrna_pbacktrack_mem_t object to be appended to the results(list),
  but prepended instead. Otherwise, a simple
  %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0));
  would have sufficed already
*/
%typemap(argout) vrna_pbacktrack_mem_t *INOUT {
  /* increase output stack if necessary */
  if (argvi >= items) {
    EXTEND(sp,1);
  }
  /* move already existing return values to the back */
  for (int i = argvi; i > 0; i--) {
    ST(i) = ST(i - 1);
  }
  /* store result as first element in the stack */
  ST(0) = SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0);
  /* increase return argument counter */
  argvi++;
}


%typemap(in) vrna_pbacktrack_mem_t *INOUT (vrna_pbacktrack_mem_t *retval)
{
  if (!SvOK($input)) {
    retval = new vrna_pbacktrack_mem_t();
    $1 = retval;
  } else {
    /* INOUT in */
    SWIG_ConvertPtr($input,SWIG_as_voidptrptr(&retval), 0, SWIG_POINTER_DISOWN);
    $1 = retval;
  }
}


%typemap(in) SV *PerlFunc {
  $1 = $input;
}

%typemap(in) SV *data {
  $1 = $input;
}
