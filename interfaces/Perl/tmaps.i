// convert between perl and C file handle
%typemap(in) FILE * {
  if (SvOK($input)) /* check for undef */
    $1 = PerlIO_findFILE(IoIFP(sv_2io($input)));
  else  $1 = NULL;
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
        };
        $result = newRV_noinc((SV*)av_make(len,svs));
        sv_2mortal( $result );
        free(svs);
        argvi++;
}

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

%typemap(in) SV *PerlFunc {
  $1 = $input;
}

%typemap(in) SV *data {
  $1 = $input;
}
