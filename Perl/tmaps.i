// Perl typemaps
 #ifdef SWIGPERL5
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


// Creates a new Perl array and places a NULL-terminated char ** into it
%typemap(out) char ** {
        AV *myav;
        SV **svs;
        int i = 0,len = 0;
        /* Figure out how many elements we have */
        while ($1[len])
           len++;
        svs = (SV **) malloc(len*sizeof(SV *));
        for (i = 0; i < len ; i++) {
            svs[i] = sv_newmortal();
            sv_setpv((SV*)svs[i],$1[i]);
        };
        myav =  av_make(len,svs);
        free(svs);
        $result = newRV((SV*)myav);
        sv_2mortal($result);
        argvi++;
}

#endif

// Typemaps that are independent of scripting language

// This cleans up the char ** array after the function call
%typemap(freearg) char ** {
         free($1);
}
// Now a few test functions
//%inline %{
//int print_args(char **argv) {
//    int i = 0;
//    while (argv[i]) {
//         printf("argv[%d] = %s\n", i,argv[i]);
//         i++;
//    }
//    return i;
//}

// Returns a char ** list
//char **get_args() {
//    static char *values[] = {"Dave","Mike","Susan","John","Michelle",0};
//    return &values[0];
//}
//%}



// Python typemaps
#ifdef SWIGPYTHON
// convert between python and C file handle
%typemap(in) FILE * {
  if (PyFile_Check($input)) /* check for undef */
        $1 = PyFile_AsFile($input);
  else  $1 = NULL;
}
#endif
