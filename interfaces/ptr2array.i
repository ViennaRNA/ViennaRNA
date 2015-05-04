//
// SWIG pointer conversion and utility library
// 
// Dave Beazley
// April 19, 1997
//
// Perl5 specific implementation.   This file is included
// by the file ../pointer.i

%{

#ifdef WIN32
#undef isspace
#define isspace(c) (c == ' ')
#endif



/*------------------------------------------------------------------
  ptr2array(ptr,index,len,type = 0)

  Attempts to dereference a pointer value, just like ptrvalue(), but
  returns a reference to a list containing len values.  If type is given,
  it will try to use that type.  Otherwise, this function will attempt to
  "guess" the proper datatype by checking against all of the builtin C
  datatypes.
  ------------------------------------------------------------------ */

#ifdef PERL_OBJECT
static AV *_ptr2array(CPerlObj *pPerl,SV *_PTRVALUE, int index, int len, char *type) {
#define ptr2array(a,b,c,d) _ptr2array(pPerl,a,b,c,d)
#else
static AV *_ptr2array(SV *_PTRVALUE, int index, int len, char *type) {
#define ptr2array(a,b,c,d) _ptr2array(a,b,c,d)
#endif

  void     *ptr;
  SV *obj;
  SV      **svs = 0;
  AV        *av = 0;

  int i;
  
  svs = (SV **) malloc(len*sizeof(SV *));
  for (i=0; i<len; i++)
    svs[i] = sv_newmortal();
  
  if (SWIG_ConvertPtr(_PTRVALUE,&ptr,0,0)) {
    croak("Type error in ptr2array. Argument is not a valid pointer value.");
  } else {
    /* If no datatype was passed, try a few common datatypes first */
    if (!type) {

      /* No datatype was passed.   Type to figure out if it's a common one */

      if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_int_p,0) >= 0) {
        type = "int";
      } else if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_double_p,0) >= 0) {
        type = "double";
      } else if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_short_p,0) >= 0) {
        type = "short";
      } else if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_long_p,0) >= 0) {
        type = "long";
      } else if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_float_p,0) >= 0) {
        type = "float";
      } else if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_char_p,0) >= 0) {
        type = "char";
      } else if (SWIG_ConvertPtr(_PTRVALUE,&ptr,SWIG_POINTER_char_pp,0) >= 0) {
        type = "char *";
      } else {
        type = "unknown";
      }
    }
    
    if (!ptr) {
      croak("Unable to dereference NULL pointer.");
      return 0;
    }

    /* Now we have a datatype.  Try to figure out what to do about it */
    if (strcmp(type,"int") == 0) {
      for (i = 0; i < len; ++i) 
	sv_setiv(svs[i], (IV) *(((int *) ptr) + index+i));
    } else if (strcmp(type,"double") == 0) {
      for (i = 0; i < len; ++i) 
	sv_setnv(svs[i], (double) *(((double *) ptr) + index+i));
    } else if (strcmp(type,"short") == 0) {
      for (i = 0; i < len; ++i) 
	sv_setiv(svs[i],(IV) *(((short *) ptr) + index+i));
    } else if (strcmp(type,"long") == 0) {
      for (i = 0; i < len; ++i) 
	sv_setiv(svs[i],(IV) *(((long *) ptr) + index+i));
    } else if (strcmp(type,"float") == 0) {
      for (i = 0; i < len; ++i) 
	sv_setnv(svs[i],(double) *(((float *) ptr)+index+i));
    } else if (strcmp(type,"char") == 0) {
      for (i = 0; i < len; ++i) 
	sv_setpv(svs[i],((char *) ptr)+index+i);
    } else if (strcmp(type,"char *") == 0) {
      for (i = 0; i < len; ++i) {
	char *c = *(((char **) ptr)+index+i);
	if (c) 
	  sv_setpv(svs[i],c);
	else 
	  sv_setpv(svs[i],"NULL");
      }
    } else {
      croak("Unable to dereference unsupported datatype.");
      len = 0;
    }
  }

  av = av_make(len,svs);
  free(svs);
  return av;
}

%}

%typemap(perl5, out) AV *ptr2array
{
  $result = newRV_noinc((SV*) $1);
  sv_2mortal($result);
  argvi++;
}

 AV *ptr2array(SV *ptr, int index = 0, int len = 1, char *type = 0); 

