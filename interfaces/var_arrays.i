/**********************************************/
/* BEGIN interface variable length arrays     */
/**********************************************/

%{
#include <sstream>

#define VAR_ARRAY_LINEAR      1U
#define VAR_ARRAY_TRI         2U
#define VAR_ARRAY_SQR         4U
#define VAR_ARRAY_ONE_BASED   8U
#define VAR_ARRAY_OWNED       16U

template <typename T>
struct var_array {
  size_t        length;
  T             *data;
  unsigned int  type;
};


template <typename T>
std::string
var_array_type_str(var_array<T> *a)
{
  std::ostringstream out;

  if (a->type & VAR_ARRAY_LINEAR)
    out << "RNA.VAR_ARRAY_LINEAR";
  else  if (a->type & VAR_ARRAY_TRI)
    out << "RNA.VAR_ARRAY_TRI";
  else  if (a->type & VAR_ARRAY_SQR)
    out << "RNA.VAR_ARRAY_SQR";

  if (a->type & VAR_ARRAY_ONE_BASED)
    out << " | RNA.VAR_ARRAY_ONE_BASED";

  return std::string(out.str());
}


/****************/
/* Constructors */
/****************/
inline var_array<unsigned char> *
var_array_uchar_new(size_t        length,
                    unsigned char *data,
                    unsigned int  type)
{
  var_array<unsigned char> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<unsigned char> *)vrna_alloc(sizeof(var_array<unsigned char>));
    a->length = length;
    a->data   = data;
    a->type   = type;
  }

  return a;
}


inline var_array<char> *
var_array_uchar_new(size_t        length,
                    char          *data,
                    unsigned int  type)
{
  var_array<char> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<char> *)vrna_alloc(sizeof(var_array<char>));
    a->length = length;
    a->data   = data;
    a->type   = type;
  }

  return a;
}


inline var_array<short> *
var_array_short_new(size_t        length,
                    short         *data,
                    unsigned int  type)
{
  var_array<short> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<short> *)vrna_alloc(sizeof(var_array<short>));
    a->length = length;
    a->data   = data;
    a->type   = type;
  }

  return a;
}


inline var_array<int> *
var_array_int_new(size_t        length,
                  int           *data,
                  unsigned int  type)
{
  var_array<int> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<int> *)vrna_alloc(sizeof(var_array<int>));
    a->length = length;
    a->data   = data;
    a->type   = type;
  }

  return a;
}


inline var_array<unsigned int> *
var_array_uint_new(size_t       length,
                   unsigned int *data,
                   unsigned int type)
{
  var_array<unsigned int> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<unsigned int> *)vrna_alloc(sizeof(var_array<unsigned int>));
    a->length = length;
    a->data   = data;
    a->type   = type;
  }

  return a;
}

inline var_array<FLT_OR_DBL> *
var_array_dbl_new(size_t        length,
                  FLT_OR_DBL    *data,
                  unsigned int  type)
{
  var_array<FLT_OR_DBL> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<FLT_OR_DBL> *)vrna_alloc(sizeof(var_array<FLT_OR_DBL>));
    a->length = length;
    a->data   = data;
    a->type   = type;
  }

  return a;
}

%}

%nodefaultctor var_array;
%nodefaultdtor var_array;

template <typename T>
struct var_array {};


%{
  inline var_array<short> *
  var_array_short_from_int(std::vector<int> d,
                           unsigned int   type) {
    var_array<short>  *a  = NULL;
    size_t            n   = d.size();

    if (n > 0) {
      short *sd = (short *)vrna_alloc(sizeof(short) * n);

      for (size_t i = 0; i < n; i++)
        sd[i] = (short)d[i];

      a = var_array_short_new(n, sd, type | VAR_ARRAY_OWNED);

    }

    return a;
  }

%}

/***************************************/
/* Constructor/Destructor              */
/***************************************/
%extend var_array {
  var_array(std::vector<T> d,
            unsigned int   type) {
    var_array<T> *a = NULL;
    size_t        n = d.size();

    if (n > 0) {
      a = (var_array<T> *)vrna_alloc(sizeof(var_array<T>));
      a->data   = (T *)vrna_alloc(sizeof(T) * n);
      memcpy(a->data, &(d[0]), sizeof(T) * n);

      a->length = n;
      a->type   = type | VAR_ARRAY_OWNED;
    }

    return a;
  }

  ~var_array() {
    if ($self->type & VAR_ARRAY_OWNED)
      free($self->data);
    free($self);
  }
};


%extend var_array<short> {
  var_array(std::vector<int> data,
            unsigned int   type = VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED) {
    var_array<short> *a = NULL;
    size_t        n = data.size();

    if (n > 0) {
      a = (var_array<short> *)vrna_alloc(sizeof(var_array<short>));

      a->length = n;

      if (type & VAR_ARRAY_ONE_BASED)
        a->length--;

      a->data   = (short *)vrna_alloc(sizeof(short) * n);

      for (size_t i = 0; i < n; i++)
        a->data[i] = (short)data[i];

      a->type   = type | VAR_ARRAY_OWNED;
    }

    return a;
  }
};


/***************************************/
/* Extensions for [] operator access   */
/***************************************/
%extend var_array {
  size_t __len__() const {
    size_t n = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      n += 1;

    if ($self->type & VAR_ARRAY_TRI)
      n = ((n - 1) * (n - 2)) / 2 + n;
    else if ($self->type & VAR_ARRAY_SQR)
      n = n * n;

    return n;
  }

  const T __getitem__(int i) const throw(std::out_of_range) {
    size_t max_i = $self->length - 1;
    
    if ($self->type & VAR_ARRAY_ONE_BASED)
      max_i += 1;

    if ($self->type & VAR_ARRAY_TRI)
      max_i = ((max_i - 1) * (max_i - 2)) / 2 + max_i;
    else if ($self->type & VAR_ARRAY_SQR)
      max_i = max_i * max_i;

    if ((i < 0) ||
        (i > max_i))
      throw std::out_of_range("out of bounds access");

    return $self->data[i];
  }

  const T __setitem__(int i, const T d) const throw(std::out_of_range) {
    size_t max_i = $self->length - 1;
    
    if ($self->type & VAR_ARRAY_ONE_BASED)
      max_i += 1;

    if ($self->type & VAR_ARRAY_TRI)
      max_i = ((max_i - 1) * (max_i - 2)) / 2 + max_i;
    else if ($self->type & VAR_ARRAY_SQR)
      max_i = max_i * max_i;

    if ((i < 0) ||
        (i > max_i))
      throw std::out_of_range("out of bounds access");

    return $self->data[i] = d;
  }
};


#ifdef SWIGPYTHON
%extend var_array {
  std::string
  __str__()
  {
    size_t n = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      n += 1;

    if ($self->type & VAR_ARRAY_TRI)
      n = ((n - 1) * (n - 2)) / 2 + n;
    else if ($self->type & VAR_ARRAY_SQR)
      n = n * n;

    std::ostringstream out;
    out << "{ data: [" << $self->data[0];
    for (size_t i = 1; i < n; i++)
      out << ", " << $self->data[i];
    out << "], ";
    out << "type: " << var_array_type_str($self);
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}

};

#endif



%template (varArrayUChar) var_array<unsigned char>;
%template (varArrayChar) var_array<char>;
%template (varArrayShort) var_array<short>;
%template (varArrayUInt) var_array<unsigned int>;
%template (varArrayInt) var_array<int>;
%template (varArrayFLTorDBL) var_array<FLT_OR_DBL>;

%constant unsigned int VAR_ARRAY_LINEAR     = VAR_ARRAY_LINEAR;
%constant unsigned int VAR_ARRAY_TRI        = VAR_ARRAY_TRI;
%constant unsigned int VAR_ARRAY_SQR        = VAR_ARRAY_SQR;
%constant unsigned int VAR_ARRAY_ONE_BASED  = VAR_ARRAY_ONE_BASED;
%constant unsigned int VAR_ARRAY_OWNED      = VAR_ARRAY_OWNED;
