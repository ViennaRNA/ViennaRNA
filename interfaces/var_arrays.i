/**********************************************/
/* BEGIN interface variable length arrays     */
/**********************************************/

%{
#define VAR_ARRAY_LINEAR    0U
#define VAR_ARRAY_TRI       1U
#define VAR_ARRAY_SQR       2U
#define VAR_ARRAY_ONE_BASED 4U

template <typename T>
struct var_array {
  size_t        length;
  T             *data;
  unsigned int  type;
};

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

template <typename T>
struct var_array {};

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
};

%template (varArrayUChar) var_array<unsigned char>;
%template (varArrayChar) var_array<char>;
%template (varArrayUInt) var_array<unsigned int>;
%template (varArrayInt) var_array<int>;
%template (varArrayFLTorDBL) var_array<FLT_OR_DBL>;
