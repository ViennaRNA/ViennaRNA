/**********************************************/
/* BEGIN interface variable length arrays     */
/**********************************************/

%{
#define VAR_ARRAY_LINEAR    0U
#define VAR_ARRAY_TRI       1U

template <typename T>
struct var_array {
  size_t      length;
  T           *data;
  unsigned int type;
};

/****************/
/* Constructors */
/****************/
inline var_array<int> *
var_array_lin_int_new(size_t  length,
                      int     *data)
{
  if ((length) &&
      (data)) {
    var_array<int> *a = (var_array<int> *)vrna_alloc(sizeof(var_array<int>));
    a->length = length;
    a->data   = data;
    a->type   = VAR_ARRAY_LINEAR;
    return a;
  }

  return NULL;
}

inline var_array<int> *
var_array_tri_int_new(size_t  length,
                      int     *data)
{
  if ((length) &&
      (data)) {
    var_array<int> *a = (var_array<int> *)vrna_alloc(sizeof(var_array<int>));
    a->length = length;
    a->data   = data;
    a->type   = VAR_ARRAY_TRI;
    return a;
  }

  return NULL;
}

inline var_array<FLT_OR_DBL> *
var_array_lin_dbl_new(size_t      length,
                      FLT_OR_DBL  *data)
{
  if ((length) &&
      (data)) {
    var_array<FLT_OR_DBL> *a = (var_array<FLT_OR_DBL> *)vrna_alloc(sizeof(var_array<FLT_OR_DBL>));
    a->length = length;
    a->data   = data;
    a->type   = VAR_ARRAY_LINEAR;
    return a;
  }

  return NULL;
}

inline var_array<FLT_OR_DBL> *
var_array_tri_dbl_new(size_t      length,
                      FLT_OR_DBL  *data)
{
  if ((length) &&
      (data)) {
    var_array<FLT_OR_DBL> *a = (var_array<FLT_OR_DBL> *)vrna_alloc(sizeof(var_array<FLT_OR_DBL>));
    a->length = length;
    a->data   = data;
    a->type   = VAR_ARRAY_TRI;
    return a;
  }

  return NULL;
}

%}

template <typename T>
struct var_array {};

/***************************************/
/* Extensions for [] operator access   */
/***************************************/
%extend var_array {
  size_t __len__() const {
    switch($self->type) {
      case VAR_ARRAY_LINEAR:
        return $self->length + 1;

      case VAR_ARRAY_TRI:
        return ($self->length * ($self->length - 1)) / 2 + self->length + 1;

      default:
        return 0;
    }
  }

  const T __getitem__(int i) const throw(std::out_of_range) {
    switch ($self->type) {
      case VAR_ARRAY_LINEAR:
        if ((i < 0) ||
            (i > $self->length))
          throw std::out_of_range("out of bounds access");
        break;

      case VAR_ARRAY_TRI:
        if ((i < 0) ||
            (i > ($self->length * ($self->length - 1)) / 2 + self->length))
          throw std::out_of_range("out of bounds access");
        break;

      default:
        break;
    }

    return $self->data[i];
  }
};

%template (varArrayInt) var_array<int>;
%template (varArrayFLTorDBL) var_array<FLT_OR_DBL>;
