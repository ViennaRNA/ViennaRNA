/**********************************************/
/* BEGIN interface variable length arrays     */
/**********************************************/

%nodefaultctor var_array_tri_int_t;
%nodefaultctor var_array_lin_int_t;
%nodefaultctor var_array_tri_dbl_colwise_t;
%nodefaultctor var_array_tri_dbl_rowwise_t;
%nodefaultctor var_array_lin_dbl_t;

%{
/* (upper) triangular matrix of type int */
typedef struct {
  size_t length;
  int    *data;
} var_array_tri_int_t;

/* linear array of type int */
typedef struct {
  size_t length;
  int    *data;
} var_array_lin_int_t;

/* upper triangular matrix of type FLT_OR_DBL (access column wise) */
typedef struct {
  size_t      length;
  FLT_OR_DBL  *data;
} var_array_tri_dbl_colwise_t;

/* upper triangular matrix of type FLT_OR_DBL (access row wise) */
typedef struct {
  size_t      length;
  FLT_OR_DBL  *data;
} var_array_tri_dbl_rowwise_t;

/* linear array of type FLT_OR_DBL */
typedef struct {
  size_t      length;
  FLT_OR_DBL  *data;
} var_array_lin_dbl_t;

/****************/
/* Constructors */
/****************/
inline var_array_lin_int_t *
var_array_lin_int_new(size_t  length,
                      int     *data)
{
  if ((length) &&
      (data)) {
    var_array_lin_int_t *a = (var_array_lin_int_t *)vrna_alloc(sizeof(var_array_lin_int_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  return NULL;
}

inline var_array_tri_int_t *
var_array_tri_int_new(size_t  length,
                      int     *data)
{
  if ((length) &&
      (data)) {
    var_array_tri_int_t *a = (var_array_tri_int_t *)vrna_alloc(sizeof(var_array_tri_int_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  return NULL;
}

inline var_array_lin_dbl_t *
var_array_lin_dbl_new(size_t      length,
                      FLT_OR_DBL  *data)
{
  if ((length) &&
      (data)) {
    var_array_lin_dbl_t *a = (var_array_lin_dbl_t *)vrna_alloc(sizeof(var_array_lin_dbl_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  return NULL;
}

inline var_array_tri_dbl_colwise_t *
var_array_tri_dbl_colwise_new(size_t      length,
                              FLT_OR_DBL  *data)
{
  if ((length) &&
      (data)) {
    var_array_tri_dbl_colwise_t *a = (var_array_tri_dbl_colwise_t *)vrna_alloc(sizeof(var_array_tri_dbl_colwise_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  return NULL;
}

inline var_array_tri_dbl_rowwise_t *
var_array_tri_dbl_rowwise_new(size_t      length,
                              FLT_OR_DBL  *data)
{
  if ((length) &&
      (data)) {
    var_array_tri_dbl_rowwise_t *a = (var_array_tri_dbl_rowwise_t *)vrna_alloc(sizeof(var_array_tri_dbl_rowwise_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  return NULL;
}

%}

typedef struct {} var_array_tri_int_t;
typedef struct {} var_array_lin_int_t;
typedef struct {} var_array_tri_dbl_colwise_t;
typedef struct {} var_array_tri_dbl_rowwise_t;
typedef struct {} var_array_lin_dbl_t;

/***************************************/
/* Extensions for [] operator access   */
/***************************************/

%extend var_array_tri_int_t {
  var_array_tri_int_t(size_t length, int *data)
  {
    var_array_tri_int_t *a = (var_array_tri_int_t *)vrna_alloc(sizeof(var_array_tri_int_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  size_t __len__() const {
    return ($self->length * ($self->length - 1)) / 2 + self->length + 1;
  }

  const int __getitem__(int i) const throw(std::out_of_range) {
    if ((i < 0) ||
        (i > ($self->length * ($self->length - 1)) / 2 + self->length))
      throw std::out_of_range("out of bounds access");

    return $self->data[i];
  }
};

%extend var_array_lin_int_t {
  var_array_lin_int_t(size_t length, int *data)
  {
    var_array_lin_int_t *a = (var_array_lin_int_t *)vrna_alloc(sizeof(var_array_lin_int_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  size_t __len__() const {
    return $self->length + 1;
  }

  const int __getitem__(int i) const throw(std::out_of_range) {
    if ((i < 0) ||
        (i > $self->length))
      throw std::out_of_range("out of bounds access");

    return $self->data[i];
  }
};

%extend var_array_tri_dbl_colwise_t {
  var_array_tri_dbl_colwise_t(size_t length, FLT_OR_DBL *data)
  {
    var_array_tri_dbl_colwise_t *a = (var_array_tri_dbl_colwise_t *)vrna_alloc(sizeof(var_array_tri_dbl_colwise_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  size_t __len__() const {
    return ($self->length * ($self->length - 1)) / 2 + self->length + 1;
  }

  const FLT_OR_DBL __getitem__(int i) const throw(std::out_of_range) {
    if ((i < 0) ||
        (i > ($self->length * ($self->length - 1)) / 2 + self->length))
      throw std::out_of_range("out of bounds access");

    return $self->data[i];
  }
};

%extend var_array_tri_dbl_rowwise_t {
  var_array_tri_dbl_rowwise_t(size_t length, FLT_OR_DBL *data)
  {
    var_array_tri_dbl_rowwise_t *a = (var_array_tri_dbl_rowwise_t *)vrna_alloc(sizeof(var_array_tri_dbl_rowwise_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  size_t __len__() const {
    return ($self->length * ($self->length - 1)) / 2 + self->length + 1;
  }

  const FLT_OR_DBL __getitem__(int i) const throw(std::out_of_range) {
    if ((i < 0) ||
        (i > ($self->length * ($self->length - 1)) / 2 + self->length))
      throw std::out_of_range("out of bounds access");

    return $self->data[i];
  }
};

%extend var_array_lin_dbl_t {
  var_array_lin_dbl_t(size_t length, FLT_OR_DBL *data)
  {
    var_array_lin_dbl_t *a = (var_array_lin_dbl_t *)vrna_alloc(sizeof(var_array_lin_dbl_t));
    a->length = length;
    a->data   = data;
    return a;
  }

  size_t __len__() const {
    return $self->length + 1;
  }

  const FLT_OR_DBL __getitem__(int i) const throw(std::out_of_range) {
    if ((i < 0) ||
        (i > $self->length))
      throw std::out_of_range("out of bounds access");

    return $self->data[i];
  }
};
