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


/* determine length n given array_size elements of the respective
   upper triangular square matrix. Here, we assume that the largest
   entry in the array resides at position (n*(n-1))/2 + n. So the
   actual size of the array is (n*(n-1))/2 + n + 1, since arrays 
   are always 0-based
*/
inline size_t
var_array_tri_size(size_t array_size)
{
  return (size_t)floor((sqrt(8 * array_size - 7) - 1) / 2);
}


/* compute number of entries in the upper-triangular n*n square
   matrix with dimension n
*/
inline size_t
var_array_data_size_tri(size_t n)
{
  return (n * (n - 1)) / 2 + n + 1;
}


/* determine length n given array_size elements of the respective
   square matrix
*/
inline size_t
var_array_sqr_size(size_t array_size)
{
  return (size_t)(sqrt(array_size - 1));
}


/* compute number of entries in the n*n square
   matrix with dimension n
*/
inline size_t
var_array_data_size_sqr(size_t n)
{
  return n * n + 1;
}


/****************/
/* Constructors */
/****************/
template <typename T>
inline var_array<T> *
var_array_new(size_t        length,
              T             *data,
              unsigned int  type)
{
  var_array<T> *a = NULL;

  if ((length) &&
      (data)) {
    a         = (var_array<T> *)vrna_alloc(sizeof(var_array<T>));
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

      if (type & VAR_ARRAY_TRI)
        n = var_array_tri_size(n);
      else if (type & VAR_ARRAY_SQR)
        n = var_array_sqr_size(n);
      else if ((type & VAR_ARRAY_LINEAR) &&
               (type & VAR_ARRAY_ONE_BASED))
        n -= 1;

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
      short *data_s   = (short *)vrna_alloc(sizeof(short) * n);

      for (size_t i = 0; i < n; i++)
        data_s[i] = (short)data[i];

      if (type & VAR_ARRAY_TRI)
        n = var_array_tri_size(n);
      else if (type & VAR_ARRAY_SQR)
        n = var_array_sqr_size(n);
      else if ((type & VAR_ARRAY_LINEAR) &&
               (type & VAR_ARRAY_ONE_BASED))
        n -= 1;

      a = var_array_new(n, data_s, type | VAR_ARRAY_OWNED);
    }

    return a;
  }
};


/***************************************/
/* Extensions for [] operator access   */
/***************************************/
%extend var_array {
  size_t size() const {
    size_t n = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      n += 1;

    if ($self->type & VAR_ARRAY_TRI)
      n = var_array_data_size_tri(n - 1);
    else if ($self->type & VAR_ARRAY_SQR)
      n = var_array_data_size_sqr(n);

    return n;
  }

  size_t __len__() const {
    size_t n = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      n += 1;

    if ($self->type & VAR_ARRAY_TRI)
      n = var_array_data_size_tri(n - 1);
    else if ($self->type & VAR_ARRAY_SQR)
      n = var_array_data_size_sqr(n);

    return n;
  }

  const T get(int i) const throw(std::out_of_range) {
    size_t max_i = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      max_i += 1;

    if ($self->type & VAR_ARRAY_TRI)
      max_i = var_array_data_size_tri(max_i - 1);
    else if ($self->type & VAR_ARRAY_SQR)
      max_i = var_array_data_size_sqr(max_i);

    /* also handle negative indices here */
    if (((i < 0) && (-i >= max_i)) ||
        ((i >= 0) && (i >= max_i)))
      throw std::out_of_range("out of bounds access");

    if (i < 0)
      i = max_i + i;

    return $self->data[i];
  }

  const T __getitem__(int i) const throw(std::out_of_range) {
    size_t max_i = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      max_i += 1;

    if ($self->type & VAR_ARRAY_TRI)
      max_i = var_array_data_size_tri(max_i - 1);
    else if ($self->type & VAR_ARRAY_SQR)
      max_i = var_array_data_size_sqr(max_i);

    /* also handle negative indices here */
    if (((i < 0) && (-i >= max_i)) ||
        ((i >= 0) && (i >= max_i)))
      throw std::out_of_range("out of bounds access");

    if (i < 0)
      i = max_i + i;

    return $self->data[i];
  }

  const T __setitem__(int i, const T d) const throw(std::out_of_range) {
    size_t max_i = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      max_i += 1;

    if ($self->type & VAR_ARRAY_TRI)
      max_i = var_array_data_size_tri(max_i - 1);
    else if ($self->type & VAR_ARRAY_SQR)
      max_i = var_array_data_size_sqr(max_i);

    /* also handle negative indices here */
    if (((i < 0) && (-i >= max_i)) ||
        ((i >= 0) && (i >= max_i)))
      throw std::out_of_range("out of bounds access");

    if (i < 0)
      i = max_i + i;

    return $self->data[i] = d;
  }
};


#ifdef SWIGPYTHON
%pythoncode %{
class var_array_Iterator:
    def __init__(self, var_arr):
        self.index = 0
        self.var_arr = var_arr
    def __iter__(self):
        return self
    def next(self):
        if self.index >= self.var_arr.size():
            raise StopIteration;
        ret = self.var_arr.get(self.index)
        self.index += 1
        return ret
    __next__ = next
%}

%extend var_array {
  var_array<T> * __getitem__(PyObject *ob) {
    size_t        n, real_i, max_i;
    unsigned int  type;
    Py_ssize_t    start, stop, step, slice_l, l;
    var_array<T>  *a;

    if (PySlice_Check(ob)) {
      a = NULL;
      max_i = $self->length;

      if ($self->type & VAR_ARRAY_ONE_BASED)
        max_i += 1;

      if ($self->type & VAR_ARRAY_TRI)
        max_i = var_array_data_size_tri(max_i - 1);
      else if ($self->type & VAR_ARRAY_SQR)
        max_i = var_array_data_size_sqr(max_i);

      l = max_i;


      if (PySlice_GetIndicesEx(ob, l, &start, &stop, &step, &slice_l) == 0) {
        real_i  = 0;
        type    = $self->type;
        n       = slice_l;

        if ((type & VAR_ARRAY_LINEAR) &&
            (type & VAR_ARRAY_ONE_BASED) &&
            (start != 0))
          type &= ~VAR_ARRAY_ONE_BASED;

        if (n > 0) {
          a = (var_array<T> *)vrna_alloc(sizeof(var_array<T>));
          a->data   = (T *)vrna_alloc(sizeof(T) * n);
          /* copy slice data */
          for (Py_ssize_t i = start; i <= stop; i += step, real_i++)
            a->data[real_i] = $self->data[i];

          if ((type & VAR_ARRAY_LINEAR) &&
               (type & VAR_ARRAY_ONE_BASED)) {
            n -= 1;
            a->data[0] = n;
          }

          a->length = n;
          a->type   = type | VAR_ARRAY_OWNED;
        }
      } else {
        return NULL;
      }

    } else {
      PyErr_SetString(PyExc_TypeError, "Expected integer or slice object");
      return 0;
    }

    return a;
  }

  std::string
  __str__()
  {
    size_t n = $self->length;

    if ($self->type & VAR_ARRAY_ONE_BASED)
      n += 1;

    if ($self->type & VAR_ARRAY_TRI)
      n = var_array_data_size_tri(n - 1);
    else if ($self->type & VAR_ARRAY_SQR)
      n = var_array_data_size_tri(n);

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

def __iter__(self):
    return var_array_Iterator(self)

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
