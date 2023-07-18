Arrays
======

Interface for an abstract implementation of an array data structure

Arrays of a particular ``Type`` are defined and initialized using the following code:

.. code:: c

  vrna_array(Type)  my_array;
  vrna_array_init(my_array);

or equivalently:

.. code:: c

  vrna_array_make(Type, my_array);

Dynamic arrays can be used like regular pointers, i.e. elements are simply
addressed using the ``[]`` operator, e.g.:

.. code:: c

  my_array[1] = 42;

Using the :c:macro:`vrna_array_append` macro, items can be safely appended and the
array will grow accordingly if required:

.. code:: c

  vrna_array_append(my_array, item);

Finally, memory occupied by an array must be released using the :c:macro:`vrna_array_free`
macro:

.. code:: c

  vrna_array_free(my_array);

Use the :c:macro:`vrna_array_size` macro to get the number of items stored in an array,
e.g. for looping over its elements:

.. code:: c

  // define and initialize
  vrna_array_make(int, my_array);

  // append some items
  vrna_array_append(my_array, 42);
  vrna_array_append(my_array, 23);
  vrna_array_append(my_array, 5);

  // loop over items and print
  for (size_t i = 0; i < vrna_array_size(my_array); i++)
    printf("%d\n", my_array[i]);

  // release memory of the array
  vrna_array_free(my_array);


Under the hood, arrays are preceded by a header that actually stores the
number of items they contain and the capacity of elements they are able
to store.  The general ideas for this implementation are taken from
Ginger Bill's C Helper Library (public domain).

.. doxygengroup:: array_utils
