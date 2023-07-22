Callback Functions
==================

Many functions in *RNAlib* support so-called *callback mechanisms* to propagate
the computation results. In essence, this means that a user defines a function
that takes computation results as input and does whatever should be done with
these results. Then, this function is provided to our algorithms in *RNAlib*
such that they can call the function and provide corresponding data as soon
as it has been computed.

Why Callbacks?
--------------

Using callback mechanisms, our library enables users not only to retrieve computed
data without the need for parsing complicated data structures, but also allows one
to tweak our implementation to do additional tasks without the requirement of a
re-implementation of basic algorithms.

Our implementation of the callback mechanisms always follows the same scheme:

The user ...

* ... defines a function that complies with the interface we've defined, and
* ... passes a pointer to said function to our implementations

In addition to the specific arguments of our callback interfaces, virtually all callbacks
receive an additional @em pass-through-pointer as their last argument. This enables one
to ...

* ... encapsulate data, and
* ... provide thread-safe operations,

since this pointer is simply passed through by our library functions. It may therefore hold
the address of an arbitrary, user-defined data structure.

Scripting Language Support
--------------------------

Our callback mechanisms also work in a cross-language specific manner. This
means that users can write functions (or methods) in a scripting language,
such as Python, and use them as callbacks for our C-library. Our scripting
language interface will take care of transforming the relevant data structures
from the target language into C and vice-versa. Whenever our algorithms trigger
the callback, they will then be calling the actual scripting language function
and provide the corresponding data directly to them.

.. warning::

  Keep in mind that the translation between the scripting language and C
  involves many extra function calls to prepare and evaluate the corresponding
  data structures. This in turn impacts the runtime of our algorithms that
  can be substantial. For instance providing callbacks for the hard- or
  soft constraints framework from a scripting language can lead to a slow-down
  of up to a factor of 10.


Available Callbacks
-------------------

Below, you find an enumeration of the individual callback functions that are available
in *RNAlib*:

.. doxygenpage:: callbacks
    :no-title:
