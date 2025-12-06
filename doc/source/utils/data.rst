Data Handling and Manipulation
==============================

It is not unusual that data needs to be tranformed before
being used in any method or algorithm (*pre-processing*) or
after it has been generated (*post-processing*). Here, we
provide a generic implementation to transform an array of
data into another array of data of the same size. Our
implementation utilizes a callback mechanism that is responsible
to transform any source value to a specific target value.

We already provide a handful of useful transformation functions
through our :ref:`utils/math:mathematical functions` API, e.g.
linear functions, logarithms, logistic functions, but also
bining and kernel density estimation (KDE). However, our
choice to implement callback mechanism based transformation
function(s) allows the users to implement any transformation function
themselfs. This enables a wide range of possible adaptations
of other methods throughout the ViennaRNA library that are build
on top of the transformation functions.

Below are all API symbols for the data transform implementation.


Data Transform API
------------------

.. doxygengroup:: data_utils
    :no-title:
