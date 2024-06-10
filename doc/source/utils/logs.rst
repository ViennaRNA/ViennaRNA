Log Messages
============

The ViennaRNA Package comes with a log message system that
allows for filtering and (re-)routing of all messages issued
by the library.

The default log message system prints messages that meet a certain
log level to a file pointer that defaults to *stderr*. This file
pointer, the log level threshold and the content of the messages
can be adapted with the functions below.

In addition, the logging system allows for adding user-defined
callbacks that can also retrieve the log messages for further
processing. Along with that, the default logging mechanism can
be turned off entirely such that only the callbacks process the
log messages.

Below are all API symbols for interacting with the ViennaRNA
logging system.


Logging System API
------------------

.. doxygengroup:: utils_log
    :no-title:
