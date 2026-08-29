Structured Domains
==================

Add and modify structured domains to the RNA folding grammar.

.. contents:: Table of Contents
    :local:
    :depth: 2


Introduction
------------

This module provides the tools to add and modify structured domains to
the production rules of the RNA folding grammar.


Usually, structured domains represent self-enclosed structural modules that
exhibit a more or less complex base pairing pattern. This can be more or less
well-defined 3D motifs, such as *G-Quadruplexes*, or loops with additional
non-canonical base pair interactions, such as *kink-turns*.

.. note::

  Currently, our implementation only provides the specialized case of *G-Quadruplexes*.


Structured Domains API
----------------------

.. doxygengroup:: domains_struc
    :no-title:
