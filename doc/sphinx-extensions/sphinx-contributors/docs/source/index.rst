Overview
========

``sphinx-contributors`` is a Sphinx extension that automatically lists the people who have contributed to your project or Git repositories, right inside your documentation.

Features
---------

- List contributors from one or multiple **GitHub** or **local** ``git`` repositories in a single directive.
- Include additional contributors from an `All Contributors <https://allcontributors.org>`_-style JSON file
- Manually include contributors not detected by the API (e.g., ``Co-authored-by`` contributors).
- Display real names, avatars, and contribution types, commit counts, etc. if such data is available.
- Exclude bots and specific users.
- Merge contributors based on login or email aliases
- Sort by name, login, email address, contributions or commits and limit the number shown.

Requirements
---------------------

- GitHub public or local git repositories or local ``.all-contributorsrc`` JSON file.
- Python 3.10+ and Sphinx 7.0+.

Get started
-----------

.. toctree::
   :maxdepth: 2

   installation
   examples
   customization
   who-is-using-it
   reference
   contribute

Contributing
------------

We encourage public contributions!
Please review `CONTRIBUTING <https://sphinx-contributors.readthedocs.io/en/latest/contribute.html>`_ for details on our code of conduct and development process.

**THANK YOU TO THE CONTRIBUTORS OF SPHINX-CONTRIBUTORS!**

.. contributors:: dgarcia360/sphinx-contributors
   :avatars:
   :exclude: dependabot[bot]

License
-------

Copyright (c) 2018 - present David Garcia (`@dgarcia360 <https://twitter.com/dgarcia360>`_).

Licensed under the `MIT License <https://github.com/dgarcia360/sphinx-contributors/blob/main/LICENSE.md>`_.
