Installation
============

.. important:: sphinx-contributors is compatible with Sphinx 7.0+ using Python 3.10+.

Getting started
---------------

Follow these steps to add sphinx-contributors to your Sphinx project:

#. Install ``sphinx-contributors`` using PIP.

   .. code-block:: bash

      pip install sphinx-contributors

#. Add the extension to your Sphinx project ``conf.py`` file.

   .. code-block:: python

      extensions = ['sphinx_contributors']

#. Add a ``contributors`` directive to any ``.rst`` file.

   .. code-block:: rst

      .. contributors:: dgarcia360/sphinx-contributors
         :avatars:

   This renders all contributors with their avatars (use ``:limit:`` to show fewer):

   .. contributors:: dgarcia360/sphinx-contributors
      :avatars:

   See :doc:`Examples <examples>` for all available options.


Configuration
-------------

Use the ``sphinx_contributors_conf`` dictionary in your ``conf.py`` file to
overwrite any of the default values. The default configuration for this extension is

.. code-block:: python

    sphinx_contributors_conf = {
        "avatar_image_default" : "sc_user_avartar_default.svg",
        "avatar_size_default" : 150,
        "all_contributors" : False,
        "all_contributors_file" : ".all-contributorsrc",
        "all_contributors_imageSize_ignore" : False,
        "all_contributors_precedence" : [ "name", "avatar_url", "email", "profile" ],
    }

General settings
^^^^^^^^^^^^^^^^

The key/value pairs ``avatar_image_default`` and ``avatar_size_default`` can be used to change the default
avatar image if none can be obtained for a contributor, and the default avatar image size in pixels, respectively.

All Contributors settings
^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to contributors taken directly from the git repository, or passed via
the ``:include:`` directive option, auxiliary contributors may also be specified in
an `All Contributors <https://allcontributors.org>`_-style JSON file. This has the
advantage to **include contributors that do not show up anywhere in github or your
local git repository**. It also allows for a more fine-grained level of detail for
the actual contribution type. See https://allcontributors.org for more details of
contribution types.

Moreover, the JSON file allows for specifying a list of ``aliases`` for each contributor
that in turn helps to **merge contributors** that are the same person or entity but
show up under different ``login`` or ``email address``.

The ``all_contributors`` flag decides, whether or not contributor data from an
`All Contributors <https://allcontributors.org>`_-style JSON file should be included
in the output. The actual path to this JSON file can be specified with the
``all_contributors_file`` and it should be relative to your ``conf.py`` location. If
the ``imageSize`` attribute is present in the all-contributors JSON file, it will by
default overwrite the ``avatar_size_default`` value. Set the
``all_contributors_imageSize_ignore = True`` to always use the ``avatar_size_default``
instead. Last but not least, the ``all_contributors_precedence`` list specifies which
entries from the all-contributors JSON file must not be overwritten by any data
queried from git repositories.

.. note::

    If ``all_contributors = True``, the contributors taken from the all-contributors JSON file are collected first! Any git repositories are queried afterwards and the
    contributors (and their contributions) taken from those repositories are merged.

GitHub API rate limits
----------------------

This extension uses the GitHub REST API. Unauthenticated requests are limited to **60 per hour**. Each ``contributors`` directive makes at least one API call, and options like ``:names:`` or ``:include:`` make additional calls per contributor. If you hit the rate limit, contributor lists will appear empty until it resets.

To increase the limit to **5,000 requests per hour**, set a ``GITHUB_TOKEN`` environment variable with a `GitHub personal access token <https://github.com/settings/tokens>`_. A fine-grained token with no extra permissions is sufficient for public repositories.

.. code-block:: bash

   # Locally
   export GITHUB_TOKEN=ghp_your_token_here

For CI/CD environments like **ReadTheDocs**, add the token as an environment variable in your project settings. For **GitHub Actions**, use a repository secret:

.. code-block:: yaml

   - name: Build docs
     env:
       GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
     run: make -C docs dirhtml
