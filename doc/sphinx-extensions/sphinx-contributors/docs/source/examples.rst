Examples
========

The ``contributors`` directive renders the contributors of a given GitHub or local git repository in a Sphinx documentation page. This page shows the different options available.

Basic usage
-----------

List the contributors from a repository:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors

.. contributors:: dgarcia360/sphinx-contributors

Local repository
----------------

List the contributors from a local repository:

.. code-block:: rst

   .. contributors:: ../../
      :provider: local

.. contributors:: ../../
    :provider: local

With avatars
------------

Display contributor avatars:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :avatars:

.. contributors:: dgarcia360/sphinx-contributors
   :avatars:

Avatars only
------------

Show only clickable avatars with no text. Hover over an avatar to see the username:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :avatars_only:

.. contributors:: dgarcia360/sphinx-contributors
   :avatars_only:

With contribution counts
------------------------

Show how many contributions each user has made:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :contributions:

.. contributors:: dgarcia360/sphinx-contributors
   :contributions:

With real names
---------------

Show real names instead of usernames (fetched from GitHub profiles):

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :names:

.. contributors:: dgarcia360/sphinx-contributors
   :names:

Excluding users
---------------

Exclude bots or specific users:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :exclude: dependabot[bot], pre-commit-ci[bot]

.. contributors:: dgarcia360/sphinx-contributors
   :exclude: dependabot[bot], pre-commit-ci[bot]

Including additional contributors
---------------------------------

Manually add contributors not detected by the API (e.g., ``Co-authored-by`` contributors):

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :include: Peque

.. contributors:: dgarcia360/sphinx-contributors
   :include: Peque

Including anonymous contributors
--------------------------------

Also add contributors that are not linked to a github account:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :anonymous:

.. contributors:: dgarcia360/sphinx-contributors
   :anonymous:

Multiple repositories
---------------------

Merge contributors from multiple repositories into a single list. Contribution counts are summed for users who appear in more than one repository:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors sphinx-doc/sphinx
      :contributions:
      :limit: 10

.. contributors:: dgarcia360/sphinx-contributors sphinx-doc/sphinx
   :contributions:
   :limit: 10

Limiting and sorting
--------------------

Control how many contributors are shown and in what order:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :contributions:
      :limit: 3
      :order: ASC

.. contributors:: dgarcia360/sphinx-contributors
   :contributions:
   :limit: 3
   :order: ASC

Combining options
-----------------

Avatars, real names, and exclude list together:

.. code-block:: rst

   .. contributors:: dgarcia360/sphinx-contributors
      :avatars:
      :names:
      :exclude: dependabot[bot]

.. contributors:: dgarcia360/sphinx-contributors
   :avatars:
   :names:
   :exclude: dependabot[bot]

Local Repo with default avatars
-------------------------------

Local repository with default avatars, real names, and contribution counts together:

.. code-block:: rst

   .. contributors:: ../../
      :provider: local
      :contributions:
      :avatars:
      :names:
      :exclude: dependabot[bot]

.. contributors:: ../../
   :provider: local
   :contributions:
   :avatars:
   :names:
   :exclude: dependabot[bot]

For more information on how to style the list, see :doc:`Customization <customization>`.
