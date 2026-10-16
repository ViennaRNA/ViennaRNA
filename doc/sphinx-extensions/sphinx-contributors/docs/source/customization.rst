Customization
-------------

The directive ``contributors`` generates an ``<ul>`` node with class ``sphinx-contributors``.
This makes the list is highly customizable through CSS.

.. tip:: See how to `add custom CSS to Sphinx Documentation <https://docs.readthedocs.io/en/stable/guides/adding-custom-css.html>`_.

For example, the following CSS snippet makes the images round:

.. code-block:: css

  .sphinx-contributors img {
      border-radius: 50%;
  }

.. rst-class:: custom-contributors

.. contributors:: dgarcia360/sphinx-contributors
    :avatars:
    :exclude: dgarcia360,sphinx
    :limit: 5
    :order: DESC

You can also adjust the spacing between avatars by overriding the ``gap`` property. Increase it for more spacing, or decrease it for a tighter layout:

.. code-block:: css

  /* More spacing */
  .sphinx-contributors--avatars .sphinx-contributors_list {
      gap: 40px;
  }

  /* Tighter layout */
  .sphinx-contributors--avatars .sphinx-contributors_list {
      gap: 10px;
  }
