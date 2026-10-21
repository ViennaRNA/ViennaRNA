Reference
=========

Contributors directive
----------------------

.. rst:directive:: .. contributors:: username/repository [username/repository ...]

    One or more GitHub or local repositories, separated by spaces. When multiple repositories are specified, contributors are merged into a single list and contribution counts are summed for users who appear in more than one repository.

    .. code-block:: rst

       .. contributors:: dgarcia360/sphinx-contributors dgarcia360/other-repo
          :contributions:

    or

    .. code-block:: rst

       .. contributors:: ../../
          :provider: local
          :contributions:

    The latter assumes that your local git repository root resides two levels up from the document this directive has been placed.

    .. rst:directive:option:: provider
      :type: string

      By default, repositories specified as first argument of the directive are considered public github repositories. Set this option to ``local`` if they point to a local git repository instead. Defaults to ``github`` if not provided.

    .. rst:directive:option:: limit
      :type: integer

      Maximum number of contributors to list. If omitted, all contributors are shown.

    .. rst:directive:option:: order
      :type: string

      Results are unordered by default. Use this option to specify how the contributors should be ordered. The string argument is a whitespace delimited list of up to two keywords, one describing the
      attribute to use while ordering, the second to specify ascending or descending order. Use the following keywords to select an order attribute:

      - ``byContribution`` to order by the total number of contributions
      - ``byCommit`` to order by the total number of git commits
      - ``byName`` to order by the cotributors full name
      - ``byLogin`` to order by the contributors login name
      - ``byEmail`` to order by the contributors email address

      The keywords ``ASC`` or ``asc`` and ``DESC`` or ``desc`` can be used to change between ascending and descending order. Default is ``DESC``.

    .. rst:directive:option:: avatars
      :type: boolean

      Whether to include an image with the avatar returned by GitHub or read from all-contributors JSON data. Use CSS or the ``avatars_size`` option to customize the image size, since the returned images might vary.

    .. rst:directive:option:: avatars_only
      :type: boolean

      Show only avatar images with no usernames or contribution counts. Each avatar links to the contributor's (GitHub) profile, and the username is shown as a tooltip on hover. Implies ``:avatars:``.

    .. rst:directive:option:: avatar_size
      :type: int

      Set the avatar size in pixels. Default is specified by the global default setting ``sphinx_contributors_conf['avatar_size_default']``.

    .. rst:directive:option:: contributions
      :type: string

      If present, individual contributions are shown. The optional string argument is a whitespace delimited
      list of keywords to change the level of detail of the contribution output. By default, a detailed list
      of contributions is created. Use the ``brief`` keyword to create a brief list of contribution types.
      Add the ``emoji`` keyword to place an emoji icon before the contribution type. Use the ``emojis-only``
      keyword to omit the contribution type text and only render the emoji icon if one is available for the
      particular contribution type. Detailed output may consist of additional text describing the contribution
      and key/value pairs, such as total number of commits in git repositories.

    .. rst:directive:option:: names
      :type: boolean

      Whether to display real names (fetched from GitHub user profiles) instead of usernames. Falls back to the username if the user has not set a name on their GitHub profile.

    .. rst:directive:option:: include
      :type: string

      Comma separated GitHub usernames to add to the list of contributors. This is useful for including people who contributed via ``Co-authored-by`` commit trailers, since the GitHub REST API does not count those as contributors. Already listed contributors are not duplicated. For example: ``dgarcia360,otheruser``.

    .. rst:directive:option:: exclude
      :type: string

      Comma separated usernames to exclude from the list of contributors, for example: ``dependabot[bot],pre-commit-ci[bot]``.

    .. rst:directive:option:: anonymous
      :type: boolean

      Also show ``anonymous`` contributors, i.e. those that appear in the git commit history but are not registered with github

    .. rst:directive:option:: class_name
      :type: string

      By default, the directive uses the CSS class prefix ``sphinx-contributors``. Use this option to over-write the prefix to allow for individual CSS styling of the generated (HTML) output.

    .. note::

       The GitHub REST API used by this extension does not include contributors added via ``Co-authored-by`` commit trailers. Use the ``:include:`` option to manually add those contributors.
