=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[0.3.0] - 20 March 2026
========================

Added
-----

- `#36 <https://github.com/dgarcia360/sphinx-contributors/issues/36>`_: New ``:names:`` option to display real names fetched from GitHub user profiles instead of usernames.
- `#35 <https://github.com/dgarcia360/sphinx-contributors/issues/35>`_: New ``:include:`` option to manually add contributors not detected by the GitHub REST API (e.g., ``Co-authored-by`` contributors).
- `#6 <https://github.com/dgarcia360/sphinx-contributors/issues/6>`_: Support for multiple repositories in a single directive. Contributors are merged and contribution counts are summed.
- `#37 <https://github.com/dgarcia360/sphinx-contributors/pull/37>`_: Support for whitespaces in the ``:exclude:`` list.
- `#34 <https://github.com/dgarcia360/sphinx-contributors/pull/34>`_: Added ``:exclude:`` option to reference docs.
- Support for ``GITHUB_TOKEN`` environment variable to authenticate API requests and increase the rate limit from 60 to 5,000 requests per hour. The token is only used at build time and is never included in the generated output.
- New ``:avatars_only:`` option to show only clickable avatar images with no text. Useful for compact contributor grids.
- GitHub API pagination: repositories with more than 100 contributors are now fully supported.
- New examples page in the documentation showcasing all directive options.

Changed
-------

- **Breaking**: The ``:limit:`` option no longer defaults to 10. If omitted, all contributors are shown. Set ``:limit: 10`` explicitly to restore the previous behavior.

Updated
-------

- `#40 <https://github.com/dgarcia360/sphinx-contributors/pull/40>`_: Python support updated to 3.10-3.14.
- Minimum Sphinx version bumped from 5 to 7.
- `#33 <https://github.com/dgarcia360/sphinx-contributors/pull/33>`_: Bumped ``actions/checkout`` from 3 to 4.
- CI matrix optimized: older Sphinx versions tested only on Python 3.12.
- Enabled parallel read/write safe for the extension.
- Avatar layout now aligns left-to-right instead of being distributed evenly.

Fixed
-----

- `#38 <https://github.com/dgarcia360/sphinx-contributors/pull/38>`_: Fixed typo in docs (whithout -> without).
- Avatar images now render at a consistent square size.

[0.2.7] - 10 January 2023
=========================

Added
-----

- [#27](https://github.com/dgarcia360/sphinx-contributors/pull/27) Sphinx extension classifier - thanks @jdillard!

Updated
-------

- [#28](https://github.com/dgarcia360/sphinx-contributors/pull/28): Black dependency for pre-commit - thanks @jdillard!

[0.2.6] - 26 January 2022
=========================

Added
-----

- Default CSS for lists.

[0.2.4] - 21 January 2022
=========================

Added
-----

- Linting with pre-commit.

Updated
-------

- Breaking change: The package is now named ``sphinx-contributors``.
- Breaking change: The required Sphinx version is at least 3.0.
- Replaced Travis CI for GitHub Actions.
- Replaced ``setup.py`` for ``pyproject.toml``.

[0.2.3] - 25 November 2021
==========================

Added
-----

* [#1](https://github.com/dgarcia360/sphinx-contributors/issues/1): The directive supports an option to render contributor's avatars - thanks @segfaultxavi!
* [#3])(https://github.com/dgarcia360/sphinx-contributors/issues/3): The docs now explain how to render contributors in a grid view.
* [#4](https://github.com/dgarcia360/sphinx-contributors/issues/4): The directive supports an option to hide the number of commits done by a contributor.
* [#13](https://github.com/dgarcia360/sphinx-contributors/issues/13): Docs for this extesion are now hosted on ReadTheDocs.

Updated
-------

* The extension now raises a warning instead of an error if the repository cannot be retrieved. This is especially useful when trying to build docs without internet connection.

[0.2.2] - 03 August 2019
========================

Added
-----

* [#2](https://github.com/dgarcia360/sphinx-contributors/issues/2): The directive supports an option to limit the number of contributors.
* The directive supports an option to sort contributors per number of commits.
* The directive supports an option to exclude usernames from the contributors' list.

Changed
--------

* Before, the repository name was passed as an option. Now, the directive accepts the repository name as an argument.


[0.1.0] - 18 July 2018
=======================

Added
-----

* Initial code release
