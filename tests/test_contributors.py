"""
Tests for sphinx-contributors extension.
(c) 2022 - present David Garcia (@dgarcia360)
# This code is licensed under MIT license (see LICENSE.md for details)
"""

from pathlib import Path
from textwrap import dedent
from unittest.mock import MagicMock, patch

from sphinx.application import Sphinx

from sphinx_contributors.contributors import Contributor, ContributorsRepository

_CLASS_NAME = "sphinx_contributors"


def test_contributor_build() -> None:
    """
    Initializing a contributor.
    """
    contributor = Contributor("dgarcia360", "http://#", contributions=10).build(
        _CLASS_NAME,
        { 'show_contributions' : True }
    )
    assert contributor.astext() == "dgarcia360\n\nCode:\n\n (commits: 10)"


def test_contributor_build_with_name() -> None:
    """
    A contributor with a name displays the name instead of login.
    """
    contributor = Contributor(
        "dgarcia360", "http://#", contributions=5, name="David Garcia"
    ).build(
        _CLASS_NAME,
        { 'show_contributions' : True }
    )
    assert contributor.astext() == "David Garcia\n\nCode:\n\n (commits: 5)"


def test_contributor_build_with_empty_name() -> None:
    """
    A contributor with an empty name falls back to login.
    """
    contributor = Contributor("dgarcia360", "http://#", contributions=5, name="").build(
        _CLASS_NAME,
        { 'show_contributions' : True }
    )
    assert contributor.astext() == "dgarcia360\n\nCode:\n\n (commits: 5)"


def test_contributor_build_with_no_contribution() -> None:
    """
    Create a contributor with no contributions.
    """
    contributor = Contributor("dgarcia360", "http://#").build(
        _CLASS_NAME,
        { 'show_contributions' : True }
    )
    assert contributor.astext() == "dgarcia360"


def test_contributor_build_with_one_contribution() -> None:
    """
    Create a contributor with one contribution.
    """
    contributor = Contributor("dgarcia360", "http://#", contributions=1).build(
        _CLASS_NAME,
        { 'show_contributions' : True }
    )
    assert contributor.astext() == "dgarcia360\n\nCode:\n\n (commits: 1)"


def test_contributor_repository_build() -> None:
    """
    Initializing a repository build.
    """
    contributors = [
        Contributor("dgarcia360", "http://#", contributions=2),
        Contributor("user", "http://#", contributions=1),
    ]
    contributor_repository = ContributorsRepository(
        contributors,
        {
          'reverse' : True,
          'show_contributions' : True
        }
    ).build(
        _CLASS_NAME
    )
    assert (
        contributor_repository.astext()
        == "dgarcia360\n\nCode:\n\n (commits: 2)\n\nuser\n\nCode:\n\n (commits: 1)"
    )


def test_contributor_repository_build_empty() -> None:
    """
    Initializing a contributor repository with no contributors.
    """
    contributors = []
    contributor_repository = ContributorsRepository(contributors, { 'reverse' : True }).build(
        _CLASS_NAME
    )
    assert contributor_repository.astext() == ""


def test_contributor_repository_build_order_desc() -> None:
    """
    The sort option (descending).
    """
    contributors = [
        Contributor("user", "http://#", contributions=1),
        Contributor("dgarcia360", "http://#", contributions=2),
    ]
    contributor_repository = ContributorsRepository(
        contributors,
        {
          'sort_by' : 'commit',
          'reverse' : True,
          'show_contributions' : True
        }
    ).build(
        _CLASS_NAME
    )
    assert (
        contributor_repository.astext()
        == "dgarcia360\n\nCode:\n\n (commits: 2)\n\nuser\n\nCode:\n\n (commits: 1)"
    )


def test_contributor_repository_build_order_asc() -> None:
    """
    The sort option works (ascending).
    """
    contributors = [
        Contributor("dgarcia360", "http://#", contributions=2),
        Contributor("user", "http://#", 1),
    ]
    contributor_repository = ContributorsRepository(
        contributors,
        {
          'sort_by' : 'commit',
          'reverse' : False,
          'show_contributions' : True
        }
    ).build(
        _CLASS_NAME
    )
    assert (
        contributor_repository.astext()
        == "user\n\nCode:\n\n (commits: 1)\n\ndgarcia360\n\nCode:\n\n (commits: 2)"
    )


def test_contributor_repository_build_with_limit() -> None:
    """
    The limit option works.
    """
    contributors = [
        Contributor("dgarcia360", "http://#", contributions=2),
        Contributor("user", "http://#", contributions=1),
    ]
    contributor_repository = ContributorsRepository(
        contributors,
        {
          'reverse' : True,
          'limit' : 1,
          'show_contributions' : True
        }
    ).build(
        _CLASS_NAME
    )
    assert contributor_repository.astext() == "dgarcia360\n\nCode:\n\n (commits: 2)"


def test_contributor_repository_build_exclude() -> None:
    """
    The exclude option works.
    """
    contributors = [
        Contributor("dgarcia360", "http://#", contributions=2),
        Contributor("sphinx", "http://#", contributions=1),
        Contributor("user", "http://#", contributions=1),
    ]
    exclude = [ "sphinx", "user" ]
    contributor_repository = ContributorsRepository(
        contributors,
        {
          'reverse' : True,
          'limit' : 10,
          'show_contributions' : True,
          'exclude' : exclude
        }
    ).build(
        _CLASS_NAME
    )
    assert contributor_repository.astext() == "dgarcia360\n\nCode:\n\n (commits: 2)"


@patch("sphinx_contributors.repository.requests.get")
def test_multiple_repositories(mock_get) -> None:
    """
    Multiple repositories merge contributors and sum contributions.
    """
    from sphinx_contributors import ContributorsDirective

    mock_repo1 = MagicMock()
    mock_repo1.json.return_value = [
        {
            "login": "shared",
            "html_url": "https://github.com/shared",
            "contributions": 3,
            "avatar_url": "",
            "type": "User"
        },
        {
            "login": "only_repo1",
            "html_url": "https://github.com/only_repo1",
            "contributions": 1,
            "avatar_url": "",
            "type": "User",
        },
    ]
    mock_repo1.links = {}
    mock_repo2 = MagicMock()
    mock_repo2.json.return_value = [
        {
            "login": "shared",
            "html_url": "https://github.com/shared",
            "contributions": 5,
            "avatar_url": "",
            "type": "User",
        },
        {
            "login": "only_repo2",
            "html_url": "https://github.com/only_repo2",
            "contributions": 2,
            "avatar_url": "",
            "type": "User",
        },
    ]
    mock_repo2.links = {}
    mock_get.side_effect = [
        mock_repo1,
        mock_repo2,
    ]

    directive = ContributorsDirective.__new__(ContributorsDirective)
    directive.arguments = ["owner/repo1 owner/repo2"]
    directive.options = {"contributions": None}
    directive.content = []
    directive.lineno = 0
    directive.state = MagicMock()
    directive.state_machine = MagicMock()

    result = directive.run()
    text = result[0].astext()
    mock_get.assert_called()
    assert "shared" in text
    assert "commits: 8" in text
    assert "only_repo1" in text
    assert "only_repo2" in text


@patch("sphinx_contributors.repository.requests.get")
def test_single_repository_unchanged(mock_get) -> None:
    """
    A single repository still works as before.
    """
    from sphinx_contributors import ContributorsDirective

    mock_response = MagicMock()
    mock_response.json.return_value = [
        {
            "login": "user1",
            "html_url": "https://github.com/user1",
            "contributions": 10,
            "avatar_url": "",
            "type": "User",
        },
    ]
    mock_response.links = {}
    mock_get.side_effect = [
        mock_response,
    ]

    directive = ContributorsDirective.__new__(ContributorsDirective)
    directive.arguments = ["owner/repo"]
    directive.options = {}
    directive.content = []
    directive.lineno = 0
    directive.state = MagicMock()
    directive.state_machine = MagicMock()

    result = directive.run()
    text = result[0].astext()
    assert "user1" in text
    mock_get.assert_called()


@patch("sphinx_contributors.repository.requests.get")
def test_include_adds_new_contributors(mock_get) -> None:
    """
    The include option adds contributors not already in the list.
    """
    from sphinx_contributors import ContributorsDirective

    mock_repo_response = MagicMock()
    mock_repo_response.json.return_value = [
        {
            "login": "existing",
            "html_url": "https://github.com/existing",
            "contributions": 5,
            "avatar_url": "",
            "type": "User",
        },
    ]
    mock_repo_response.links = {}
    mock_user_response = MagicMock()
    mock_user_response.json.return_value = {
        "html_url": "https://github.com/newuser",
        "avatar_url": "https://github.com/newuser.png",
    }
    mock_get.side_effect = [
        mock_repo_response,
        mock_user_response
    ]

    directive = ContributorsDirective.__new__(ContributorsDirective)
    directive.arguments = ["owner/repo"]
    directive.options = {"include": "newuser"}
    directive.content = []
    directive.lineno = 0
    directive.state = MagicMock()
    directive.state_machine = MagicMock()

    result = directive.run()
    text = result[0].astext()
    assert "existing" in text
    assert "newuser" in text


@patch("sphinx_contributors.repository.requests.get")
def test_include_does_not_duplicate_existing(mock_get) -> None:
    """
    The include option does not duplicate a contributor already returned by the API.
    """
    from sphinx_contributors import ContributorsDirective

    mock_repo_response = MagicMock()
    mock_repo_response.json.return_value = [
        {
            "login": "existing",
            "html_url": "https://github.com/existing",
            "contributions": 5,
            "avatar_url": "",
            "type": "User",
        },
    ]
    mock_repo_response.links = {}
    mock_get.return_value = mock_repo_response

    directive = ContributorsDirective.__new__(ContributorsDirective)
    directive.arguments = ["owner/repo"]
    directive.options = {"include": "existing"}
    directive.content = []
    directive.lineno = 0
    directive.state = MagicMock()
    directive.state_machine = MagicMock()

    result = directive.run()
    text = result[0].astext()
    assert text.count("existing") == 1


@patch("sphinx_contributors.repository.requests.get")
def test_contributor_directive(mock_get, tmp_path: Path) -> None:
    """
    The ``contributors`` directive runs with no errors.

    Uses a mocked GitHub API response so the test works without network access.
    """
    mock_response = MagicMock()
    mock_response.json.return_value = [
        {
            "login": "testuser",
            "html_url": "https://github.com/testuser",
            "contributions": 42,
            "avatar_url": "https://github.com/testuser.png",
            "type": "User",
        },
    ]
    mock_response.links = {}
    mock_get.return_value = mock_response

    source_directory = tmp_path / "source"
    source_directory.mkdir()
    source_file = source_directory / "index.rst"
    conf_py = source_directory / "conf.py"
    conf_py.touch()
    source_file.touch()
    conf_py_content = dedent(
        """\
        extensions = ['sphinx_contributors']
        """,
    )
    conf_py.write_text(conf_py_content)
    source_file_content = dedent(
        """\
        Test
        ====

        .. contributors:: sphinx-doc/sphinx
        """,
    )
    source_file.write_text(source_file_content)
    destination_directory = tmp_path / "destination"
    doctree_directory = tmp_path / "doctrees"

    app = Sphinx(
        srcdir=str(source_directory),
        confdir=str(source_directory),
        outdir=str(destination_directory),
        doctreedir=str(doctree_directory),
        buildername="html",
        warningiserror=True,
    )
    app.build()

    assert app.statuscode == 0
    mock_get.assert_called()


@patch("sphinx_contributors.repository._git_shortlog")
def test_contributor_directive_local(mock_git, tmp_path: Path) -> None:
    """
    The ``contributors`` directive runs with no errors.

    Uses a local repository 
    """
    mock_git.return_value = "  4584\tRonny Lorenz <ronny@tbi.univie.ac.at>"

    source_directory = tmp_path / "source"
    source_directory.mkdir()
    source_file = source_directory / "index.rst"
    conf_py = source_directory / "conf.py"
    conf_py.touch()
    source_file.touch()
    conf_py_content = dedent(
        """\
        extensions = ['sphinx_contributors']
        """,
    )
    conf_py.write_text(conf_py_content)
    source_file_content = dedent(
        """\
        Test
        ====

        .. contributors:: .
            :provider: local
        """,
    )
    source_file.write_text(source_file_content)
    destination_directory = tmp_path / "destination"
    doctree_directory = tmp_path / "doctrees"

    app = Sphinx(
        srcdir=str(source_directory),
        confdir=str(source_directory),
        outdir=str(destination_directory),
        doctreedir=str(doctree_directory),
        buildername="html",
        warningiserror=False, # Can't instantiate Sphinx() twice without warnings, see https://github.com/sphinx-doc/sphinx/issues/13608
        freshenv = True
    )
    app.build(force_all = True)

    assert app.statuscode == 0
    mock_git.assert_called()
