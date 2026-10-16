"""
Contributors extension for Sphinx.
(c) 2018 - present David Garcia (@dgarcia360)
# This code is licensed under MIT license (see LICENSE.md for details)
"""

__version__ = "0.3.0"

import copy
import os
from os import path
from pathlib import Path
from docutils.parsers.rst import Directive, directives
from sphinx.util import logging, osutil

from .contributors import Contributor, ContributorsRepository
from .repository import Repository
from .all_contributors import AllContributors

logger = logging.getLogger(__name__)

DEFAULT_CONTRIBUTORS_CONF = {
    "avatar_image_default" : "sc_user_avartar_default.svg",
    "avatar_size_default" : 150,
    "all_contributors" : False,
    "all_contributors_file" : ".all-contributorsrc",
    "all_contributors_imageSize_ignore" : False,
    "all_contributors_precedence" : [ "name", "avatar_url", "email", "profile" ],
}


def sc_path_static():
    """Returns path to packaged static files"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "_static"))


def sc_path_static_append(app) :
    contributors_config = app.config.sphinx_contributors_conf

    if contributors_config['avatar_image_default'] == DEFAULT_CONTRIBUTORS_CONF['avatar_image_default']:
        contributors_config['avatar_image_default'] = Path(
                                                os.path.relpath(
                                                  os.path.join(sc_path_static(), contributors_config['avatar_image_default']),
                                                  app.builder.srcdir
                                                )
                                             ).as_posix()

    # manually register the default image so that it works with latex builder
    app.env.images.add_file('', contributors_config['avatar_image_default'])


def fill_contributors_conf_defaults(app, config, check_keys = True):
    """Handle user config"""
    contributors_conf = copy.deepcopy(DEFAULT_CONTRIBUTORS_CONF)
    contributors_conf.update(config.sphinx_contributors_conf)

    config.sphinx_contributors_conf = contributors_conf

    config.html_static_path.append(sc_path_static())


def parse_all_contributors(app, docname, source):
    contributors_config = app.config.sphinx_contributors_conf

    if contributors_config.get("all_contributors", False) and \
        "all_contributors_file" in contributors_config:
        ac = AllContributors(contributors_config["all_contributors_file"])
        contributors_config['all_contributors_data'] = ac

class ContributorsDirective(Directive):
    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {
        "avatars": directives.flag,
        "avatars_only": directives.flag,
        "class_name": directives.unchanged,
        "contributions": directives.unchanged,
        "exclude": directives.unchanged,
        "include": directives.unchanged,
        "limit": directives.positive_int,
        "names": directives.flag,
        "anonymous": directives.flag,
        "order": directives.unchanged,
        "provider": directives.unchanged,
        "avatar_size" : directives.positive_int,
    }

    def run(self):
        config = self.state.document.settings.env.config.sphinx_contributors_conf

        order = self.options.get("order", "DESC").split()
        reverse = False if "ASC" in order or 'asc' in order else True
        order = [ o for o in order if o not in ('ASC', 'DESC', 'asc', 'desc') ]

        if "byContribution" in order:
            sort_by = 'contribution'
        if "byCommit" in order:
            sort_by = 'commit'
        elif "byName" in order:
            sort_by = 'name'
        elif "byLogin" in order:
            sort_by = 'login'
        elif "byEmail" in order:
            sort_by = 'email'
        else:
            sort_by = None

        contributors_options = {
            # pass over and pre-process options taken from the directive call
            'provider' : self.options.get("provider", "github"),
            # compile list of additional users to exclude/exclude
            'exclude' : [
                _exclude.strip() for _exclude in self.options.get("exclude", "").split(",")
            ],
            'include' : [
                _include.strip()
                for _include in self.options.get("include", "").split(",")
                if _include.strip()
            ],
            'sort_by' : sort_by,
            'reverse' : reverse,
            'rst_dir' : os.path.dirname(self.state_machine.document.attributes['source'])
                        if self.state_machine.document.attributes['source'] else ".",
            'anonymous' :  "anonymous" in self.options,
            'show_names' : "names" in self.options,
            'show_contributions' : self.options.get('contributions', False),
            'avatars_only' : "avatars_only" in self.options,
            'show_avatars' : "avatars" in self.options or "avatars_only" in self.options,
            'default_avatar' : config.get('avatar_image_default', ""),
            # store path to invoking rst source file in addition
            # to the other options set by this call of the contributors
            # directive
            'rst_source' : self.state_machine.document.attributes['source'],
            'limit': self.options.get('limit', None),
            'avatar_size': self.options.get('avatar_size', 0)
        }

        contributors_by = {}

        if  config.get('all_contributors', False) and \
            'all_contributors_data' in config:

            if not config.get('all_contributors_imageSize_ignore', False) and \
               contributors_options['avatar_size'] == 0:
                contributors_options['avatar_size'] = config['all_contributors_data'].imageSize

            if False:
                key = 'email' if contributors_options['provider'] == 'local' else 'login'
            else:
                # always use login from all_contributors as key
                key = 'login'

            contributors_by = { c.get(key, f"__acid_{i}") : Contributor(
                                   login = c.get('login', ""),
                                   url = c.get('profile', ""),
                                   avatar_url = c.get('avatar_url',""),
                                   email = c.get('email',""),
                                   name = c.get('name',""),
                                   aliases = c.get('aliases', []),
                                   contributions = c.get('contributions', []),
                                   source = '.all-contributorsrc'
                                 ) for i, c in enumerate(config['all_contributors_data'].contributors)
                               }

            if config['all_contributors_data'].contributorsSortAlphabetically:
                contributors_options['sort_by'] = 'login'

        if contributors_options['avatar_size'] == 0:
            contributors_options['avatar_size'] = config['avatar_size_default']

        for r in self.arguments[0].split():
            repo = Repository(r, options = contributors_options)

            for k, v in repo.get_contributors().items():
                if k not in contributors_by:
                    contributors_by[k] = v
                else:
                    if '.all-contributorsrc' in contributors_by[k].source:
                        contributors_by[k].update(v, keep = config['all_contributors_precedence'])
                    else:
                        contributors_by[k].update(v)

        # attempt to merge entries
        merged = []
        for k, v in contributors_by.items():
            if k not in merged:
                aliases = v.aliases
                if v.email != "":
                    aliases += [ v.email ]
                for alias in aliases:
                    if alias in [ c for c in contributors_by.keys() if c != k ]:
                        if '.all-contributorsrc' in v.source:
                            v.update(contributors_by[alias], keep = config['all_contributors_precedence'])
                        else:
                            v.update(contributors_by[alias])

                        merged.append(alias)

        for m in merged:
            del contributors_by[m]

        contributors = list(contributors_by.values())

        repo = ContributorsRepository(
            contributors,
            options = contributors_options
        )

        class_name = self.options.get("class_name", "sphinx-contributors")

        return [repo.build(class_name)]


def setup(app):
    app.add_config_value('sphinx_contributors_conf', DEFAULT_CONTRIBUTORS_CONF, 'html')

    # Early filling of sphinx_contributors_conf defaults at config-inited
    app.connect("config-inited", fill_contributors_conf_defaults, priority=10)
    app.connect("builder-inited", sc_path_static_append)
    app.connect('env-before-read-docs', parse_all_contributors)

    # Add directive
    app.add_directive("contributors", ContributorsDirective)

    # Add CSS
    app.add_css_file("sphinx_contributors.css")

    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
