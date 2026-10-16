import json


def default_types(repo = None):
    return {
      "a11y":{
        "symbol": '️️️️♿️',
        "description": 'Accessibility',
      },
      "audio":{
        "symbol": '🔊',
        "description": 'Audio',
      },
      "blog":{
        "symbol": '📝',
        "description": 'Blogposts',
      },
      "bug":{
        "symbol": '🐛',
        "description": 'Bug reports',
        "link": repo.issue_link if repo else None,
      },
      "business":{
        "symbol": '💼',
        "description": 'Business development',
      },
      "code":{
        "symbol": '💻',
        "description": 'Code',
        "link": repo.commits_link if repo else None,
      },
      "content":{
        "symbol": '🖋',
        "description": 'Content',
      },
      "data":{
        "symbol": '🔣',
        "description": 'Data',
      },
      "design":{
        "symbol": '🎨',
        "description": 'Design',
      },
      "doc":{
        "symbol": '📖',
        "description": 'Documentation',
        "link": repo.commits_link if repo else None,
      },
      "eventOrganizing":{
        "symbol": '📋',
        "description": 'Event Organizing',
      },
      "example":{
        "symbol": '💡',
        "description": 'Examples',
      },
      "financial":{
        "symbol": '💵',
        "description": 'Financial',
      },
      "fundingFinding":{
        "symbol": '🔍',
        "description": 'Funding Finding',
      },
      "ideas":{
        "symbol": '🤔',
        "description": 'Ideas, Planning, & Feedback',
      },
      "infra":{
        "symbol": '🚇',
        "description": 'Infrastructure (Hosting, Build-Tools, etc)',
      },
      "maintenance":{
        "symbol": '🚧',
        "description": 'Maintenance',
      },
      "mentoring":{
        "symbol": '🧑‍🏫',
        "description": 'Mentoring',
      },
      "platform":{
        "symbol": '📦',
        "description": 'Packaging/porting to new platform',
      },
      "plugin":{
        "symbol": '🔌',
        "description": 'Plugin/utility libraries',
      },
      "projectManagement":{
        "symbol": '📆',
        "description": 'Project Management',
      },
      "question":{
        "symbol": '💬',
        "description": 'Answering Questions',
      },
      "research":{
        "symbol": '🔬',
        "description": 'Research',
      },
      "review":{
        "symbol": '👀',
        "description": 'Reviewed Pull Requests',
        "link": repo.revoews_link if repo else None,
      },
      "security":{
        "symbol": '🛡️',
        "description": 'Security',
      },
      "talk":{
        "symbol": '📢',
        "description": 'Talks',
      },
      "test":{
        "symbol": '⚠️',
        "description": 'Tests',
        "link": repo.commits_link if repo else None,
      },
      "tool":{
        "symbol": '🔧',
        "description": 'Tools',
      },
      "translation":{
        "symbol": '🌍',
        "description": 'Translation',
      },
      "tutorial":{
        "symbol": '✅',
        "description": 'Tutorials',
      },
      "userTesting":{
        "symbol": '📓',
        "description": 'User Testing',
      },
      "video":{
        "symbol": '📹',
        "description": 'Videos',
      },
      "promotion":{
        "symbol": '📣',
        "description": 'Promotion',
      }
    }


def _parse_allcontributors_json(filename):

    from . import logger as log

    try:
        with open(filename, "r", encoding='utf-8') as f:
            return json.load(f)

    except Exception as err:
        log.warning(f"Can't access all-contributors json file \"{filename}\" {err=}, {type(err)=}")
        pass

    return {}


class AllContributors:
    def __init__(self, config_file):
        self.data = _parse_allcontributors_json(config_file)

    @property
    def contributorsPerLine(self):
        return self.data.get('contributorsPerLine', 7)

    @property
    def imageSize(self):
        return int(self.data['imageSize']) if 'imageSize' in self.data and int(self.data['imageSize']) > 0 else 0

    @property
    def contributorsSortAlphabetically(self):
        return self.data.get('contributorsSortAlphabetically', False)

    @property
    def contributors(self):
        return self.data['contributors'] if 'contributors' in self.data else []

