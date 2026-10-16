from docutils import nodes

class Contributor:
    def __init__(
                  self,
                  login,
                  url,
                  contributions = [],
                  avatar_url="",
                  name="",
                  email="",
                  aliases = [],
                  anonymous = True,
                  source = []
                ):
        self.login = login
        self.url = url
        self.avatar_url = avatar_url
        self.name = name
        self.email = email
        self.anonymous = anonymous
        self.aliases = []
        self.source = []

        from . import logger as log

        if aliases:
            if isinstance(aliases, str):
                self.aliases.append(aliases)
            elif isinstance(aliases, list) or isinstance(aliases, tuple):
                self.aliases.extend([a for a in aliases])

        if source:
            if isinstance(source, str):
                self.source.append(source)
            elif isinstance(source, list) or isinstance(source, tuple):
                self.source.extend([s for s in source])

        self.contributions = contributions
        if isinstance(self.contributions, int):
            # assume contribution is a code commit
            self.contributions = [ { 'type' : "code", 'commits': self.contributions } ]
        elif not isinstance(self.contributions, list):
            if isinstance(self.contributions, str) or isinstance(self.contributions, dict):
                self.contributions = [ self.contributions ]
            else:
                log.warning(f"contributions \"{self.contributions}\" of wrong type! Expected str, dict, or list")
                self.contributions = []

    @property
    def display_name(self):
        return self.name or self.login
    def _merge_contributions(self, data):
        for c in data:
            ctype = c if isinstance(c, str) else c['type']
            found = None
            for i, cc in enumerate(self.contributions):
                if isinstance(cc, str):
                    if cc == ctype:
                        found = i
                        break
                elif cc['type'] == ctype:
                    found = i
                    break

            if found != None:
                if isinstance(c, str):
                    continue
                if isinstance(self.contributions[found], str):
                    self.contributions[found] = c
                else:
                    # attempt to merge contributions
                    for k, v in c.items():
                        if k == 'type':
                            continue

                        if k in self.contributions[found]:
                            if isinstance(v, str):
                                self.contributions[found][k] = ", ".join(self.contributions[found][k], v)
                            else:
                                self.contributions[found][k] += v
                        else:
                            self.contributions[found][k] = v
            else:
                self.contributions.append(c)

    def update(self, data, keep = []):
        if isinstance(data, Contributor):
            self.source.extend([ s for s in data.source ])
            self.aliases.extend([ a for a in data.aliases ])
            if 'contributions' not in keep:
                self._merge_contributions(data.contributions)
            if data.url and 'profile' not in keep:
                self.url = data.url
            if  data.avatar_url and 'avatar_url' not in keep:
                self.avatar_url = data.avatar_url 
            if  data.name and 'name' not in keep:
                self.name = data.name
            if  data.email and 'email' not in keep:
                self.email = data.email
            if  data.login and 'login' not in keep:
                self.login = data.login
            if  data.anonymous and 'anonymous' not in keep:
                self.anonymous = data.anonymous


    def build(self, class_name, options = {}):
        container_class = class_name + "_contributor"
        image_class = container_class + "__image"
        username_class = container_class + "__username"
        contributions_class = container_class + "__contributions"
        show_login = options.get('show_login', False)
        show_contributions = options.get('show_contributions', False)

        # make a list out of the 'show_contributions' option
        if show_contributions == "" or show_contributions == None:
            show_contributions = [ "" ]
        if show_contributions == False:
            show_contributions = []
        elif isinstance(show_contributions, str):
            show_contributions = show_contributions.split()
        elif not isinstance(show_contributions, list):
            show_contributions = [ show_contributions ]

        node_container = nodes.container(classes=[container_class])

        if options.get('show_avatars', False):
            if self.avatar_url or 'default_avatar' in options:
                node_image = nodes.image(
                    uri=self.avatar_url if self.avatar_url else options['default_avatar'],
                    alt=self.display_name if not show_login else self.login,
                    classes=[image_class, "no-scaled-link"],
                )
                node_image['alt']     = f"{self.display_name if not show_login else self.login} avatar"
                node_image['width']   = f"{options.get('avatar_size', 100)}px"
                node_image['height']  = f"{options.get('avatar_size', 100)}px"

                if self.url:
                    node_image_link = nodes.reference("", refuri=self.url)
                    node_image_link += node_image
                    node_container += node_image_link
                else:
                    node_container += node_image

            if options.get('avatars_only', False):
                return node_container

        node_username = nodes.paragraph(classes=[username_class])

        login = self.login.replace('@', ' [at] ')  # replace @ by [at]

        if self.url or self.email:
            url = self.url if self.url else f"mailto:{self.email}"
            node_username += nodes.reference(text=self.display_name if not show_login else login,
                                         refuri=url)
        else:
            node_username += nodes.Text(self.display_name if not show_login else login)

        node_container += node_username

        if show_contributions and self.contributions:
            from .all_contributors import default_types as ac_types

            contributions = []

            if 'brief' in show_contributions:
                # show only the contribution types without any details
                contrib = []
                for c in self.contributions:
                    emoji = None
                    ctype = c if isinstance(c, str) else c['type']

                    known_contributions = ac_types()
                    if ctype in known_contributions:
                        emoji = known_contributions[ctype]['symbol']
                        ctype = known_contributions[ctype]['description']

                    text = ""
                    if emoji and [ e for e in show_contributions if e in ('emoji', 'emojis-only') ]:
                        text += emoji + " "

                    if emoji == None or 'emojis-only' not in show_contributions:
                        text += ctype

                    contrib.append(text)

                if [ e for e in show_contributions if e in ('emoji', 'emojis-only') ]:
                    contributions.append(nodes.Text(" ".join(contrib)))
                else:
                    contributions.append(nodes.Text(", ".join(contrib)))
            else: 
                # show detailed contributions
                if 'deflist' in show_contributions:
                    # create definition_list for each contribution
                    defs = []

                    for c in self.contributions:
                        emoji = None
                        ctype = c if isinstance(c, str) else c['type']

                        known_contributions = ac_types()
                        if ctype in known_contributions:
                            emoji = known_contributions[ctype]['symbol']
                            ctype = known_contributions[ctype]['description']
                        elif 'all' not in show_contributions:
                            # skip unknown / non-standard contributions if not requested otherwise
                            continue

                        items = []

                        if isinstance(c, str):
                            items.append(nodes.paragraph("", "", nodes.Text(" ")))
                        else:
                            if 'text' in c:
                                items.append(nodes.paragraph("", "", nodes.Text(str(c['text']))))
                            # collect items
                            for k, v in c.items():
                                if k in [ 'type', 'text' ]:
                                    continue
                                else:
                                    if not (isinstance(v, int) or isinstance(v, float)) or v > 0:
                                        items.append(nodes.paragraph("", "", nodes.Text(str(k) + ": " + f"{v}")))

                        term = nodes.term("", nodes.Text(ctype))
                        if emoji and [ e for e in show_contributions if e in ('emoji', 'emojis-only') ]:
                            term = nodes.term("", nodes.Text(emoji + " " + ctype))

                        definition = nodes.definition("", *items)
                        df = nodes.definition_list_item("", term, definition)
                        defs.append(df)

                    if defs:
                        contributions.append(nodes.definition_list("", *defs))
                else:
                    # create a bullet_list instead
                    items = []
                    for c in self.contributions:
                        emoji = None
                        ctype = c if isinstance(c, str) else c['type']

                        known_contributions = ac_types()
                        if ctype in known_contributions:
                            emoji = known_contributions[ctype]['symbol']
                            ctype = known_contributions[ctype]['description']
                        elif 'all' not in show_contributions:
                            # skip unknown / non-standard contributions if not requested otherwise
                            continue

                        details = ""
                        text = ""
                        if isinstance(c, dict):
                            text = c.get('text', "")

                            # collect further key value pairs
                            for k, v in c.items():
                                if k in [ 'type', 'text' ]:
                                    continue
                                else:
                                    if not (isinstance(v, int) or isinstance(v, float)) or v > 0:
                                        if details:
                                            details = details + ", " + str(k) + ": " + str(v)
                                        else:
                                            details = str(k) + ": " + str(v)

                        title = nodes.strong("", nodes.Text(f"{ctype}{':' if text or details else ''}"))

                        children = []

                        if emoji and [ e for e in show_contributions if e in ('emoji', 'emojis-only') ]:
                            children.append(nodes.Text(emoji + " "))

                        if emoji == None or 'emojis-only' not in show_contributions:
                            children.append(title)

                        if text:
                            children.append(nodes.Text(" " + text))
                        if details:
                            children.append(nodes.emphasis("", " (" + details + ")"))

                        items.append(nodes.list_item("", *children))

                    if items:
                        contributions.append(nodes.bullet_list("", *items))

            node_container += nodes.container("", *contributions, classes=[contributions_class])

        return node_container


def _get_commits(contributor):
    for c in contributor.contributions:
        if isinstance(c, dict) and (c.get('type', None) == 'code'):
            return c.get('commits', 0)

    return 0

class ContributorsRepository:
    def __init__(
        self,
        contributors,
        options = {}
    ):

        sort_by = options.get('sort_by', None)
        reverse = options.get('reverse', False)
        exclude = options.get('exclude', [])
        limit = options.get('limit', None)

        sorted_contributors = [c for c in contributors if c.login not in exclude]

        if sort_by in ['login', 'name', 'email']:
            sorted_contributors = sorted(
                sorted_contributors,
                key = attrgetter(sort_by),
                reverse=reverse,
            )
        elif sort_by == 'contribution':
            sorted_contributors = sorted(
                sorted_contributors,
                key = lambda c: len(c.contributions),
                reverse=reverse,
            )
        elif sort_by == 'commit':
            sorted_contributors = sorted(
                sorted_contributors,
                key=_get_commits,
                reverse=reverse,
            )

        self._contributors = (
            sorted_contributors[:limit] if limit else sorted_contributors
        )
        self.options = options


    def build(self, class_name):
        list_class = class_name + "_list"
        item_class = list_class + "__item"

        node_container = nodes.container(classes=[class_name])

        if self.options.get('show_avatars', False):
            node_container["classes"].append(class_name + "--avatars")

        node_list = nodes.bullet_list(classes=[list_class])

        for contributor in self._contributors:
            node_item = nodes.list_item(classes=[item_class])
            node_item += contributor.build(class_name, self.options)
            node_list += node_item
        node_container += node_list

        return node_container


