from docutils import nodes
from sphinx.util.docutils import SphinxDirective


class ShowBinaryPackagesDirective(SphinxDirective):
    has_content = True
    optional_arguments = 2
    option_spec = {
      'dist': str,
      'dist-version': str
    }

    def run(self) -> list:
        self.options['dist'] = self.options.get('dist', '')
        self.options['dist-version'] = self.options.get('dist-version', '')

        return [nodes.paragraph("", f"Hello World. {self.options['dist']} - {self.options['dist-version']}")]


def setup(app: object) -> dict:
    app.add_directive("binary-packages", ShowBinaryPackagesDirective)
