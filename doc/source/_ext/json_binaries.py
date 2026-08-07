import json
import os, re

from docutils import nodes
from sphinx.util.docutils import SphinxDirective


def detect_version_number(string, default = ""):
    try:
        return re.search(r"([\d.]+\w+)", string).group(1)
    except:
        return default


def detect_release_number(string, default = ""):
    try:
        return str(re.search(r"([\d.]+)\-([\d.]+)([\.\_\-])\w+", string).group(2))
    except:
        return str(default)


def detect_architecture(filename):
    architecture = "NA"
    archs = ['64bit', 'x86_64', 'amd64', '32bit', 'i686', 'i586', 'i486', 'i386', 'universal', 'arm64', 'src']
    # detect true architecture
    for a in archs:
        if re.search(a, filename):
            architecture = a
            break

    return architecture


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

        if not hasattr(self.env, 'jb_json_data'):
            self.env.jb_json_data = []

        for d in self.env.jb_json_data:
            print(d)

        return [nodes.paragraph("", f"Hello World. {self.options['dist']} - {self.options['dist-version']}")]


def parse_json_files(app, docname, source):
    if app.config.jb_json_dir:
        print(f"parsing json files in {app.config.jb_json_dir}")
        env = app.env

        if not hasattr(env, 'jb_json_data'):
            env.jb_json_data = []

        if not hasattr(env, 'jb_packages'):
            env.jb_packages = {}

        for flavor, p in app.config.jb_packages.items():
            for pp, desc in p.items():
                # attempt to open and parse corresponding json data
                filename = f"{app.config.jb_json_dir}/{pp}.json"
                try:
                    with open(filename, "r") as f:
                        d = json.load(f)
                        for k, v in d.items():
                            dist = k
                            if flavor.lower() == v['flavor'].lower():
                                for package, fname in v['package'].items():
                                    version = detect_version_number(package)
                                    release = detect_release_number(package)
                                    architecture = detect_architecture(package)
                                    if flavor not in env.jb_packages:
                                        env.jb_packages[flavor] = {}
                                    if dist not in env.jb_packages[flavor]:
                                        env.jb_packages[flavor][dist] = []
                                    
                                    fname = fname.split("/")[-1]
                                    env.jb_packages[flavor][dist].append(
                                    {
                                      'version': version,
                                      'release': release,
                                      'architecture': architecture,
                                      'filename': fname,
                                      'description': desc 
                                    })


                except:
                    print(f"Can't access json file \"{filename}\"")
                    pass

        print(env.jb_packages)


def setup(app: object) -> dict:
    app.add_config_value('jb_json_dir', '', 'html')
    app.add_config_value('jb_packages', {}, 'html')
    app.add_config_value('jb_archs', [], 'html')
    app.add_config_value('jb_latest_version', '', 'html')
    
    app.add_directive("binary-packages", ShowBinaryPackagesDirective)
    app.connect('env-before-read-docs', parse_json_files)
