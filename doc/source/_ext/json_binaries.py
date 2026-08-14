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
      'flavor': str,
      'flavor-version': str
    }

    def run(self) -> list:
        flavor          = self.options.get('flavor', '')
        flavor_version  = self.options.get('flavor-version', '')

        if not hasattr(self.env, 'jb_packages'):
            self.env.jb_packages = {}

        print("packages available:")
        # construct RST code for the packages
        rst = ""

        if (flavor != ''):
          if flavor in self.env.jb_packages:
              label = self.env.jb_package_labels.get(flavor, flavor)
              data = {}

              if flavor_version != '':
                  # search and collect all entries for requested version
                  for flavor_v, d in self.env.jb_packages[flavor].items():
                      for dd in d:
                          if dd['version'] == flavor_version:
                              if flavor_v not in data:
                                  data[flavor_v] = []

                              data[flavor_v].append(dd)
              else:
                  data = self.env.jb_packages[flavor]

              rst += self._get_rst_flavor(flavor, label, data, standalone = True, indent = 0) + "\n"
          else:
              print(f"flavor \"{flavor}\" not in available packages")
        else:
            rst +=  ".. tab-set::\n"
            rst +=  "    :class: sd-align-major-justify\n\n"
            for flavor, data in self.env.jb_packages.items():
                if flavor in ['sourcecode']:
                    continue

                label = self.env.jb_package_labels.get(flavor, flavor)
                rst += self._get_rst_flavor(flavor, label, data) + "\n"

        print(rst)

        return self.parse_text_to_nodes(rst)

    def _get_rst_flavor(self, flavor, label, data, standalone = False, indent = 4) -> str:
        if not label or label == "":
            label = flavor
        mindent = " " * indent
        lindent = ""

        rst = []
        if not standalone:
            rst.append(f".. tab-item:: |{flavor}-logo| |br| :sub:`{label}`")
            rst.append(f"    :class-content: dl-{flavor}")
            rst.append("    :class-label: sd-text-center\n")
            lindent = "    "

        if len(data) == 1:
            flavor_v, d = next(iter( data.items() ))
            for v in d:
                rst.append(lindent + f":octicon:`download;1em;sd-text-danger sd-font-weight-bold` " \
                       f"`{v['description']} <https://www.tbi.univie.ac.at/RNA/download/{flavor}/{flavor_v}/{v['filename']}>`__ " \
                       f":bdg-primary:`{v['version']}`" \
                       f":bdg-info:`{v['release'] if v['release'] != "" else '1'}`" \
                       f":bdg-danger:`{v['architecture']} ; Target Architecture`" \
                       f":bdg-secondary:`{flavor_v}`\n")
        else:
            rst.append(lindent + ".. tab-set::")
            rst.append(lindent + "    :class: sd-align-major-start\n")
            for flavor_v, d in data.items():
                rst.append(lindent + f"    .. tab-item:: {flavor_v}")
                rst.append(lindent + f"        :class-label: sd-text-center\n")
                for v in d:
#                    rst.append(lindent * 2 + f"    .. div:: \n")
                    rst.append(lindent + f"        :octicon:`download;1em;sd-text-danger sd-font-weight-bold` " \
                           f"`{v['description']} <https://www.tbi.univie.ac.at/RNA/download/{flavor}/{flavor_v}/{v['filename']}>`__ " \
                           f":bdg-primary:`{v['version']}`" \
                           f":bdg-info:`{v['release'] if v['release'] != "" else '1'}`" \
                           f":bdg-danger:`{v['architecture']} ; Target Architecture`" \
                           f":bdg-secondary:`{flavor_v}`\n")

        text = ""
        for line in rst:
            text += mindent + line + "\n"

        return text


def parse_json_files(app, docname, source):
    if app.config.jb_json_dir:
        print(f"parsing json files in {app.config.jb_json_dir}")
        env = app.env

        env.jb_packages = {}
        env.jb_package_labels = {}

        for flavor, package in app.config.jb_packages.items():
            env.jb_package_labels[flavor] = package['label']

            for pp, desc in package['package-items'].items():
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
    app.add_config_value('jb_latest_version', '', 'html')
    
    app.add_directive("binary-packages", ShowBinaryPackagesDirective)
    app.connect('env-before-read-docs', parse_json_files)
