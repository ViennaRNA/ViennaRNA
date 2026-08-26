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

              rst += self._get_rst_flavor(flavor, label, data, standalone = True, indent = 0)
              rst += self._add_flavor_extras(flavor, standalone = True, indent = 0)

              rst += "\n"
          else:
              print(f"flavor \"{flavor}\" not in available packages")
        else:
            rst +=  ".. tab-set::\n"
            rst +=  "    :class: sd-align-major-justify\n\n"
            for flavor, data in self.env.jb_packages.items():
                if flavor in ['sourcecode']:
                    continue

                label = self.env.jb_package_labels.get(flavor, flavor)
                rst += self._get_rst_flavor(flavor, label, data)
                rst += self._add_flavor_extras(flavor)

                rst += "\n"

        if self.env.jb_debug:
            print(rst)

        return self.parse_text_to_nodes(rst)


    def _add_flavor_extras(self, flavor, standalone = False, indent = 4):
        rst = ""
        for d in self.env.jb_package_text.get(flavor, []):
            for atype, text in d.items():
                if atype in ["warning", "note", "tip"]:
                    rst += self._format_admonition(atype, text, standalone, indent) + "\n"

        return rst


    def _format_admonition(self, atype, tip, standalone = False, indent = 4):
        mindent = " " * indent
        lindent = ""
        rst = []

        if not standalone:
            lindent = "    "

        rst.append(lindent + f".. {atype}::\n")
        rst.append(lindent + "    " + tip)

        text = ""
        for line in rst:
            text += mindent + line + "\n"

        return text


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
            flavor_v_text = flavor_v.replace("-", " ").replace("_", " ")
            flavor_v_link = flavor_v.replace("-", "_").replace(".", "_").lower()
            for v in d:
                version_string = f"v{v['version']}"
                if v['release'] != "":
                    version_string += f"-{v['release']}"
#                rst.append(lindent + f":octicon:`download;1em;sd-text-danger sd-font-weight-bold` " \
#                       f"`{v['description']} <https://www.tbi.univie.ac.at/RNA/download/{flavor}/{flavor_v_link}/{v['filename']}>`__ - " \
#                       f":bdg-primary:`{v['version']}`" \
#                       f":bdg-info:`{v['release'] if v['release'] != "" else '1'}`" \
#                       f":bdg-danger:`{v['architecture']}`" \
#                       f":bdg-secondary:`{flavor_v_text}`\n")
                rst.append(lindent + f".. button-link:: https://www.tbi.univie.ac.at/RNA/download/{flavor}/{flavor_v_link}/{v['filename']}")
                rst.append(lindent + f"    :color: primary")
                rst.append(lindent + f"    :outline:")
                rst.append(lindent + f"    :tooltip: Download {v['description']} - v{v['version']}-{v['release'] if v['release'] != "" else '1'} for {flavor_v_text} and {v['architecture']} architecture")
                rst.append(lindent + f"    :expand:\n")
                rst.append(lindent + f"    :octicon:`download;1em;sd-text-danger sd-font-weight-bold` " \
                                     f"{v['description']} " \
                                     f":bdg-info:`{version_string}`" \
                                     f":bdg-danger:`{v['architecture']}`\n")
        else:
            rst.append(lindent + ".. tab-set::")
            rst.append(lindent + "    :class: sd-align-major-start\n")
            for flavor_v, d in data.items():
                flavor_v_text = flavor_v.replace("-", " ").replace("_", " ")
                flavor_v_link = flavor_v.replace("-", "_").replace(".", "_").lower()
                rst.append(lindent + f"    .. tab-item:: {flavor_v_text}")
                rst.append(lindent + f"        :class-label: sd-text-center\n")
                for v in d:
                    version_string = f"v{v['version']}"
                    if v['release'] != "":
                        version_string += f"-{v['release']}"
#                    rst.append(lindent + f"        :octicon:`download;1em;sd-text-danger sd-font-weight-bold` " \
#                           f"`{v['description']} <https://www.tbi.univie.ac.at/RNA/download/{flavor}/{flavor_v_link}/{v['filename']}>`__ - " \
#                           f":bdg-primary:`{v['version']}`" \
#                           f":bdg-info:`{v['release'] if v['release'] != "" else '1'}`" \
#                           f":bdg-danger:`{v['architecture']} ; Target Architecture`" \
#                           f":bdg-secondary:`{flavor_v_text}`\n")
                    rst.append(lindent + f"        .. button-link:: https://www.tbi.univie.ac.at/RNA/download/{flavor}/{flavor_v_link}/{v['filename']}")
                    rst.append(lindent + f"            :color: primary")
                    rst.append(lindent + f"            :outline:")
                    rst.append(lindent + f"            :tooltip: Download {v['description']} - v{v['version']}-{v['release'] if v['release'] != "" else '1'} for {flavor_v_text} and {v['architecture']} architecture")
                    rst.append(lindent + f"            :expand:\n")
                    rst.append(lindent + f"            :octicon:`download;1em;sd-text-danger sd-font-weight-bold` " \
                                         f"{v['description']} " \
                                         f":bdg-info:`{version_string}`" \
                                         f":bdg-danger:`{v['architecture']}`\n")

        text = ""
        for line in rst:
            text += mindent + line + "\n"

        return text


def parse_json_files(app, docname, source):
    if app.config.jb_json_dir:
        env = app.env

        env.jb_packages = {}
        env.jb_package_labels = {}
        env.jb_package_text = {}
        env.jb_debug = app.config.jb_debug

        if env.jb_debug:
            print(f"parsing json files in {app.config.jb_json_dir}")

        for flavor, package in app.config.jb_packages.items():
            env.jb_package_labels[flavor] = package['label']
            if 'text' in package:
                env.jb_package_text[flavor] = package['text']

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

        if env.jb_debug:
            print(env.jb_packages)


def setup(app: object) -> dict:
    app.add_config_value('jb_json_dir', '', 'html')
    app.add_config_value('jb_packages', {}, 'html')
    app.add_config_value('jb_latest_version', '', 'html')
    app.add_config_value('jb_debug', False, 'html')

    app.add_directive("binary-packages", ShowBinaryPackagesDirective)
    app.connect('env-before-read-docs', parse_json_files)
