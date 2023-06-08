#!/usr/bin/env python
"""doxy2swig.py [options] index.xml output.i

Doxygen XML to SWIG docstring converter (improved version).

Converts Doxygen generated XML files into a file containing docstrings
for use by SWIG.

index.xml is your doxygen generated XML file and output.i is where the
output will be written (the file will be clobbered).
"""
#
# The current version of this code is hosted on a github repository:
#   https://github.com/m7thon/doxy2swig
#
# This code is implemented using Mark Pilgrim's code as a guideline:
#   http://www.faqs.org/docs/diveintopython/kgp_divein.html
#
# Original Author: Prabhu Ramachandran
# Modified by:     Michael Thon (June 2015)
# Modified by:     Ronny Lorenz (May 2023)
# License: BSD style
#
# Thanks:
#   Johan Hake:  the include_function_definition feature
#   Bill Spotz:  bug reports and testing.
#   Sebastian Henschel:   Misc. enhancements.
#
# Changes:
# May 2023 (Ronny Lorenz):
#   - 
# June 2015 (Michael Thon):
#   - class documentation:
#     -c: add constructor call signatures and a "Constructors" section
#         collecting the respective docs (e.g. for python)
#     -a: add "Attributes" section collecting the documentation for member
#         variables (e.g. for python)
#   - overloaded functions:
#     -o: collect all documentation into one "Overloaded function" section
#   - option to include function definition / signature renamed to -f
#   - formatting:
#     + included function signatures slightly reformatted
#     + option (-t) to turn off/on type information for funciton signatures
#     + lists (incl. nested and ordered)
#     + attempt to produce docstrings that render nicely as markdown
#     + translate code, emphasis, bold, linebreak, hruler, blockquote,
#       verbatim, heading tags to markdown
#     + new text-wrapping and option -w to specify the text width
#

from xml.dom import minidom
import re
import textwrap
import sys
import os.path
import optparse
import io
import csv

def my_open_read(source):
    if hasattr(source, "read"):
        return source
    else:
        try:
            return open(source, encoding='utf-8')
        except TypeError:
            return open(source)

def my_open_write(dest):
    if hasattr(dest, "write"):
        return dest
    else:
        try:
            return open(dest, 'w', encoding='utf-8')
        except TypeError:
            return open(dest, 'w')

# MARK: Text handling:
def shift(txt, indent = '    ', prepend = ''):
    """Return a list corresponding to the lines of text in the `txt` list
    indented by `indent`. Prepend instead the string given in `prepend` to the
    beginning of the first line. Note that if len(prepend) > len(indent), then
    `prepend` will be truncated (doing better is tricky!). This preserves a 
    special '' entry at the end of `txt` (see `do_para` for the meaning).
    """
    if type(indent) is int:
        indent = indent * ' '
    special_end = txt[-1:] == ['']
    lines = ''.join(txt).splitlines(True)
    for i in range(1,len(lines)):
        if lines[i].strip() or indent.strip():
            lines[i] = indent + lines[i]
    if not lines:
        return prepend
    prepend = prepend[:len(indent)]
    indent = indent[len(prepend):]
    lines[0] = prepend + indent + lines[0]
    ret = [''.join(lines)]
    if special_end:
        ret.append('')
    return ret

class Doxy2SWIG:
    """Converts Doxygen generated XML files into a file containing
    docstrings that can be used by SWIG-1.3.x that have support for
    feature("docstring").  Once the data is parsed it is stored in
    self.pieces.

    """

    def __init__(self, src,
                 with_function_signature = False,
                 with_type_info = False,
                 with_constructor_list = False,
                 with_attribute_list = False,
                 with_overloaded_functions = False,
                 textwidth = 80,
                 quiet = False,
                 object_binding = dict(),
                 prefixes = None,
                 prefix_name = "",
                 suffixes = None,
                 constructor_suffixes = None,
                 typedef_mapping = dict(),
                 ignore_list = dict(),
                 deprecated_version = None):
        """Initialize the instance given a source object.  `src` can
        be a file or filename.  If you do not want to include function
        definitions from doxygen then set
        `include_function_definition` to `False`.  This is handy since
        this allows you to use the swig generated function definition
        using %feature("autodoc", [0,1]).

        """
        # options:
        self.with_function_signature = with_function_signature
        self.with_type_info = with_type_info
        self.with_constructor_list = with_constructor_list
        self.with_attribute_list = with_attribute_list
        self.with_overloaded_functions = with_overloaded_functions
        self.textwidth = textwidth
        self.quiet = quiet
        self.object_binding = object_binding
        self.prefixes = [] if not prefixes else prefixes
        self.constructor_suffixes = [] if not constructor_suffixes else constructor_suffixes
        self.prefix_name = prefix_name
        self.suffixes = [] if not suffixes else suffixes
        self.typedef_mapping = typedef_mapping
        self.ignore_list = ignore_list
        self.deprecated_version = '2.6.0' if not deprecated_version else deprecated_version

        # state:
        self.indent = 0
        self.listitem = ''
        self.pieces = []

        f = my_open_read(src)
        self.my_dir = os.path.dirname(f.name)
        self.xmldoc = minidom.parse(f).documentElement
        f.close()

        self.pieces.append('\n// File: %s\n' %
                           os.path.basename(f.name))

        self.space_re = re.compile(r'\s+')
        self.lead_spc = re.compile(r'^(%feature\S+\s+\S+\s*?)"\s+(\S)')
        self.math_block_re = re.compile(r'\\\[(.*?)\\\]', re.MULTILINE | re.DOTALL)
        self.math_inline_re = re.compile(r'\$(.*?)\$')
        self.since_version_re = re.compile(r'\(\s*[Ss]ince\s+v?(\d+\.\d+\.\d+)\s*\)')

        self.multi = 0
        self.ignores = ['inheritancegraph', 'param', 'listofallmembers',
                        'innerclass', 'name', 'declname', 'incdepgraph',
                        'invincdepgraph', 'programlisting', 'type',
                        'references', 'referencedby', 'location',
                        'collaborationgraph', 'reimplements',
                        'reimplementedby', 'derivedcompoundref',
                        'basecompoundref',
                        'argsstring', 'definition', 'exceptions']
        #self.generics = []

    def generate(self):
        """Parses the file set in the initialization.  The resulting
        data is stored in `self.pieces`.

        """
        self.enforce_numpydoc_order(self.xmldoc)
        self.parse(self.xmldoc)

    def write(self, fname):
        o = my_open_write(fname)
        o.write(''.join(self.pieces))
        o.write('\n')
        o.close()

    def make_last_child(self, nodes):
        for n in nodes:
            if (len(n.parentNode.childNodes) > 1):
                # remove node from its parent
                for i, nn in enumerate(n.parentNode.childNodes):
                    if n == nn:
                        n.parentNode.childNodes.pop(i)
                        break
                # append node as last in parent
                n.parentNode.childNodes.append(n)

    def enforce_numpydoc_order(self, node):
        """
        Numpydoc wants us to order function descriptions
        """
        if node.nodeType == node.ELEMENT_NODE and node.tagName == "memberdef" and node.getAttribute("kind") == "function":
            for nn in node.childNodes:
                if nn.nodeType == nn.ELEMENT_NODE and nn.tagName == "detaileddescription":
                    # move parameterlist to the end
                    nnn = self.get_specific_subnodes(nn, "parameterlist", 3)
                    self.make_last_child(nnn)

                    # move return values to the end
                    nnn = self.get_specific_subnodes(nn, "simplesect", 3, "return")
                    self.make_last_child(nnn)

                    # move warnings to the end
                    nnn = self.get_specific_subnodes(nn, "simplesect", 3, "warning")
                    self.make_last_child(nnn)

                    # move see-also to the end
                    nnn = self.get_specific_subnodes(nn, "simplesect", 3, "see")
                    self.make_last_child(nnn)

                    # move notes to the end
                    nnn = self.get_specific_subnodes(nn, "simplesect", 3, "note")
                    self.make_last_child(nnn)

        else:
            for n in node.childNodes:
                self.enforce_numpydoc_order(n)


    def parse(self, node):
        """Parse a given node.  This function in turn calls the
        `parse_<nodeType>` functions which handle the respective
        nodes.

        """
        pm = getattr(self, "parse_%s" % node.__class__.__name__)
        pm(node)

    def parse_Document(self, node):
        self.parse(node.documentElement)

    def parse_Text(self, node):
        txt = node.data
        if txt == ' ':
        #    # this can happen when two tags follow in a text, e.g.,
        #    # " ...</emph> <formaula>$..." etc.
        #    # here we want to keep the space.
        #    self.add_text(txt)
            return
        txt = txt.replace('\\', r'\\')
        txt = txt.replace('"', r'\"')

        # ignore pure whitespace
        m = self.space_re.match(txt)

        if not (m and len(m.group()) == len(txt)):
            self.add_text(txt)

    def parse_Comment(self, node):
        """Parse a `COMMENT_NODE`.  This does nothing for now."""
        return

    def parse_Element(self, node):
        """Parse an `ELEMENT_NODE`.  This calls specific
        `do_<tagName>` handers for different elements.  If no handler
        is available the `subnode_parse` method is called.  All
        tagNames specified in `self.ignores` are simply ignored.

        """
        name = node.tagName
        ignores = self.ignores
        if name in ignores:
            return
        attr = "do_%s" % name
        if hasattr(self, attr):
            handlerMethod = getattr(self, attr)
            handlerMethod(node)
        else:
            self.subnode_parse(node)
            #if name not in self.generics: self.generics.append(name)

# MARK: Special format parsing
    def subnode_parse(self, node, pieces=None, indent=0, ignore=[], restrict=None):
        """Parse the subnodes of a given node. Subnodes with tags in the
        `ignore` list are ignored. If pieces is given, use this as target for
        the parse results instead of self.pieces. Indent all lines by the amount
        given in `indent`. Note that the initial content in `pieces` is not
        indented. The final result is in any case added to self.pieces."""
        if pieces is not None:
            old_pieces, self.pieces = self.pieces, pieces
        else:
            old_pieces = []
        if type(indent) is int:
            indent = indent * ' '
        if len(indent) > 0:
            pieces = ''.join(self.pieces)
            i_piece = pieces[:len(indent)]
            if self.pieces[-1:] == ['']:
                self.pieces = [pieces[len(indent):]] + ['']
            elif self.pieces != []:
                self.pieces = [pieces[len(indent):]]
        self.indent += len(indent)
        for n in node.childNodes:
            if restrict is not None:
                if n.nodeType == n.ELEMENT_NODE and n.tagName in restrict:
                    self.parse(n)
            elif n.nodeType != n.ELEMENT_NODE or n.tagName not in ignore:
                self.parse(n)
        if len(indent) > 0:
            self.pieces = shift(self.pieces, indent, i_piece)
        self.indent -= len(indent)
        old_pieces.extend(self.pieces)
        self.pieces = old_pieces

    def surround_parse(self, node, pre_char, post_char):
        """Parse the subnodes of a given node. Subnodes with tags in the
        `ignore` list are ignored. Prepend `pre_char` and append `post_char` to
        the output in self.pieces."""
        self.add_text(pre_char)
        self.subnode_parse(node)
        self.add_text(post_char)
    
# MARK: Helper functions
    def get_specific_subnodes(self, node, name, recursive=0, kind = None):
        """Given a node and a name, return a list of child `ELEMENT_NODEs`, that
        have a `tagName` matching the `name`. Search recursively for `recursive`
        levels.
        """
        children = [x for x in node.childNodes if x.nodeType == x.ELEMENT_NODE]
        if kind:
            ret = [x for x in children if x.tagName == name and x.getAttribute("kind")  == kind]
        else:
            ret = [x for x in children if x.tagName == name]
        if recursive > 0:
            for x in children:
                ret.extend(self.get_specific_subnodes(x, name, recursive-1, kind))
        return ret

    def get_specific_parentnode(self, node, name, recursive = 0):
        if node.parentNode and \
           node.parentNode.nodeType == node.ELEMENT_NODE and \
           node.parentNode.tagName == name:
            return node.parentNode
        elif recursive > 1:
            return self.get_specific_parentnode(node.parentNode, name, recursive - 1)
        else:
            return None

    def get_specific_nodes(self, node, names):
        """Given a node and a sequence of strings in `names`, return a
        dictionary containing the names as keys and child
        `ELEMENT_NODEs`, that have a `tagName` equal to the name.

        """
        nodes = [(x.tagName, x) for x in node.childNodes
                 if x.nodeType == x.ELEMENT_NODE and
                 x.tagName in names]
        return dict(nodes)

    def add_text(self, value):
        """Adds text corresponding to `value` into `self.pieces`."""
        if isinstance(value, (list, tuple)):
            self.pieces.extend(value)
        else:
            self.pieces.append(value)

    def start_new_paragraph(self):
        """Make sure to create an empty line. This is overridden, if the previous
        text ends with the special marker ''. In that case, nothing is done.
        """
        if self.pieces[-1:] == ['']: # respect special marker
            return
        elif self.pieces == []: # first paragraph, add '\n', override with ''
            self.pieces = ['\n']
        elif self.pieces[-1][-1:] != '\n': # previous line not ended
            self.pieces.extend(['  \n' ,'\n'])
        else: #default
            self.pieces.append('\n')

    def add_line_with_subsequent_indent(self, line, indent=4):
        """Add line of text and wrap such that subsequent lines are indented
        by `indent` spaces.
        """
        if isinstance(line, (list, tuple)):
            line = ''.join(line)
        line = line.strip()
        width = self.textwidth-self.indent-indent
        wrapped_lines = textwrap.wrap(line[indent:], width=width)
        for i in range(len(wrapped_lines)):
            if wrapped_lines[i] != '':
                wrapped_lines[i] = indent * ' ' + wrapped_lines[i]
        self.pieces.append(line[:indent] + '\n'.join(wrapped_lines)[indent:] + '  \n')

    def extract_text(self, node):
        """Return the string representation of the node or list of nodes by parsing the
        subnodes, but returning the result as a string instead of adding it to `self.pieces`.
        Note that this allows extracting text even if the node is in the ignore list.
        """
        if not isinstance(node, (list, tuple)):
            node = [node]
        pieces, self.pieces = self.pieces, ['']
        for n in node:
            for sn in n.childNodes:
                self.parse(sn)
        ret = ''.join(self.pieces)
        self.pieces = pieces
        return ret

    def get_function_signature(self, node):
        """Returns the function signature string for memberdef nodes."""
        omit_param = None
        name = self.extract_text(self.get_specific_subnodes(node, 'name'))

        if name in self.object_binding:
            omit_param = self.object_binding[name]['self_arg_name']
            name = self.object_binding[name]['method']

        if self.with_type_info:
            argsstring = self.extract_text(self.get_specific_subnodes(node, 'argsstring'))
        else:
            argsstring = []
            param_id = 1
            for n_param in self.get_specific_subnodes(node, 'param'):
                declname = self.extract_text(self.get_specific_subnodes(n_param, 'declname'))
                if not declname:
                    declname = 'arg' + str(param_id)
                if declname == omit_param:
                    declname = 'self'

                defval = self.extract_text(self.get_specific_subnodes(n_param, 'defval'))
                if defval:
                    defval = '=' + defval
                argsstring.append(declname + defval)
                param_id = param_id + 1
            argsstring = '(' + ', '.join(argsstring) + ')'
        type = self.extract_text(self.get_specific_subnodes(node, 'type'))
        function_definition = name + argsstring
        if type != '' and type != 'void':
            function_definition = function_definition + ' -> ' + type
        return '`' + function_definition + '`  '

# MARK: Special parsing tasks (need to be called manually)
    def make_constructor_list(self, constructor_nodes, classname):
        """Produces the "Constructors" section and the constructor signatures
        (since swig does not do so for classes) for class docstrings."""
        if constructor_nodes == []:
            return
        self.add_text(['\n', 'Constructors',
                       '\n', '------------'])
        for n in constructor_nodes:
            self.add_text('\n')
            self.add_line_with_subsequent_indent('* ' + self.get_function_signature(n))
            self.subnode_parse(n, pieces = [], indent=4, ignore=['definition', 'name'])

    def make_attribute_list(self, node):
        """Produces the "Attributes" section in class docstrings for public
        member variables (attributes).
        """
        atr_nodes = []
        for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
            if n.attributes['kind'].value == 'variable' and n.attributes['prot'].value == 'public':
                atr_nodes.append(n)
        if not atr_nodes:
            return
        self.add_text(['\n', 'Attributes',
                       '\n', '----------'])
        for n in atr_nodes:
            name = self.extract_text(self.get_specific_subnodes(n, 'name'))
            if 1:
                self.add_text(['\n', name, ' : '])
                self.add_text([self.extract_text(self.get_specific_subnodes(n, 'type'))])
                self.add_text('  \n')
                restrict = ['briefdescription', 'detaileddescription']
                self.subnode_parse(n, pieces=[''], indent=4, restrict=restrict)
            else:
                self.add_text(['\n* ', '`', name, '`', ' : '])
                self.add_text(['`', self.extract_text(self.get_specific_subnodes(n, 'type')), '`'])
                self.add_text('  \n')
                restrict = ['briefdescription', 'detaileddescription']
                self.subnode_parse(n, pieces=[''], indent=4, restrict=restrict)

    def get_memberdef_nodes_and_signatures(self, node, kind):
        """Collects the memberdef nodes and corresponding signatures that
        correspond to public function entries that are at most depth 2 deeper
        than the current (compounddef) node. Returns a dictionary with 
        function signatures (what swig expects after the %feature directive)
        as keys, and a list of corresponding memberdef nodes as values."""
        sig_dict = {}
        sig_prefix = ''
        
        if kind in ('file'):
            # TODO: This is still not correct as files may have multiple innernamespace tags, but I don't know
            # if any information is actually extracted from <compounddef kind="file"... tags?
            ns_node = node.getElementsByTagName('innernamespace')
            if ns_node:
                sig_prefix = self.extract_text(ns_node[0]) + '::'
        elif kind in ('namespace'):
            # Namespace name is stored in compoundname tag.
            ns_node = node.getElementsByTagName('compoundname')
            if ns_node:
                sig_prefix = self.extract_text(ns_node[0]) + '::'
        elif kind in ('class', 'struct'):
            # Get the full function name.
            cn_node = node.getElementsByTagName('compoundname')
            sig_prefix = self.extract_text(cn_node[0]) + '::'

        md_nodes = self.get_specific_subnodes(node, 'memberdef', recursive=2)
        for n in md_nodes:
            if n.attributes['prot'].value != 'public':
                continue
            if n.attributes['kind'].value in ['variable', 'typedef']:
                continue
            if n.attributes['kind'].value not in ['define'] and not self.get_specific_subnodes(n, 'definition'):
                continue
            name = self.extract_text(self.get_specific_subnodes(n, 'name'))
            if name[:8] == 'operator':
                continue
            sig = sig_prefix + name
            if sig in sig_dict:
                sig_dict[sig].append(n)
            else:
                sig_dict[sig] = [n]
        return sig_dict
    
    def handle_typical_memberdefs_no_overload(self, signature, memberdef_nodes):
        """Produce standard documentation for memberdef_nodes."""
        for n in memberdef_nodes:
            if signature in self.object_binding:
                # prepend the class name if function is bound to object
                signature = self.object_binding[signature]['class'] + \
                            "::" + \
                            self.object_binding[signature]['method']
            elif signature in self.ignore_list:
                # rename signature if the swig interface defines a replacement
                if self.ignore_list[signature]['replacement']:
                    signature = self.ignore_list[signature]['replacement']
                else:
                    # skip signature if it should really be ignored
                    continue
            elif n.attributes['kind'].value == 'define' and len(self.prefixes) > 0:
                # automatically remove prefixes from #define cpp macros
                prefix = r''
                for p in self.prefixes:
                    signature = signature.replace(p, prefix)
            self.add_text(['\n', '%feature("docstring") ', signature, ' "', '\n'])
            if self.with_function_signature:
                self.add_line_with_subsequent_indent(self.get_function_signature(n))
            self.subnode_parse(n, pieces=[], ignore=['definition', 'name', 'initializer'])
            self.add_text(['";', '\n'])

    def handle_typical_memberdefs(self, signature, memberdef_nodes):
        """Produces docstring entries containing an "Overloaded function"
        section with the documentation for each overload, if the function is
        overloaded and self.with_overloaded_functions is set. Else, produce
        normal documentation.
        """
        if len(memberdef_nodes) == 1 or not self.with_overloaded_functions:
            self.handle_typical_memberdefs_no_overload(signature, memberdef_nodes)
            return
        self.add_text(['\n', '%feature("docstring") ', signature, ' "', '\n'])
        if self.with_function_signature:
            for n in memberdef_nodes:
                self.add_line_with_subsequent_indent(self.get_function_signature(n))
        self.add_text('\n')
        self.add_text(['Overloaded function', '\n',
                       '-------------------'])
        for n in memberdef_nodes:
            self.add_text('\n')
            self.add_line_with_subsequent_indent('* ' + self.get_function_signature(n))
            self.subnode_parse(n, pieces=[], indent=4, ignore=['definition', 'name'])
        self.add_text(['";', '\n'])
    

# MARK: Tag handlers
    def do_linebreak(self, node):
        self.add_text('  ')
    
    def do_ndash(self, node):
        self.add_text('--')

    def do_mdash(self, node):
        self.add_text('---')

    def do_emphasis(self, node):
        self.surround_parse(node, '*', '*')

    def do_bold(self, node):
        self.surround_parse(node, '**', '**')
    
    def do_computeroutput(self, node):
        self.surround_parse(node, '`', '`')

    def do_heading(self, node):
        self.start_new_paragraph()
        pieces, self.pieces = self.pieces, ['']
        level = int(node.attributes['level'].value)
        self.subnode_parse(node)
        if level == 1:
            self.pieces.insert(0, '\n')
            self.add_text(['\n', len(''.join(self.pieces).strip()) * '='])
        elif level == 2:
            self.add_text(['\n', len(''.join(self.pieces).strip()) * '-'])
        elif level >= 3:
            self.pieces.insert(0, level * '#' + ' ')
        # make following text have no gap to the heading:
        pieces.extend([''.join(self.pieces) + '  \n', ''])
        self.pieces = pieces
    
    def do_verbatim(self, node):
        self.start_new_paragraph()
        self.subnode_parse(node, pieces=[''], indent=4)
    
    def do_blockquote(self, node):
        self.start_new_paragraph()
        self.subnode_parse(node, pieces=[''], indent='> ')
    
    def do_hruler(self, node):
        self.start_new_paragraph()
        self.add_text('* * * * *  \n')
    
    def do_includes(self, node):
        self.add_text('\nC++ includes: ')
        self.subnode_parse(node)
        self.add_text('\n')

    def do_formula(self, node):
        if len(node.childNodes) == 1:
            f = node.childNodes[0].data
            # check whether this is an inline formula or block formula
            m = self.math_block_re.match(f)
            if m:
                ff = m.group(1)
                # remove any trailing whitespaces
                ff = re.sub(r'^\s+', '', ff)
                ff = re.sub(r'\s+$', '', ff)
                # workaround for suffix replacement in do_para()
                ff = re.sub(r'_([a-zA-Z0-9]|\\\w+)', r'_{\1}', ff)
                self.start_new_paragraph()
                node.childNodes[0].data = ff
                self.add_text('.. math::\n\n')
                self.subnode_parse(node, indent = 2)
                self.add_text('\n')
                return

            m = self.math_inline_re.match(f)
            if m:
                ff = m.group(1)
                # remove any trailing whitespaces
                ff = re.sub(r'^\s+', '', ff)
                ff = re.sub(r'\s+$', '', ff)
                # workaround for suffix replacement in do_para()
                ff = re.sub(r'_([a-zA-Z0-9]|\\\w+)', r'_{\1}', ff)
                node.childNodes[0].data = ff
                self.surround_parse(node, ':math:`', '`')
                return

        self.subnode_parse(node)

# MARK: Para tag handler
    def do_para(self, node):
        """This is the only place where text wrapping is automatically performed.
        Generally, this function parses the node (locally), wraps the text, and
        then adds the result to self.pieces. However, it may be convenient to
        allow the previous content of self.pieces to be included in the text
        wrapping. For this, use the following *convention*:
        If self.pieces ends with '', treat the _previous_ entry as part of the
        current paragraph. Else, insert new-line and start a new paragraph
        and "wrapping context".
        Paragraphs always end with '  \n', but if the parsed content ends with
        the special symbol '', this is passed on.
        """
        if self.pieces[-1:] == ['']:
            pieces, self.pieces = self.pieces[:-2], self.pieces[-2:-1]
        else:
            self.add_text('\n')
            pieces, self.pieces = self.pieces, ['']
        self.subnode_parse(node)
        dont_end_paragraph = self.pieces[-1:] == ['']
        # Now do the text wrapping:
        width = self.textwidth - self.indent
        wrapped_para = []
        for line in ''.join(self.pieces).splitlines():
            keep_markdown_newline = line[-2:] == '  '
            w_line = textwrap.wrap(line, width=width, break_long_words=False)
            if w_line == []:
                w_line = ['']
            if keep_markdown_newline:
                w_line[-1] = w_line[-1] + '  '
            for wl in w_line:
                wrapped_para.append(wl + '\n')
        if wrapped_para:
            if wrapped_para[-1][-3:] != '  \n':
                wrapped_para[-1] = wrapped_para[-1][:-1] + '  \n'
            if dont_end_paragraph:
                wrapped_para.append('')

        for i in range(len(wrapped_para)):
            if len(self.suffixes) > 0:
                # remove suffixes
                for s in self.suffixes:
                    if s in self.constructor_suffixes:
                        repl = r'()'
                    else:
                        repl = r''

                    wrapped_para[i] = re.sub(re.escape(s) + r'\b', repl, wrapped_para[i])
            if len(self.prefixes) > 0:
                prefix = r''
                if self.prefix_name != "":
                    prefix += self.prefix_name + r'.'
                for p in self.prefixes:
                    wrapped_para[i] = wrapped_para[i].replace(p, prefix)

        pieces.extend(wrapped_para)
        self.pieces = pieces

# MARK: List tag handlers
    def do_itemizedlist(self, node):
        if self.listitem == '':
            self.start_new_paragraph()
        elif self.pieces != [] and self.pieces[-1:] != ['']:
            self.add_text('\n')
        listitem = self.listitem
        if self.listitem in ['*', '-']:
            self.listitem = '-'
        else:
            self.listitem = '*'
        self.subnode_parse(node)
        self.listitem = listitem

    def do_orderedlist(self, node):
        if self.listitem == '':
            self.start_new_paragraph()
        elif self.pieces != [] and self.pieces[-1:] != ['']:
            self.add_text('\n')
        listitem = self.listitem
        self.listitem = 0
        self.subnode_parse(node)
        self.listitem = listitem

    def do_listitem(self, node):
        try:
            self.listitem = int(self.listitem) + 1
            item = str(self.listitem) + '. '
        except:
            item = str(self.listitem) + ' '
        self.subnode_parse(node, item, indent=4)

    def do_xrefsect(self, node):
        self.start_new_paragraph()
        title = self.extract_text(self.get_specific_subnodes(node, 'xreftitle'))
        if title == 'Deprecated':
            # try to find out since when its deprecated
            since = self.deprecated_version # default
            m = self.since_version_re.search(self.extract_text(self.get_specific_subnodes(node, 'xrefdescription')))
            if m:
                since = m.group(1)

            self.add_text(['.. deprecated:: ', since, '\n'])
            for nn in self.get_specific_subnodes(node, 'xrefdescription'):
                self.subnode_parse(nn, pieces=[''], indent=4)
        else:
            self.add_text(['**' + title + '**', '\n'])
            for nn in self.get_specific_subnodes(node, 'xrefdescription'):
                self.subnode_parse(nn, pieces=[''], indent=4)

# MARK: Parameter list tag handlers
    def do_parameterlist(self, node):
        name = None
        if node.parentNode and \
           node.parentNode.parentNode and \
           node.parentNode.parentNode.parentNode:
            name = self.extract_text(self.get_specific_subnodes(node.parentNode.parentNode.parentNode, 'name'))

        self.start_new_paragraph()
        text = 'unknown'
        for key, val in node.attributes.items():
            if key == 'kind':
                if val == 'param':
                    text = 'Parameters'
                elif val == 'exception':
                    text = 'Exceptions'
                elif val == 'retval':
                    text = 'Returns'
                else:
                    text = val
                break

#        if name in self.object_binding:
#            print(name, text)
#            return

        if self.indent == 0:
            self.add_text([text, '\n', len(text) * '-', '\n'])
        else:
            self.add_text([text, ':  \n'])
        self.subnode_parse(node)

    def do_parameteritem(self, node):
        member_root = node.parentNode
        for i in range(3):
            member_root = member_root.parentNode
        name = self.extract_text(self.get_specific_subnodes(member_root, 'name'))
        if name in self.object_binding:
            param_name = self.get_specific_subnodes(node, 'parametername', 2)
            try:
                if param_name[0].firstChild.data == self.object_binding[name]['self_arg_name']:
                    return
            except:
                name = 'bla'

        self.subnode_parse(node, pieces=[''])

    def do_parameternamelist(self, node):
        self.subnode_parse(node)

        # try identifying the type of the parameter
        ptype = ""
        if self.with_type_info:
            try:
                pname = self.extract_text(self.get_specific_subnodes(node, 'parametername')[0])
                for p in self.get_specific_subnodes(node.parentNode.parentNode.parentNode.parentNode.parentNode, 'param'):
                    ppname = self.extract_text(self.get_specific_subnodes(p, 'declname'))
                    pptype = self.extract_text(self.get_specific_subnodes(p, 'type'))
                    if pname == ppname:
                        ptype = pptype
                        break
            except:
                ptype = ""

        self.add_text([' : ', ptype, '\n'])
    
    def do_parametername(self, node):
        if self.pieces != [] and self.pieces != ['']:
            self.add_text(', ')
        data = self.extract_text(node)
        self.add_text([data])

    def do_parameterdescription(self, node):
        self.subnode_parse(node, pieces=[''], indent=4)


# MARK: Section tag handler
    def do_simplesect(self, node):
        kind = node.attributes['kind'].value
        if kind in ('date', 'rcs', 'version'):
            return
        self.start_new_paragraph()
        if kind == 'note':
            self.subnode_parse(node, pieces=['Notes', '\n', len('Notes') * '-','\n', ''], indent=0)
            #self.subnode_parse(node, pieces=['**Notes**', '\n',''], indent=4)
        elif kind == 'warning':
            #self.subnode_parse(node, pieces=['**Warnings**', '\n',''], indent=4)
            self.subnode_parse(node, pieces=['Warnings', '\n', len('Warnings') * '-','\n', ''], indent=0)
        elif kind == 'see':
            self.subnode_parse(node, pieces=['See Also', '\n', len('See Also') * '-','\n', ''], indent=0)
        elif kind == 'pre':
            self.subnode_parse(node, pieces=['**Precondition**', '\n',''], indent=4)
        elif kind == 'post':
            self.subnode_parse(node, pieces=['**Postcondition**', '\n',''], indent=4)
        elif kind == 'return':
            pieces=[]
            # determine return type if possible
            p = self.get_specific_parentnode(node, "memberdef", 3)
            if p:
                rettype = self.extract_text(self.get_specific_subnodes(p, 'type'))
            else:
                rettype = "object"

            # remove any occurences pointer marker '*'
            rettype = re.sub(r'\s*\*','', rettype)

            if self.indent == 0:
                self.add_text(['Returns', '\n', len('Returns') * '-', '\n', rettype])
                #pieces = ['Returns', '\n', len('Returns') * '-', '\n', '']
            else:
                self.add_text(['Returns:', '\n', rettype])
                #pieces = ['Returns:', '\n', '']
            self.subnode_parse(node, pieces=pieces, indent = 4)
        else:
            self.subnode_parse(node, pieces=[kind + ': ',''], indent=4)

# MARK: %feature("docstring") producing tag handlers
    def do_compounddef(self, node):
        """This produces %feature("docstring") entries for classes, and handles
        class, namespace and file memberdef entries specially to allow for 
        overloaded functions. For other cases, passes parsing on to standard
        handlers (which may produce unexpected results).
        """
        kind = node.attributes['kind'].value
        if kind in ('class', 'struct'):
            prot = node.attributes['prot'].value
            if prot != 'public':
                return
            self.add_text('\n\n')
            classdefn = self.extract_text(self.get_specific_subnodes(node, 'compoundname'))
            classname = classdefn.split('::')[-1]

            # rename symbol in case we have a type -> typdef mapping
            if classdefn in self.typedef_mapping:
                classdefn = self.typedef_mapping[classdefn]

            self.add_text('%%feature("docstring") %s "\n' % classdefn)

            if self.with_constructor_list:
                constructor_nodes = []
                for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
                    if n.attributes['prot'].value == 'public':
                        if self.extract_text(self.get_specific_subnodes(n, 'definition')) == classdefn + '::' + classname:
                            constructor_nodes.append(n)
                for n in constructor_nodes:
                    self.add_line_with_subsequent_indent(self.get_function_signature(n))

            names = ('briefdescription','detaileddescription')
            sub_dict = self.get_specific_nodes(node, names)
            for n in ('briefdescription','detaileddescription'):
                if n in sub_dict:
                    self.parse(sub_dict[n])
            if self.with_constructor_list:
                self.make_constructor_list(constructor_nodes, classname)
            if self.with_attribute_list:
                self.make_attribute_list(node)

            sub_list = self.get_specific_subnodes(node, 'includes')
            if sub_list:
                self.parse(sub_list[0])
            self.add_text(['";', '\n'])

            names = ['compoundname', 'briefdescription','detaileddescription', 'includes']
            self.subnode_parse(node, ignore = names)

        elif kind in ('file', 'namespace'):
            nodes = node.getElementsByTagName('sectiondef')
            for n in nodes:
                self.parse(n)

        # now explicitely handle possibly overloaded member functions.
        if kind in ['class', 'struct','file', 'namespace', 'group']:
            md_nodes = self.get_memberdef_nodes_and_signatures(node, kind)
            for sig in md_nodes:
                self.handle_typical_memberdefs(sig, md_nodes[sig])
    
    def do_memberdef(self, node):
        """Handle cases outside of class, struct, file or namespace. These are
        now dealt with by `handle_overloaded_memberfunction`.
        Do these even exist???
        """
        prot = node.attributes['prot'].value
        id = node.attributes['id'].value
        kind = node.attributes['kind'].value
        tmp = node.parentNode.parentNode.parentNode
        compdef = tmp.getElementsByTagName('compounddef')[0]
        cdef_kind = compdef.attributes['kind'].value
        if cdef_kind in ('file', 'namespace', 'class', 'struct', 'group'):
            # These cases are now handled by `handle_typical_memberdefs`
            return
        if prot != 'public':
            return
        first = self.get_specific_nodes(node, ('definition', 'name'))
        name = self.extract_text(first['name'])
        if name[:8] == 'operator':  # Don't handle operators yet.
            return
        if not 'definition' in first or kind in ['variable', 'typedef']:
            return

        data = self.extract_text(first['definition'])
        self.add_text('\n')
        self.add_text(['/* where did this entry come from??? */', '\n'])
        self.add_text('%feature("docstring") %s "\n%s' % (data, data))

        for n in node.childNodes:
            if n not in first.values():
                self.parse(n)
        self.add_text(['";', '\n'])
    
# MARK: Entry tag handlers (dont print anything meaningful)
    def do_sectiondef(self, node):
        kind = node.attributes['kind'].value
        if kind in ('public-func', 'func', 'user-defined', 'define'):
            self.subnode_parse(node)

    def do_header(self, node):
        """For a user defined section def a header field is present
        which should not be printed as such, so we comment it in the
        output."""
        data = self.extract_text(node)
        self.add_text('\n/*\n %s \n*/\n' % data)
        # If our immediate sibling is a 'description' node then we
        # should comment that out also and remove it from the parent
        # node's children.
        parent = node.parentNode
        idx = parent.childNodes.index(node)
        if len(parent.childNodes) >= idx + 2:
            nd = parent.childNodes[idx + 2]
            if nd.nodeName == 'description':
                nd = parent.removeChild(nd)
                self.add_text('\n/*')
                self.subnode_parse(nd)
                self.add_text('\n*/\n')

    def do_member(self, node):
        kind = node.attributes['kind'].value
        refid = node.attributes['refid'].value
        if kind == 'function' and refid[:9] == 'namespace':
            self.subnode_parse(node)

    def do_doxygenindex(self, node):
        self.multi = 1
        comps = node.getElementsByTagName('compound')
        for c in comps:
            refid = c.attributes['refid'].value
            fname = refid + '.xml'
            if not os.path.exists(fname):
                fname = os.path.join(self.my_dir,  fname)
            if not self.quiet:
                print("parsing file: %s" % fname)
            p = Doxy2SWIG(fname,
                          with_function_signature = self.with_function_signature,
                          with_type_info = self.with_type_info,
                          with_constructor_list = self.with_constructor_list,
                          with_attribute_list = self.with_attribute_list,
                          with_overloaded_functions = self.with_overloaded_functions,
                          textwidth = self.textwidth,
                          quiet = self.quiet,
                          object_binding = self.object_binding,
                          prefixes = self.prefixes,
                          prefix_name = self.prefix_name,
                          suffixes = self.suffixes,
                          constructor_suffixes = self.constructor_suffixes,
                          typedef_mapping = self.typedef_mapping,
                          ignore_list = self.ignore_list,
                          deprecated_version = self.deprecated_version)
            p.generate()
            self.pieces.extend(p.pieces)

# MARK: main
def main():
    usage = __doc__
    parser = optparse.OptionParser(usage)
    parser.add_option("-f", '--function-signature',
                      action='store_true',
                      default=False,
                      dest='f',
                      help='include function signature in the documentation. This is handy when not using swig auto-generated function definitions %feature("autodoc", [0,1])')
    parser.add_option("-t", '--type-info',
                      action='store_true',
                      default=False,
                      dest='t',
                      help='include type information for arguments in function signatures. This is similar to swig autodoc level 1')
    parser.add_option("-c", '--constructor-list',
                      action='store_true',
                      default=False,
                      dest='c',
                      help='generate a constructor list for class documentation. Useful for target languages where the object construction should be documented in the class documentation.')
    parser.add_option("-a", '--attribute-list',
                      action='store_true',
                      default=False,
                      dest='a',
                      help='generate an attributes list for class documentation. Useful for target languages where class attributes should be documented in the class documentation.')
    parser.add_option("-o", '--overloaded-functions',
                      action='store_true',
                      default=False,
                      dest='o',
                      help='collect all documentation for overloaded functions. Useful for target languages that have no concept of overloaded functions, but also to avoid having to attach the correct docstring to each function overload manually')
    parser.add_option("-w", '--width', type="int",
                      action='store',
                      dest='w',
                      default=80,
                      help='textwidth for wrapping (default: 80). Note that the generated lines may include 2 additional spaces (for markdown).')
    parser.add_option("-q", '--quiet',
                      action='store_true',
                      default=False,
                      dest='q',
                      help='be quiet and minimize output')
    parser.add_option("-b", '--object-binding',
                      type = "string",
                      action = 'store',
                      dest = 'b',
                      default = None,
                      help = 'load file that specifies which functions are bound as method to which objects')

    parser.add_option("-p", '--prefix',
                      type = "string",
                      action = 'append',
                      dest = 'p',
                      default = None,
                      help = 'Prefixes to remove or substitute')

    parser.add_option('--prefix-name',
                      type = "string",
                      action = 'store',
                      dest = 'prefix_name',
                      default = None,
                      help = 'Prefix name')

    parser.add_option("-s", '--suffix',
                      type = "string",
                      action = 'append',
                      dest = 's',
                      default = None,
                      help = 'Suffixes to remove')

    parser.add_option('--constructor-suffix',
                      type = "string",
                      action = 'append',
                      dest = 'cs',
                      default = None,
                      help = 'Constructor suffixes')

    parser.add_option('-m', '--typedef-mapping',
                      type = "string",
                      action = 'store',
                      dest = 'm',
                      default = None,
                      help = 'load file that specifies which types should be replaced by a typedef')

    parser.add_option('-d', '--deprecated-version',
                      type = "string",
                      action = 'store',
                      dest = 'd',
                      default = '2.6.0',
                      help = 'default version for deprecation warnings')

    parser.add_option("-i", '--ignored',
                      type = "string",
                      action = 'store',
                      dest = 'i',
                      default = None,
                      help = 'load file that specifies which functions are ignored and possibly overwritten')

    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("no input and output specified")

    object_binding  = dict()
    if options.b:
        with io.open(options.b, "r", newline = '') as f:
            data = csv.reader(f, delimiter = ",")
            header = next(data, None)
            for row in data:
                object_binding[row[0]] = {'method': row[1],
                                          'class': row[2],
                                          'self_arg_name': row[3]}

    typedef_mapping  = dict()
    if options.m:
        with io.open(options.m, "r", newline = '') as f:
            data = csv.reader(f, delimiter = ",")
            header = next(data, None)
            for row in data:
                typedef_mapping[row[0]] = row[1]

    ignore_list = dict()
    if options.i:
        with io.open(options.i, "r", newline = '') as f:
            data = csv.reader(f, delimiter = ",")
            header = next(data, None)
            for row in data:
                ignore_list[row[0]] = {'replacement': row[1] if row[1] != "None" else None }

    p = Doxy2SWIG(args[0],
                  with_function_signature = options.f,
                  with_type_info = options.t,
                  with_constructor_list = options.c,
                  with_attribute_list = options.a,
                  with_overloaded_functions = options.o,
                  textwidth = options.w,
                  quiet = options.q,
                  object_binding = object_binding,
                  prefixes = options.p,
                  prefix_name = options.prefix_name,
                  suffixes = options.s,
                  constructor_suffixes = options.cs,
                  typedef_mapping = typedef_mapping,
                  ignore_list = ignore_list,
                  deprecated_version = options.d)
    p.generate()
    p.write(args[1])

if __name__ == '__main__':
    main()
