#!/usr/bin/env python3
#

import re
import argparse

indent = 4

def process_lines(filehandle):
    comment_line_pat = re.compile(r'^.\\\"')
    comment_inline_pat = re.compile(r'^(.*?)\\\"(.*)$')
    section_pat = re.compile(r'^.SH\s*(.*)$')
    subsection_pat = re.compile(r'^.SS\s*(.*)$')
    bold_pat = [
        re.compile(r'\\fB(.*?)\\fR'),
        re.compile(r'\\fB(.*?)\\fP'),
        re.compile(r'\\f3(.*?)\\fR'),
        re.compile(r'\\f3(.*?)\\fR'),
    ]
    ital_pat = [
        re.compile(r'\\fI\\,(.*?)\\/\\fR'),
        re.compile(r'\\fI(.*?)\\fR'),
        re.compile(r'\\fI\\,(.*?)\\/\\fP'),
        re.compile(r'\\fI(.*?)\\fP'),
        re.compile(r'\\f2(.*?)\\fR'),
        re.compile(r'\\f2(.*?)\\fR'),
    ]

    cw_pat = [
        re.compile(r'\\f\(CW(.*?)\\fR'),
    ]

    ital_line_pat = re.compile(r'.I\s+(.*)$')
    dblquoted_text_pat = re.compile(r'`(.*?)\'')
    quoted_text_pat = re.compile(r'\'(.*?)\'')
    shell_esc_pat = re.compile(r'\$\'\\e(.*?)\'')

    current_indent = 0
    indent_next_line = 0
    content = []
    # go through man page line by line
    for line in filehandle:
        line = line.strip()

        # reset indent whenever we encounter a command
        if line.startswith("."):
            current_indent = 0

        # skip comment lines
        if comment_line_pat.match(line):
            continue

        # skip table header
        if line.startswith(".TH"):
            continue

        if line.startswith(".PP"):
            content.append("")
            continue

        # remove trailing comments from line
        m = comment_inline_pat.match(line)
        if m:
            line = str(m.group(1))

        # reset indent whenever we have a blank line
        if line.strip() == "":
            current_indent = 0

        m = section_pat.match(line)
        if m:
            section_head = m.group(1).strip("\"").strip()
            content.append("\n{}\n{}\n".format(section_head, "-" * len(section_head)))
            continue

        m = subsection_pat.match(line)
        if m:
            subsection_head = m.group(1).strip("\"").strip()
            content.append("\n{}\n{}\n".format(subsection_head, "^" * len(subsection_head)))
            continue

        # do a lot of replacements
        line = line.replace('\\-', '-').replace('\n\r', '\n')

        for bp in bold_pat:
            line = bp.sub(r'**\1**', line)

        for ip in ital_pat:
            line = ip.sub(r'*\1*', line)

        # usually, inline constant width is some error from help2man
        # so we simply remove it
        for cwp in cw_pat:
            line = cwp.sub(r'\1', line)

        if indent_next_line == 1:
            current_indent = indent

        if indent_next_line > 0:
            indent_next_line -= 1

        # indent paragraphs if requested
        if line.strip().startswith(".IP"):
            current_indent = indent
            content.append("\n")
            continue

        if line.startswith(".TP"):
            indent_next_line = 2
            content.append(line)
            continue

        m = ital_line_pat.match(line)
        if m:
            line = "*" + m.group(1) + "*"

        line = shell_esc_pat.sub(r'##-##\1##-##', line)
        line = dblquoted_text_pat.sub(r'"\1"', line)
        line = quoted_text_pat.sub(r'``\1``', line)

        line = re.sub(r'\\&', r'', line)
        line = re.sub(r'\\e', r'\\', line)
        line = re.sub(r'(\(default=(.*?)\))', r'*\1*', line)

        content.append(' ' * current_indent + line)

    return "\n".join(content)

def stage_2(content):
    # Process tagged paragraphs
    #
    # Here, we assume tagged paragraphs (.TP) only appear for
    # options.
    tagged_par = re.compile(r'^\.TP\s*\n\s*(.+?)\n(.*?)(?=\n\s*\n|\n\.)',
                            re.MULTILINE|re.DOTALL)

    def prepare_option(m):
        return "\n.. option:: {}\n\n{}{}".format(m.group(1),
                                                 ' ' * indent,
                                                 " ".join([t for t in m.group(2).split() if t]))

    content = tagged_par.sub(lambda m: prepare_option(m), content)

    # Repair bold-face and italic parts in .. options::
    #
    pat = re.compile(r'^\.\. option::\s+(.*)$', re.MULTILINE)
    def repair_opt(m):
        opt = m.group(1)
        bold_pat = re.compile(r'\*\*(.*?)\*\*')
        ital_pat = re.compile(r'\*(.*?)\*')
        opt = bold_pat.sub(r'\1', opt)
        opt = ital_pat.sub(r'\1', opt)
        return ".. option:: " + opt

    content = pat.sub(lambda m: repair_opt(m), content)


    # Replace constant-width/code blocks
    pat = re.compile(r'\.nf\n\.ft\s+CW\n(.*?)\.ft\n\.fi', re.MULTILINE | re.DOTALL)

    def code_blocks(m):
        block = ".. code::\n\n"
        for line in m.group(1).split('\n'):
            block += ' ' * 4 + line + '\n'
        return block

    content = pat.sub(lambda m: code_blocks(m), content)

    # xref inline options
    # those should always be bold-faced in our input
    pat = [
        re.compile(r'\*\*(\-\-(.*?))\*\*', re.MULTILINE | re.DOTALL), # --opt
        re.compile(r'\*\*(\-((\S)(.*?)))\*\*', re.MULTILINE | re.DOTALL) # -o
    ]

    def xref_opt(m):
        opt = ":option:"
#        try:
#            if m.group(4):
#                opt += f"`-{m.group(3)}`{m.group(4)}"
#            else:
#                opt += f"`{m.group(1)}`"
#        except:
#            opt += f"`{m.group(1)}`"
        opt += f"`{m.group(1)}`"
        return opt

    for p in pat:
        content = p.sub(lambda m: xref_opt(m), content)

    # make page title and brief description from
    #     NAME
    #     ----
    # entry
    pat = re.compile(r'\s*NAME\n----\n\s*(.*?)\s+(.*?)\n\s*\n', re.MULTILINE | re.DOTALL)

    def page_title(m):
        prog = m.group(1)
        desc = m.group(2)
        return "{}\n{}\n{}\n\n:program:`{}` {}\n\n".format('#' * len(prog), prog, '#' * len(prog), prog, desc)

    content = pat.sub(lambda m: page_title(m), content)

    # undo shell escape safety replacement
    pat = re.compile(r'##-##(.*?)##-##')
    content = pat.sub(r"``$'\\\1'``", content)

    # Finally, lets re-format the synopsis
    pat = re.compile(r'synopsis\n--------\n+\.B\s+(.*?)\n(.*?)\n\n', re.I | re.MULTILINE | re.DOTALL)

    def synopsis(m):
        prog = m.group(1)
        opt  = m.group(2)
        opt = re.sub(r'\*(.*?)\*', r'\1', opt)
        return "Synopsis\n--------\n\n.. code:: bash\n\n    {} {}\n\n".format(prog, opt)

    content = pat.sub(lambda m: synopsis(m), content)

    return content


parser = argparse.ArgumentParser(description="Convert manpages to ReStructuredText")
parser.add_argument("-i", "--input", type=str, help="The input man page file", required=True)
parser.add_argument("-o", "--output", type=str, help="The output .rst file", required=True)

args = parser.parse_args()


with open(args.input, 'r', encoding="utf8", errors='ignore') as f:
    content = process_lines(f)
    content = stage_2(content)

    with open(args.output, "w") as of:
        of.write(content)

