#...............................................................................
#
#  This file is part of the Doxyrest toolkit.
#
#  Doxyrest is distributed under the MIT license.
#  For details see accompanying license.txt file,
#  the public copy of which is also available at:
#  http://tibbo.com/downloads/archive/doxyrest/license.txt
#
#...............................................................................

import os
import re
import warnings
from packaging import version
from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.transforms import Transform
from sphinx import __version__ as sphinx_version_string, roles, addnodes, config
from sphinx.directives.other import Include
from sphinx.domains import Domain

#...............................................................................
#
#  utils
#

sphinx_version = version.parse(sphinx_version_string)
this_dir = os.path.dirname(os.path.realpath(__file__))
url_re_prog = re.compile('(ftp|https?)://')
crefdb = {}
cref_w_target_re_prog = re.compile('(.+?)\s*<([^<>]*)>$')

def get_cref_target(text, target=None):
    if not target:
        target = text
    elif url_re_prog.match(target):
        return target

    if target in crefdb:
        return crefdb[target]

    warnings.warn('target not found for cref: ' + target, Warning, 2)
    return None

def get_cref_target_ex(text):
    match = cref_w_target_re_prog.match(text)
    if match:
        text = match.group(1)
        target = match.group(2)
    else:
        target = text

    return get_cref_target(text, target), text

if sphinx_version >= version.parse('1.8.0'):
    def add_js_file(app, filename):
        app.add_js_file(filename)

    def add_css_file(app, filename):
        app.add_css_file(filename)
else:
    def add_js_file(app, filename):
        app.add_javascript(filename)

    def add_css_file(app, filename):
        app.add_stylesheet(filename)

#...............................................................................
#
#  Sphinx nodes
#

class HighlightedText(nodes.Inline, nodes.TextElement):
    """Custom inline node for syntax-highlighted text fragments."""
    pass

class DoxyrestBlock(nodes.General, nodes.Element):
    """Container for the entire code block (<pre> wrapper)."""
    pass

class DoxyrestLine(nodes.Inline, nodes.TextElement):
    """Container for a single line (<span data-line> wrapper)."""
    pass

def visit_doxyrest_block(self, node):
    self.body.append(self.starttag(node, 'div', CLASS='highlight'))
    self.body.append('<pre>')

def depart_doxyrest_block(self, node):
    self.body.append('</pre></div>\n')

def visit_doxyrest_line(self, node):
    lineno = node.get('lineno', '')
    self.body.append(f'<span data-line="{lineno}">')

def depart_doxyrest_line(self, node):
    # CRITICAL FIX: Append \n so the <pre> block renders a line break
    self.body.append('</span>\n')

def cleanup_fragment(html):
    """
    Strips block-level wrappers (div, pre, data-line spans) from syntax-highlighted
    fragments to ensure they render as inline elements.
    """
    if not html:
        return html

    html = re.sub(r'^<div[^>]*>\s*<pre[^>]*>', '', html)
    html = re.sub(r'</pre>\s*</div>$', '', html)
    html = html.replace('<span></span>', '')

    match = re.match(r'^\s*(<span\s+[^>]*data-line=["\']\d+["\'][^>]*>)(.*)(</span>)\s*$', html, re.DOTALL)
    if match:
        html = match.group(2)

    if html.endswith('\n'):
        html = html[:-1]

    return html

def visit_highlighted_text_node(self, node):
    text_node = node.children[0]
    language = node['language']

    if language == 'none':
        self.body.append(text_node)
    else:
        options = {'nowrap': True}
        text = text_node
        highlighted = self.highlighter.highlight_block(text, language, options)
        clean_html = cleanup_fragment(highlighted)
        self.body.append(clean_html)

    raise nodes.SkipNode

def create_ref_node(raw_text, text, target, highlight_language=None):
    if not target:
        if highlight_language:
            return HighlightedText(text, text, language=highlight_language)
        return nodes.Text(text, text)

    style_opts = {'style': 'text-decoration: underline'}

    if url_re_prog.match(target):
        node = nodes.reference(raw_text, '', refuri=target, **style_opts)
    else:
        node = addnodes.pending_xref(raw_text)
        node['reftype'] = 'ref'
        node['refdomain'] = 'std'
        node['reftarget'] = target
        node['refwarn'] = True
        node['refexplicit'] = True
        node.attributes.update(style_opts)

    node['classes'] += ['doxyrest-code-link']

    if highlight_language:
        node += HighlightedText(text, text, language=highlight_language)
    else:
        node += nodes.Text(text, text)

    return node

def create_target_node(raw_text, text, target, highlight_language, lineno, document, extra_classes=[]):
    node = nodes.target(raw_text, '', ids=[target], names=[target])
    node['classes'] += extra_classes
    node.line = lineno
    document.note_explicit_target(node)

    result = [node]

    if text:
        if highlight_language:
            result.append(HighlightedText(text, text, language=highlight_language))
        else:
            result.append(nodes.Text(text, text))

    return result


#...............................................................................
#
#  Sphinx directives
#

class RefCodeBlock(Directive):
    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = False

    option_spec = {
        'linenos': directives.flag,
        'dedent': int,
        'lineno-start': int,
        'emphasize-lines': directives.unchanged_required,
        'caption': directives.unchanged_required,
        'class': directives.class_option,
        'name': directives.unchanged,
    }

    def __init__(self, *args, **kwargs):
        Directive.__init__(self, *args, **kwargs)

        role_re_src = '(:ref:|:cref:|:target:)'
        if self.state.document.settings.env.config.default_role == 'cref':
            role_re_src += '?' # explicit role is optional

        role_re_src += '`(.+?)(\s*<([^<>]*)>)?`'
        self.role_re_prog = re.compile(role_re_src, re.DOTALL)

    def run(self):
        config = self.state.document.settings.env.config
        code_lines = self.content
        base_lineno = self.options.get('lineno-start', 1)

        if len(self.arguments) >= 1:
            language = self.arguments[0]
        else:
            language = config.highlight_language

        block_node = DoxyrestBlock()
        block_node['classes'] += self.options.get('class', [])

        for i, line in enumerate(code_lines):
            lineno = base_lineno + i
            line_node = DoxyrestLine(lineno=lineno)

            # Preserve empty line height
            if not line:
                line_node += nodes.Text(" ")
                block_node += line_node
                continue

            pos = 0
            while True:
                match = self.role_re_prog.search(line, pos)
                if match is None:
                    chunk = line[pos:]
                    if chunk:
                        line_node += HighlightedText(chunk, chunk, language=language)
                    break

                pre_text = line[pos:match.start()]
                if pre_text:
                    line_node += HighlightedText(pre_text, pre_text, language=language)

                raw_text = match.group(0)
                role = match.group(1)
                text = match.group(2)
                target = match.group(4)

                if text:
                    text = text.replace('\\<', '<')

                if role == ':target:':
                    if not target: target = text
                    line_node += create_target_node(raw_text, text, target, language, None, self.state.document, ['doxyrest-code-target'])
                else:
                    if not role or role == ':cref:':
                        target = get_cref_target(text, target)
                    elif not target:
                        target = text
                    line_node += create_ref_node(raw_text, text, target, language)

                pos = match.end()

            block_node += line_node

        self.add_name(block_node)
        return [block_node]


#...............................................................................
#
#  Sphinx transforms
#

class RefTransform(Transform):
    default_priority = 100

    node_classes = {
        nodes.literal,
        nodes.strong,
        nodes.emphasis
    }

    def __init__(self, document, startnode=None):
        Transform.__init__(self, document, startnode)

        re_src = '(:c?ref:)'
        if document.settings.env.config.default_role == 'cref':
            re_src += '?' # explicit role is optional

        re_src += '`(.+?)(\s*<([^<>]*)>)?`'
        self.re_prog = re.compile(re_src)

    @staticmethod
    def node_filter(node):
        for node_class in RefTransform.node_classes:
            if isinstance (node, node_class):
                return node['classes'] == []

        return False

    def apply(self):
        for node in self.document.traverse(RefTransform.node_filter):
            code = node.astext()
            node.children = []
            pos = 0

            while True:
                match = self.re_prog.search(code, pos)
                if match is None:
                    plain_text = code[pos:]
                    if plain_text != "":
                        node += nodes.Text(plain_text, plain_text)
                    break

                plain_text = code[pos:match.start()]
                if plain_text != "":
                    node += nodes.Text(plain_text, plain_text)

                raw_text = match.group(0)
                role = match.group(1)
                text = match.group(2)
                target = match.group(4)

                if not role or role == ':cref:':
                    target = get_cref_target(text, target)

                node += create_ref_node(raw_text, text, target)[0] # Take first node
                pos = match.end()


#...............................................................................
#
#  Sphinx roles
#

def cref_role(typ, raw_text, text, lineno, inliner, options={}, content=[]):
    target, text = get_cref_target_ex(text)
    node = nodes.literal(raw_text, '')
    node['classes'] += ['doxyrest-cref']
    ref_node = create_ref_node(raw_text, text, target)
    if isinstance(ref_node, list):
         node += ref_node[0]
    else:
         node += ref_node
    return [node], []

def target_role(typ, raw_text, text, lineno, inliner, options={}, content=[]):
    nodes_list = create_target_node(raw_text, None, text, None, lineno, inliner.document)
    return nodes_list, []


#...............................................................................
#
#  Doxyrest domain
#

class DoxyrestDomain(Domain):
    name = 'doxyrest'
    label = 'Doxyrest'

    def merge_domaindata(self, docnames, otherdata):
        pass

    def resolve_any_xref(self, env, fromdocname, builder, target, node, contnode):
        cref_target = get_cref_target(target, target)
        if not cref_target:
            return []

        std = env.get_domain('std')
        node['refexplicit'] = True
        resolved_node = std.resolve_xref(env, fromdocname, builder, 'ref', cref_target, node, contnode)
        if not resolved_node:
            return []

        result_node = nodes.literal(target, '')
        result_node += resolved_node
        return [('std:ref', result_node)]


#...............................................................................
#
#  Sphinx source inputs
#

is_sphinx_tab_aware = sphinx_version >= version.parse('2.1.0')

if not is_sphinx_tab_aware:
    from sphinx.io import SphinxBaseFileInput, SphinxRSTFileInput
    from docutils.statemachine import StringList, string2lines


    class TabAwareSphinxRSTFileInput(SphinxRSTFileInput):
        def read(self):
            # type: () -> StringList
            inputstring = SphinxBaseFileInput.read(self)
            tab_width = self.env.config.doxyrest_tab_width
            lines = string2lines(inputstring, convert_whitespace=True, tab_width=tab_width)

            content = StringList()
            for lineno, line in enumerate(lines):
                content.append(line, self.source_path, lineno)

            if self.env.config.rst_prolog:
                self.prepend_prolog(content, self.env.config.rst_prolog)
            if self.env.config.rst_epilog:
                self.append_epilog(content, self.env.config.rst_epilog)

            return content


    class TabAwareInclude(Include):
        def run(self):
            # update tab_width setting
            self.state.document.settings.tab_width = \
            self.state.document.settings.env.config.doxyrest_tab_width
            return Include.run(self)

#...............................................................................
#
#  Sphinx events
#

def on_builder_inited(app):
    app.config.html_static_path += [
        this_dir + '/css/doxyrest-pygments.css',
        this_dir + '/js/target-highlight.js'
    ]

    add_css_file(app, 'doxyrest-pygments.css')
    add_js_file(app, 'target-highlight.js')

    supported_themes = {
        'sphinx_rtd_theme',
        'sphinxdoc',
        'shibuya'
    }

    if app.config.html_theme in supported_themes:
        css_file = 'doxyrest-' + app.config.html_theme + '.css'
        app.config.html_static_path += [this_dir + '/css/' + css_file];
        add_css_file(app, css_file);

    for basedir, dirnames, filenames in os.walk(app.srcdir):
        if 'crefdb.py' in filenames:
            crefdb_path = os.path.join(basedir, 'crefdb.py')
            src = open(crefdb_path).read()
            ns = {}
            exec(src, ns)
            new_crefdb = ns['crefdb']
            if isinstance(new_crefdb, dict):
                global crefdb
                crefdb.update(new_crefdb)

def on_config_inited(app, config):
    docutils_conf_in_path = this_dir + '/conf/doxyrest-docutils.conf.in'
    docutils_conf_path = app.doctreedir + '/doxyrest-docutils.conf'

    src_file = open(docutils_conf_in_path, 'r')
    contents = src_file.read()
    contents = contents.replace('%tab_width%', str(config.doxyrest_tab_width))
    src_file.close()

    if not os.path.exists(app.doctreedir):
        os.makedirs(app.doctreedir)

    dst_file = open(docutils_conf_path, 'w')
    dst_file.write(contents)
    dst_file.close()

    if 'DOCUTILSCONFIG' in os.environ:
        prev_docutils_conf = os.environ['DOCUTILSCONFIG']
        os.environ['DOCUTILSCONFIG'] = docutils_conf_path + os.pathsep + prev_docutils_conf
    else:
        os.environ['DOCUTILSCONFIG'] = docutils_conf_path

#...............................................................................
#
#  Doxyrest extenstion setup
#

def setup(app):
    app.add_domain(DoxyrestDomain)

    app.add_node(DoxyrestBlock,
        html=(visit_doxyrest_block, depart_doxyrest_block))
    app.add_node(DoxyrestLine,
        html=(visit_doxyrest_line, depart_doxyrest_line))
    app.add_node(HighlightedText,
        html=(visit_highlighted_text_node, None),
        latex=(visit_highlighted_text_node, None))

    app.add_role('cref', cref_role)
    app.add_role('target', target_role)
    app.add_config_value('doxyrest_cref_file', default=None, rebuild=True)
    app.add_config_value('doxyrest_tab_width', default=4, rebuild=True)
    directives.register_directive('ref-code-block', RefCodeBlock)
    app.add_transform(RefTransform)
    app.connect('builder-inited', on_builder_inited)
    app.connect('config-inited', on_config_inited)

    if not is_sphinx_tab_aware:
        app.registry.source_inputs['restructuredtext'] = TabAwareSphinxRSTFileInput
        directives.register_directive('include', TabAwareInclude)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
