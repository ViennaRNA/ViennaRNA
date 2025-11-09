"""
Python module to parse Doxygen's XML output.

This module defines the following types:

- TaggedValue

    TaggedValue is a name/value pair for tagged union values.

- Node

    Node is an empty class used as a base type for all classes that start with
    "Node_".

- Node_X

    These classes are generated according the input schema.

- ListItem_X

    Types that have "kind" equal to "tuple_list_element" in the schema also have
    a companion class for their elements. It will have the same name as the main
    class except it starts with ListItem_ instead of Node_. These are named
    tuples.

- ParseError

    The exception raised when there is a problem with the XML input that cannot
    be ignored.

- ParseWarning

    The warning class for possible problems in the XML input. Currently this
    is only issued for unexpected elements and attributes.


Each non-simple type has some or all of the following entities:

- _node_class_attr__X

    Attribute handlers for element X.

    This is a mapping of attribute names to functions that handle the
    attributes.

- def _node_class_attr_end__X(state: _ParseState, obj)

    This is called after all attributes are handled. It is used to check for
    unset fields and to fill them with default values or raise an exception.
    This is a separate function so that derived elements can use it.

- _node_class_child__X

    Element handlers for element X.

    This is a mapping of element names to functions that handle the elements.

- def _node_class_finish_fields__X(state: _ParseState, obj)

    This is called after all child elements are handled. It is used to check for
    unset fields and to fill them with default values or raise an exception.
    This is a separate function so that derived elements can use it.


If the type is used directly, it will have its own class and the following
function:

- def _node_class_start__X(state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]]):

    This has three responsibilities:
    - To create the object corresponding to element X.
    - To handle the XML attributes.
    - To add the new object and the appropriate XML event handlers to the top of
      the parser stack.

    This function doesn't return a value immediately, instead "setter" is called
    with the value when it's ready.

"""

from __future__ import annotations

import enum
import functools
import warnings
from typing import (
    TYPE_CHECKING,
    Literal,
    NamedTuple,
)
from xml.parsers import expat

if TYPE_CHECKING:
    import sys
    from collections.abc import Iterable, Sequence
    from typing import (
        Any,
        Callable,
        ClassVar,
        Generic,
        NoReturn,
        TypeVar,
        Union,
    )

    if sys.version_info >= (3, 11):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

    T = TypeVar("T")
    T_covar = TypeVar("T_covar", covariant=True)
    U_covar = TypeVar("U_covar", covariant=True)


class ParseError(RuntimeError):
    @property
    def message(self, /) -> str:
        return self.args[0]

    @property
    def lineno(self, /) -> int:
        return self.args[1]

    def __str__(self, /) -> str:
        if self.lineno is None:
            return "Error: " + self.message
        return f"Error on line {self.lineno}: {self.message}"


class ParseWarning(UserWarning):
    pass


class Node:
    __slots__ = ()

    _fields: ClassVar[tuple[str, ...]]


# This needs to run on Python 3.8, where built-in types don't implement
# __class_getitem__, and Python 3.9 and 3.10, which don't allow
# multiple-inheritance with NamedTuple.
if TYPE_CHECKING:

    class ListNode(list[T], Node, Generic[T]): ...

    class TaggedValue(NamedTuple, Generic[T_covar, U_covar]):
        name: T_covar
        value: U_covar

else:

    class TaggedValue(NamedTuple):
        name: str
        value: Any

        __class_getitem__ = classmethod(lambda cls, x: cls)

    class ListNode(list, Node):
        __slots__ = ()

        __class_getitem__ = classmethod(lambda cls, x: cls)


class Node_DoxygenTypeIndex(Node):
    __slots__ = (
        "version",
        "compound",
    )

    _fields = __slots__

    def __init__(
        self,
        version: str,
        compound: Iterable[Node_CompoundType] = (),
    ):  # pragma: no cover
        self.version = version
        self.compound = (
            compound if _GLOBAL_type(compound) is _GLOBAL_list else _GLOBAL_list(compound)
        )


class Node_CompoundType(Node):
    __slots__ = (
        "refid",
        "kind",
        "name",
        "member",
    )

    _fields = __slots__

    def __init__(
        self,
        refid: str,
        kind: CompoundKind,
        name: str,
        member: Iterable[Node_MemberType] = (),
    ):  # pragma: no cover
        self.refid = refid
        self.kind = kind
        self.name = name
        self.member = member if _GLOBAL_type(member) is _GLOBAL_list else _GLOBAL_list(member)


class Node_DoxygenType(Node):
    __slots__ = (
        "version",
        "compounddef",
    )

    _fields = __slots__

    def __init__(
        self,
        version: str,
        compounddef: Iterable[Node_compounddefType] = (),
    ):  # pragma: no cover
        self.version = version
        self.compounddef = (
            compounddef if _GLOBAL_type(compounddef) is _GLOBAL_list else _GLOBAL_list(compounddef)
        )


class Node_compounddefType(Node):
    __slots__ = (
        "abstract",
        "final",
        "id",
        "inline",
        "kind",
        "language",
        "prot",
        "sealed",
        "basecompoundref",
        "briefdescription",
        "collaborationgraph",
        "compoundname",
        "derivedcompoundref",
        "detaileddescription",
        "exports",
        "incdepgraph",
        "includedby",
        "includes",
        "inheritancegraph",
        "initializer",
        "innerclass",
        "innerconcept",
        "innerdir",
        "innerfile",
        "innergroup",
        "innermodule",
        "innernamespace",
        "innerpage",
        "invincdepgraph",
        "listofallmembers",
        "location",
        "programlisting",
        "qualifier",
        "requiresclause",
        "sectiondef",
        "tableofcontents",
        "templateparamlist",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        id: str,
        kind: DoxCompoundKind,
        compoundname: str,
        abstract: bool | None = None,
        final: bool | None = None,
        inline: bool | None = None,
        language: DoxLanguage | None = None,
        prot: DoxProtectionKind | None = None,
        sealed: bool | None = None,
        basecompoundref: Iterable[Node_compoundRefType] = (),
        briefdescription: Node_descriptionType | None = None,
        collaborationgraph: Node_graphType | None = None,
        derivedcompoundref: Iterable[Node_compoundRefType] = (),
        detaileddescription: Node_descriptionType | None = None,
        exports: Node_exportsType | None = None,
        incdepgraph: Node_graphType | None = None,
        includedby: Iterable[Node_incType] = (),
        includes: Iterable[Node_incType] = (),
        inheritancegraph: Node_graphType | None = None,
        initializer: Node_linkedTextType | None = None,
        innerclass: Iterable[Node_refType] = (),
        innerconcept: Iterable[Node_refType] = (),
        innerdir: Iterable[Node_refType] = (),
        innerfile: Iterable[Node_refType] = (),
        innergroup: Iterable[Node_refType] = (),
        innermodule: Iterable[Node_refType] = (),
        innernamespace: Iterable[Node_refType] = (),
        innerpage: Iterable[Node_refType] = (),
        invincdepgraph: Node_graphType | None = None,
        listofallmembers: Node_listofallmembersType | None = None,
        location: Node_locationType | None = None,
        programlisting: Node_listingType | None = None,
        qualifier: Iterable[str] = (),
        requiresclause: Node_linkedTextType | None = None,
        sectiondef: Iterable[Node_sectiondefType] = (),
        tableofcontents: Node_tableofcontentsType | None = None,
        templateparamlist: Node_templateparamlistType | None = None,
        title: str | None = None,
    ):  # pragma: no cover
        self.abstract = abstract
        self.final = final
        self.id = id
        self.inline = inline
        self.kind = kind
        self.language = language
        self.prot = prot
        self.sealed = sealed
        self.basecompoundref = (
            basecompoundref
            if _GLOBAL_type(basecompoundref) is _GLOBAL_list
            else _GLOBAL_list(basecompoundref)
        )
        self.briefdescription = briefdescription
        self.collaborationgraph = collaborationgraph
        self.compoundname = compoundname
        self.derivedcompoundref = (
            derivedcompoundref
            if _GLOBAL_type(derivedcompoundref) is _GLOBAL_list
            else _GLOBAL_list(derivedcompoundref)
        )
        self.detaileddescription = detaileddescription
        self.exports = exports
        self.incdepgraph = incdepgraph
        self.includedby = (
            includedby if _GLOBAL_type(includedby) is _GLOBAL_list else _GLOBAL_list(includedby)
        )
        self.includes = (
            includes if _GLOBAL_type(includes) is _GLOBAL_list else _GLOBAL_list(includes)
        )
        self.inheritancegraph = inheritancegraph
        self.initializer = initializer
        self.innerclass = (
            innerclass if _GLOBAL_type(innerclass) is _GLOBAL_list else _GLOBAL_list(innerclass)
        )
        self.innerconcept = (
            innerconcept
            if _GLOBAL_type(innerconcept) is _GLOBAL_list
            else _GLOBAL_list(innerconcept)
        )
        self.innerdir = (
            innerdir if _GLOBAL_type(innerdir) is _GLOBAL_list else _GLOBAL_list(innerdir)
        )
        self.innerfile = (
            innerfile if _GLOBAL_type(innerfile) is _GLOBAL_list else _GLOBAL_list(innerfile)
        )
        self.innergroup = (
            innergroup if _GLOBAL_type(innergroup) is _GLOBAL_list else _GLOBAL_list(innergroup)
        )
        self.innermodule = (
            innermodule if _GLOBAL_type(innermodule) is _GLOBAL_list else _GLOBAL_list(innermodule)
        )
        self.innernamespace = (
            innernamespace
            if _GLOBAL_type(innernamespace) is _GLOBAL_list
            else _GLOBAL_list(innernamespace)
        )
        self.innerpage = (
            innerpage if _GLOBAL_type(innerpage) is _GLOBAL_list else _GLOBAL_list(innerpage)
        )
        self.invincdepgraph = invincdepgraph
        self.listofallmembers = listofallmembers
        self.location = location
        self.programlisting = programlisting
        self.qualifier = (
            qualifier if _GLOBAL_type(qualifier) is _GLOBAL_list else _GLOBAL_list(qualifier)
        )
        self.requiresclause = requiresclause
        self.sectiondef = (
            sectiondef if _GLOBAL_type(sectiondef) is _GLOBAL_list else _GLOBAL_list(sectiondef)
        )
        self.tableofcontents = tableofcontents
        self.templateparamlist = templateparamlist
        self.title = title


class Node_graphType(Node):
    __slots__ = ("node",)

    _fields = __slots__

    def __init__(
        self,
        node: Iterable[Node_nodeType] = (),
    ):  # pragma: no cover
        self.node = node if _GLOBAL_type(node) is _GLOBAL_list else _GLOBAL_list(node)


class Node_templateparamlistType(Node):
    __slots__ = ("param",)

    _fields = __slots__

    def __init__(
        self,
        param: Iterable[Node_paramType] = (),
    ):  # pragma: no cover
        self.param = param if _GLOBAL_type(param) is _GLOBAL_list else _GLOBAL_list(param)


class Node_sectiondefType(Node):
    __slots__ = (
        "kind",
        "description",
        "header",
        "member",
        "memberdef",
    )

    _fields = __slots__

    def __init__(
        self,
        kind: DoxSectionKind,
        description: Node_descriptionType | None = None,
        header: str | None = None,
        member: Iterable[Node_MemberType] = (),
        memberdef: Iterable[Node_memberdefType] = (),
    ):  # pragma: no cover
        self.kind = kind
        self.description = description
        self.header = header
        self.member = member if _GLOBAL_type(member) is _GLOBAL_list else _GLOBAL_list(member)
        self.memberdef = (
            memberdef if _GLOBAL_type(memberdef) is _GLOBAL_list else _GLOBAL_list(memberdef)
        )


class Node_tableofcontentsType(Node):
    __slots__ = (
        "tocsect",
        "tableofcontents",
    )

    _fields = __slots__

    def __init__(
        self,
        tocsect: Iterable[Node_tableofcontentsKindType] = (),
        tableofcontents: Iterable[Node_tableofcontentsType] = (),
    ):  # pragma: no cover
        self.tocsect = tocsect if _GLOBAL_type(tocsect) is _GLOBAL_list else _GLOBAL_list(tocsect)
        self.tableofcontents = (
            tableofcontents
            if _GLOBAL_type(tableofcontents) is _GLOBAL_list
            else _GLOBAL_list(tableofcontents)
        )


if TYPE_CHECKING:
    ListItem_linkedTextType: TypeAlias = Union[str, TaggedValue[Literal["ref"], "Node_refTextType"]]


class Node_linkedTextType(ListNode["ListItem_linkedTextType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_descriptionType: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalType"],
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect1"], "Node_docSect1Type"],
        TaggedValue[Literal["sect2"], "Node_docSect2Type"],
        TaggedValue[Literal["sect3"], "Node_docSect3Type"],
    ]


class Node_descriptionType(ListNode["ListItem_descriptionType"]):
    __slots__ = ("title",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        title: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.title = title


class Node_exportsType(Node):
    __slots__ = ("export",)

    _fields = __slots__

    def __init__(
        self,
        export: Iterable[Node_exportType] = (),
    ):  # pragma: no cover
        self.export = export if _GLOBAL_type(export) is _GLOBAL_list else _GLOBAL_list(export)


class Node_listingType(Node):
    __slots__ = (
        "filename",
        "codeline",
    )

    _fields = __slots__

    def __init__(
        self,
        filename: str | None = None,
        codeline: Iterable[Node_codelineType] = (),
    ):  # pragma: no cover
        self.filename = filename
        self.codeline = (
            codeline if _GLOBAL_type(codeline) is _GLOBAL_list else _GLOBAL_list(codeline)
        )


class Node_locationType(Node):
    __slots__ = (
        "bodyend",
        "bodyfile",
        "bodystart",
        "column",
        "declcolumn",
        "declfile",
        "declline",
        "file",
        "line",
    )

    _fields = __slots__

    def __init__(
        self,
        file: str,
        bodyend: int | None = None,
        bodyfile: str | None = None,
        bodystart: int | None = None,
        column: int | None = None,
        declcolumn: int | None = None,
        declfile: str | None = None,
        declline: int | None = None,
        line: int | None = None,
    ):  # pragma: no cover
        self.bodyend = bodyend
        self.bodyfile = bodyfile
        self.bodystart = bodystart
        self.column = column
        self.declcolumn = declcolumn
        self.declfile = declfile
        self.declline = declline
        self.file = file
        self.line = line


class Node_listofallmembersType(Node):
    __slots__ = ("member",)

    _fields = __slots__

    def __init__(
        self,
        member: Iterable[Node_memberRefType] = (),
    ):  # pragma: no cover
        self.member = member if _GLOBAL_type(member) is _GLOBAL_list else _GLOBAL_list(member)


class Node_memberRefType(Node):
    __slots__ = (
        "ambiguityscope",
        "prot",
        "refid",
        "virt",
        "name",
        "scope",
    )

    _fields = __slots__

    def __init__(
        self,
        prot: DoxProtectionKind,
        refid: str,
        virt: DoxVirtualKind,
        name: str,
        scope: str,
        ambiguityscope: str | None = None,
    ):  # pragma: no cover
        self.ambiguityscope = ambiguityscope
        self.prot = prot
        self.refid = refid
        self.virt = virt
        self.name = name
        self.scope = scope


class Node_memberdefType(Node):
    __slots__ = (
        "accessor",
        "add",
        "attribute",
        "bound",
        "const",
        "constexpr",
        "consteval",
        "constinit",
        "constrained",
        "explicit",
        "extern",
        "final",
        "gettable",
        "id",
        "initonly",
        "inline",
        "kind",
        "maybeambiguous",
        "maybedefault",
        "maybevoid",
        "mutable",
        "new",
        "nodiscard",
        "noexcept",
        "noexceptexpression",
        "optional",
        "privategettable",
        "privatesettable",
        "property",
        "prot",
        "protectedgettable",
        "protectedsettable",
        "raise_",
        "readable",
        "readonly",
        "refqual",
        "removable",
        "remove",
        "required",
        "sealed",
        "settable",
        "static",
        "strong",
        "transient",
        "virt",
        "volatile",
        "writable",
        "argsstring",
        "bitfield",
        "briefdescription",
        "definition",
        "detaileddescription",
        "enumvalue",
        "exceptions",
        "inbodydescription",
        "initializer",
        "location",
        "name",
        "param",
        "qualifiedname",
        "qualifier",
        "read",
        "referencedby",
        "references",
        "reimplementedby",
        "reimplements",
        "requiresclause",
        "templateparamlist",
        "type",
        "write",
    )

    _fields = __slots__

    def __init__(
        self,
        id: str,
        kind: DoxMemberKind,
        prot: DoxProtectionKind,
        static: bool,
        location: Node_locationType,
        name: str,
        accessor: DoxAccessor | None = None,
        add: bool | None = None,
        attribute: bool | None = None,
        bound: bool | None = None,
        const: bool | None = None,
        constexpr: bool | None = None,
        consteval: bool | None = None,
        constinit: bool | None = None,
        constrained: bool | None = None,
        explicit: bool | None = None,
        extern: bool | None = None,
        final: bool | None = None,
        gettable: bool | None = None,
        initonly: bool | None = None,
        inline: bool | None = None,
        maybeambiguous: bool | None = None,
        maybedefault: bool | None = None,
        maybevoid: bool | None = None,
        mutable: bool | None = None,
        new: bool | None = None,
        nodiscard: bool | None = None,
        noexcept: bool | None = None,
        noexceptexpression: str | None = None,
        optional: bool | None = None,
        privategettable: bool | None = None,
        privatesettable: bool | None = None,
        property: bool | None = None,
        protectedgettable: bool | None = None,
        protectedsettable: bool | None = None,
        raise_: bool | None = None,
        readable: bool | None = None,
        readonly: bool | None = None,
        refqual: DoxRefQualifierKind | None = None,
        removable: bool | None = None,
        remove: bool | None = None,
        required: bool | None = None,
        sealed: bool | None = None,
        settable: bool | None = None,
        strong: bool | None = None,
        transient: bool | None = None,
        virt: DoxVirtualKind | None = None,
        volatile: bool | None = None,
        writable: bool | None = None,
        argsstring: str | None = None,
        bitfield: str | None = None,
        briefdescription: Node_descriptionType | None = None,
        definition: str | None = None,
        detaileddescription: Node_descriptionType | None = None,
        enumvalue: Iterable[Node_enumvalueType] = (),
        exceptions: Node_linkedTextType | None = None,
        inbodydescription: Node_descriptionType | None = None,
        initializer: Node_linkedTextType | None = None,
        param: Iterable[Node_paramType] = (),
        qualifiedname: str | None = None,
        qualifier: Iterable[str] = (),
        read: str | None = None,
        referencedby: Iterable[Node_referenceType] = (),
        references: Iterable[Node_referenceType] = (),
        reimplementedby: Iterable[Node_reimplementType] = (),
        reimplements: Iterable[Node_reimplementType] = (),
        requiresclause: Node_linkedTextType | None = None,
        templateparamlist: Node_templateparamlistType | None = None,
        type: Node_linkedTextType | None = None,
        write: str | None = None,
    ):  # pragma: no cover
        self.accessor = accessor
        self.add = add
        self.attribute = attribute
        self.bound = bound
        self.const = const
        self.constexpr = constexpr
        self.consteval = consteval
        self.constinit = constinit
        self.constrained = constrained
        self.explicit = explicit
        self.extern = extern
        self.final = final
        self.gettable = gettable
        self.id = id
        self.initonly = initonly
        self.inline = inline
        self.kind = kind
        self.maybeambiguous = maybeambiguous
        self.maybedefault = maybedefault
        self.maybevoid = maybevoid
        self.mutable = mutable
        self.new = new
        self.nodiscard = nodiscard
        self.noexcept = noexcept
        self.noexceptexpression = noexceptexpression
        self.optional = optional
        self.privategettable = privategettable
        self.privatesettable = privatesettable
        self.property = property
        self.prot = prot
        self.protectedgettable = protectedgettable
        self.protectedsettable = protectedsettable
        self.raise_ = raise_
        self.readable = readable
        self.readonly = readonly
        self.refqual = refqual
        self.removable = removable
        self.remove = remove
        self.required = required
        self.sealed = sealed
        self.settable = settable
        self.static = static
        self.strong = strong
        self.transient = transient
        self.virt = virt
        self.volatile = volatile
        self.writable = writable
        self.argsstring = argsstring
        self.bitfield = bitfield
        self.briefdescription = briefdescription
        self.definition = definition
        self.detaileddescription = detaileddescription
        self.enumvalue = (
            enumvalue if _GLOBAL_type(enumvalue) is _GLOBAL_list else _GLOBAL_list(enumvalue)
        )
        self.exceptions = exceptions
        self.inbodydescription = inbodydescription
        self.initializer = initializer
        self.location = location
        self.name = name
        self.param = param if _GLOBAL_type(param) is _GLOBAL_list else _GLOBAL_list(param)
        self.qualifiedname = qualifiedname
        self.qualifier = (
            qualifier if _GLOBAL_type(qualifier) is _GLOBAL_list else _GLOBAL_list(qualifier)
        )
        self.read = read
        self.referencedby = (
            referencedby
            if _GLOBAL_type(referencedby) is _GLOBAL_list
            else _GLOBAL_list(referencedby)
        )
        self.references = (
            references if _GLOBAL_type(references) is _GLOBAL_list else _GLOBAL_list(references)
        )
        self.reimplementedby = (
            reimplementedby
            if _GLOBAL_type(reimplementedby) is _GLOBAL_list
            else _GLOBAL_list(reimplementedby)
        )
        self.reimplements = (
            reimplements
            if _GLOBAL_type(reimplements) is _GLOBAL_list
            else _GLOBAL_list(reimplements)
        )
        self.requiresclause = requiresclause
        self.templateparamlist = templateparamlist
        self.type = type
        self.write = write


class Node_MemberType(Node):
    __slots__ = (
        "kind",
        "refid",
        "name",
    )

    _fields = __slots__

    def __init__(
        self,
        kind: MemberKind,
        refid: str,
        name: str,
    ):  # pragma: no cover
        self.kind = kind
        self.refid = refid
        self.name = name


class Node_paramType(Node):
    __slots__ = (
        "array",
        "attributes",
        "briefdescription",
        "declname",
        "defname",
        "defval",
        "type",
        "typeconstraint",
    )

    _fields = __slots__

    def __init__(
        self,
        array: str | None = None,
        attributes: str | None = None,
        briefdescription: Node_descriptionType | None = None,
        declname: str | None = None,
        defname: str | None = None,
        defval: Node_linkedTextType | None = None,
        type: Node_linkedTextType | None = None,
        typeconstraint: Node_linkedTextType | None = None,
    ):  # pragma: no cover
        self.array = array
        self.attributes = attributes
        self.briefdescription = briefdescription
        self.declname = declname
        self.defname = defname
        self.defval = defval
        self.type = type
        self.typeconstraint = typeconstraint


class Node_enumvalueType(Node):
    __slots__ = (
        "id",
        "prot",
        "briefdescription",
        "detaileddescription",
        "initializer",
        "name",
    )

    _fields = __slots__

    def __init__(
        self,
        id: str,
        prot: DoxProtectionKind,
        name: str,
        briefdescription: Node_descriptionType | None = None,
        detaileddescription: Node_descriptionType | None = None,
        initializer: Node_linkedTextType | None = None,
    ):  # pragma: no cover
        self.id = id
        self.prot = prot
        self.briefdescription = briefdescription
        self.detaileddescription = detaileddescription
        self.initializer = initializer
        self.name = name


class Node_referenceType(ListNode["str"]):
    __slots__ = (
        "compoundref",
        "endline",
        "refid",
        "startline",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        refid: str,
        compoundref: str | None = None,
        endline: int | None = None,
        startline: int | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.compoundref = compoundref
        self.endline = endline
        self.refid = refid
        self.startline = startline


if TYPE_CHECKING:
    ListItem_docInternalType: TypeAlias = Union[
        str,
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect1"], "Node_docSect1Type"],
    ]


class Node_docInternalType(ListNode["ListItem_docInternalType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docSect1Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalS1Type"],
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect2"], "Node_docSect2Type"],
        TaggedValue[Literal["sect3"], "Node_docSect3Type"],
    ]


class Node_docSect1Type(ListNode["ListItem_docSect1Type"]):
    __slots__ = (
        "id",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
        title: Node_docTitleType | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id
        self.title = title


class Node_nodeType(Node):
    __slots__ = (
        "id",
        "childnode",
        "label",
        "link",
    )

    _fields = __slots__

    def __init__(
        self,
        id: str,
        label: str,
        childnode: Iterable[Node_childnodeType] = (),
        link: Node_linkType | None = None,
    ):  # pragma: no cover
        self.id = id
        self.childnode = (
            childnode if _GLOBAL_type(childnode) is _GLOBAL_list else _GLOBAL_list(childnode)
        )
        self.label = label
        self.link = link


class Node_linkType(Node):
    __slots__ = (
        "external",
        "refid",
    )

    _fields = __slots__

    def __init__(
        self,
        refid: str,
        external: str | None = None,
    ):  # pragma: no cover
        self.external = external
        self.refid = refid


class Node_childnodeType(Node):
    __slots__ = (
        "refid",
        "relation",
        "edgelabel",
    )

    _fields = __slots__

    def __init__(
        self,
        refid: str,
        relation: DoxGraphRelation,
        edgelabel: Iterable[str] = (),
    ):  # pragma: no cover
        self.refid = refid
        self.relation = relation
        self.edgelabel = (
            edgelabel if _GLOBAL_type(edgelabel) is _GLOBAL_list else _GLOBAL_list(edgelabel)
        )


class Node_codelineType(Node):
    __slots__ = (
        "external",
        "lineno",
        "refid",
        "refkind",
        "highlight",
    )

    _fields = __slots__

    def __init__(
        self,
        external: bool | None = None,
        lineno: int | None = None,
        refid: str | None = None,
        refkind: DoxRefKind | None = None,
        highlight: Iterable[Node_highlightType] = (),
    ):  # pragma: no cover
        self.external = external
        self.lineno = lineno
        self.refid = refid
        self.refkind = refkind
        self.highlight = (
            highlight if _GLOBAL_type(highlight) is _GLOBAL_list else _GLOBAL_list(highlight)
        )


if TYPE_CHECKING:
    ListItem_highlightType: TypeAlias = Union[str, TaggedValue[Literal["ref"], "Node_refTextType"]]


class Node_highlightType(ListNode["ListItem_highlightType"]):
    __slots__ = ("class_",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        class_: DoxHighlightClass,
    ):  # pragma: no cover
        super().__init__(__children)
        self.class_ = class_


if TYPE_CHECKING:
    ListItem_docInternalS1Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect2"], "Node_docSect2Type"],
    ]


class Node_docInternalS1Type(ListNode["ListItem_docInternalS1Type"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docSect2Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalS2Type"],
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect3"], "Node_docSect3Type"],
    ]


class Node_docSect2Type(ListNode["ListItem_docSect2Type"]):
    __slots__ = (
        "id",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
        title: Node_docTitleType | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id
        self.title = title


if TYPE_CHECKING:
    ListItem_docSect3Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalS3Type"],
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect4"], "Node_docSect4Type"],
    ]


class Node_docSect3Type(ListNode["ListItem_docSect3Type"]):
    __slots__ = (
        "id",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
        title: Node_docTitleType | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id
        self.title = title


if TYPE_CHECKING:
    ListItem_docInternalS2Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect3"], "Node_docSect3Type"],
    ]


class Node_docInternalS2Type(ListNode["ListItem_docInternalS2Type"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docSect4Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalS4Type"],
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect5"], "Node_docSect5Type"],
    ]


class Node_docSect4Type(ListNode["ListItem_docSect4Type"]):
    __slots__ = (
        "id",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
        title: Node_docTitleType | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id
        self.title = title


if TYPE_CHECKING:
    ListItem_docSect5Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalS5Type"],
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect6"], "Node_docSect6Type"],
    ]


class Node_docSect5Type(ListNode["ListItem_docSect5Type"]):
    __slots__ = (
        "id",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
        title: Node_docTitleType | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id
        self.title = title


if TYPE_CHECKING:
    ListItem_docSect6Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["internal"], "Node_docInternalS6Type"],
        TaggedValue[Literal["para"], "Node_docParaType"],
    ]


class Node_docSect6Type(ListNode["ListItem_docSect6Type"]):
    __slots__ = (
        "id",
        "title",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
        title: Node_docTitleType | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id
        self.title = title


if TYPE_CHECKING:
    ListItem_docInternalS3Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect3"], "Node_docSect4Type"],
    ]


class Node_docInternalS3Type(ListNode["ListItem_docInternalS3Type"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docInternalS4Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect5"], "Node_docSect4Type"],
    ]


class Node_docInternalS4Type(ListNode["ListItem_docInternalS4Type"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docInternalS5Type: TypeAlias = Union[
        str,
        TaggedValue[Literal["para"], "Node_docParaType"],
        TaggedValue[Literal["sect6"], "Node_docSect4Type"],
    ]


class Node_docInternalS5Type(ListNode["ListItem_docInternalS5Type"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docInternalS6Type: TypeAlias = Union[
        str, TaggedValue[Literal["para"], "Node_docParaType"]
    ]


class Node_docInternalS6Type(ListNode["ListItem_docInternalS6Type"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docListItemType(ListNode["Node_docParaType"]):
    __slots__ = (
        "override",
        "value",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        override: DoxCheck | None = None,
        value: int | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.override = override
        self.value = value


class Node_docCaptionType(ListNode["str"]):
    __slots__ = ("id",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id


class Node_docRowType(Node):
    __slots__ = ("entry",)

    _fields = __slots__

    def __init__(
        self,
        entry: Iterable[Node_docEntryType] = (),
    ):  # pragma: no cover
        self.entry = entry if _GLOBAL_type(entry) is _GLOBAL_list else _GLOBAL_list(entry)


class Node_docEntryType(Node):
    __slots__ = (
        "align",
        "class_",
        "colspan",
        "rowspan",
        "thead",
        "valign",
        "width",
        "para",
    )

    _fields = __slots__

    def __init__(
        self,
        thead: bool,
        align: DoxAlign | None = None,
        class_: str | None = None,
        colspan: int | None = None,
        rowspan: int | None = None,
        valign: DoxVerticalAlign | None = None,
        width: str | None = None,
        para: Iterable[Node_docParaType] = (),
    ):  # pragma: no cover
        self.align = align
        self.class_ = class_
        self.colspan = colspan
        self.rowspan = rowspan
        self.thead = thead
        self.valign = valign
        self.width = width
        self.para = para if _GLOBAL_type(para) is _GLOBAL_list else _GLOBAL_list(para)


class Node_docTocItemType(ListNode["str"]):
    __slots__ = ("id",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id


class Node_docParamListItem(Node):
    __slots__ = (
        "parameterdescription",
        "parameternamelist",
    )

    _fields = __slots__

    def __init__(
        self,
        parameterdescription: Node_descriptionType,
        parameternamelist: Iterable[Node_docParamNameList] = (),
    ):  # pragma: no cover
        self.parameterdescription = parameterdescription
        self.parameternamelist = (
            parameternamelist
            if _GLOBAL_type(parameternamelist) is _GLOBAL_list
            else _GLOBAL_list(parameternamelist)
        )


class Node_docParamNameList(Node):
    __slots__ = (
        "parametername",
        "parametertype",
    )

    _fields = __slots__

    def __init__(
        self,
        parametername: Iterable[Node_docParamName] = (),
        parametertype: Iterable[Node_docParamType] = (),
    ):  # pragma: no cover
        self.parametername = (
            parametername
            if _GLOBAL_type(parametername) is _GLOBAL_list
            else _GLOBAL_list(parametername)
        )
        self.parametertype = (
            parametertype
            if _GLOBAL_type(parametertype) is _GLOBAL_list
            else _GLOBAL_list(parametertype)
        )


if TYPE_CHECKING:
    ListItem_docParamType: TypeAlias = Union[str, TaggedValue[Literal["ref"], "Node_refTextType"]]


class Node_docParamType(ListNode["ListItem_docParamType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


if TYPE_CHECKING:
    ListItem_docParamName: TypeAlias = Union[str, TaggedValue[Literal["ref"], "Node_refTextType"]]


class Node_docParamName(ListNode["ListItem_docParamName"]):
    __slots__ = ("direction",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        direction: DoxParamDir | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.direction = direction


class Node_tableofcontentsKindType(Node):
    __slots__ = (
        "name",
        "reference",
        "tableofcontents",
    )

    _fields = __slots__

    def __init__(
        self,
        name: str,
        reference: str,
        tableofcontents: Iterable[Node_tableofcontentsType] = (),
    ):  # pragma: no cover
        self.name = name
        self.reference = reference
        self.tableofcontents = (
            tableofcontents
            if _GLOBAL_type(tableofcontents) is _GLOBAL_list
            else _GLOBAL_list(tableofcontents)
        )


class Node_incType(ListNode["str"]):
    __slots__ = (
        "refid",
        "local",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        local: bool,
        refid: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid
        self.local = local


class Node_compoundRefType(ListNode["str"]):
    __slots__ = (
        "refid",
        "prot",
        "virt",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        prot: DoxProtectionKind,
        virt: DoxVirtualKind,
        refid: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid
        self.prot = prot
        self.virt = virt


class Node_refType(ListNode["str"]):
    __slots__ = (
        "refid",
        "prot",
        "inline",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        refid: str,
        prot: DoxProtectionKind | None = None,
        inline: bool | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid
        self.prot = prot
        self.inline = inline


class Node_exportType(ListNode["str"]):
    __slots__ = ("refid",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        refid: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid


class Node_refTextType(ListNode["str"]):
    __slots__ = (
        "refid",
        "kindref",
        "external",
        "tooltip",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        refid: str,
        kindref: DoxRefKind,
        external: str | None = None,
        tooltip: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid
        self.kindref = kindref
        self.external = external
        self.tooltip = tooltip


class Node_reimplementType(ListNode["str"]):
    __slots__ = ("refid",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        refid: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid


class DoxRefKind(enum.Enum):
    compound = "compound"
    member = "member"


class MemberKind(enum.Enum):
    define = "define"
    property = "property"
    event = "event"
    variable = "variable"
    typedef = "typedef"
    enum = "enum"
    enumvalue = "enumvalue"
    function = "function"
    signal = "signal"
    prototype = "prototype"
    friend = "friend"
    dcop = "dcop"
    slot = "slot"


class DoxMemberKind(enum.Enum):
    define = "define"
    property = "property"
    event = "event"
    variable = "variable"
    typedef = "typedef"
    enum = "enum"
    function = "function"
    signal = "signal"
    prototype = "prototype"
    friend = "friend"
    dcop = "dcop"
    slot = "slot"
    interface = "interface"
    service = "service"


if TYPE_CHECKING:
    ListItem_docTitleCmdGroup: TypeAlias = Union[
        str,
        TaggedValue[Literal["anchor"], "Node_docAnchorType"],
        TaggedValue[
            Literal[
                "bold",
                "center",
                "cite",
                "computeroutput",
                "del",
                "emphasis",
                "ins",
                "s",
                "small",
                "strike",
                "subscript",
                "superscript",
                "underline",
            ],
            "Node_docMarkupType",
        ],
        TaggedValue[Literal["docbookonly", "latexonly", "manonly", "rtfonly", "xmlonly"], "str"],
        TaggedValue[Literal["dot", "msc"], "Node_docDotMscType"],
        TaggedValue[Literal["emoji"], "Node_docEmojiType"],
        TaggedValue[Literal["formula"], "Node_docFormulaType"],
        TaggedValue[Literal["htmlonly"], "Node_docHtmlOnlyType"],
        TaggedValue[Literal["image"], "Node_docImageType"],
        TaggedValue[Literal["plantuml"], "Node_docPlantumlType"],
        TaggedValue[Literal["ref"], "Node_docRefTextType"],
        TaggedValue[Literal["ulink"], "Node_docURLLink"],
    ]

if TYPE_CHECKING:
    ListItem_docCmdGroup: TypeAlias = Union[
        ListItem_docTitleCmdGroup,
        TaggedValue[Literal["blockquote"], "Node_docBlockQuoteType"],
        TaggedValue[Literal["copydoc"], "Node_docCopyType"],
        TaggedValue[Literal["details"], "Node_docDetailsType"],
        TaggedValue[
            Literal["diafile", "dotfile", "mscfile", "plantumlfile"], "Node_docImageFileType"
        ],
        TaggedValue[Literal["heading"], "Node_docHeadingType"],
        TaggedValue[Literal["hruler"], "None"],
        TaggedValue[Literal["indexentry"], "Node_docIndexEntryType"],
        TaggedValue[Literal["itemizedlist", "orderedlist"], "Node_docListType"],
        TaggedValue[Literal["javadoccode", "javadocliteral", "verbatim"], "str"],
        TaggedValue[Literal["language"], "Node_docLanguageType"],
        TaggedValue[Literal["parameterlist"], "Node_docParamListType"],
        TaggedValue[Literal["parblock"], "Node_docParBlockType"],
        TaggedValue[Literal["preformatted"], "Node_docMarkupType"],
        TaggedValue[Literal["programlisting"], "Node_listingType"],
        TaggedValue[Literal["simplesect"], "Node_docSimpleSectType"],
        TaggedValue[Literal["table"], "Node_docTableType"],
        TaggedValue[Literal["title"], "Node_docTitleType"],
        TaggedValue[Literal["toclist"], "Node_docTocListType"],
        TaggedValue[Literal["variablelist"], "Node_docVariableListType"],
        TaggedValue[Literal["xrefsect"], "Node_docXRefSectType"],
    ]


class Node_docParaType(ListNode["ListItem_docCmdGroup"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docMarkupType(ListNode["ListItem_docCmdGroup"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docTitleType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docSummaryType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docURLLink(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = ("url",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        url: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.url = url


class Node_docHtmlOnlyType(ListNode["str"]):
    __slots__ = ("block",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        block: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.block = block


class Node_docImageType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = (
        "type",
        "name",
        "width",
        "height",
        "alt",
        "inline",
        "caption",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        type: DoxImageKind | None = None,
        name: str | None = None,
        width: str | None = None,
        height: str | None = None,
        alt: str | None = None,
        inline: bool | None = None,
        caption: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.type = type
        self.name = name
        self.width = width
        self.height = height
        self.alt = alt
        self.inline = inline
        self.caption = caption


class Node_docDotMscType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = (
        "name",
        "width",
        "height",
        "caption",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        name: str | None = None,
        width: str | None = None,
        height: str | None = None,
        caption: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.name = name
        self.width = width
        self.height = height
        self.caption = caption


class Node_docPlantumlType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = (
        "name",
        "width",
        "height",
        "caption",
        "engine",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        name: str | None = None,
        width: str | None = None,
        height: str | None = None,
        caption: str | None = None,
        engine: DoxPlantumlEngine | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.name = name
        self.width = width
        self.height = height
        self.caption = caption
        self.engine = engine


class Node_docRefTextType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = (
        "refid",
        "kindref",
        "external",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        refid: str,
        kindref: str,
        external: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.refid = refid
        self.kindref = kindref
        self.external = external


class Node_docHeadingType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = ("level",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        level: int,
    ):  # pragma: no cover
        super().__init__(__children)
        self.level = level


class Node_docImageFileType(ListNode["ListItem_docTitleCmdGroup"]):
    __slots__ = (
        "name",
        "width",
        "height",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        name: str | None = None,
        width: str | None = None,
        height: str | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.name = name
        self.width = width
        self.height = height


class Node_docAnchorType(ListNode["str"]):
    __slots__ = ("id",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id


class Node_docFormulaType(ListNode["str"]):
    __slots__ = ("id",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        id: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.id = id


class Node_docEmojiType(Node):
    __slots__ = (
        "name",
        "unicode",
    )

    _fields = __slots__

    def __init__(
        self,
        name: str,
        unicode: str,
    ):  # pragma: no cover
        self.name = name
        self.unicode = unicode


class Node_docIndexEntryType(Node):
    __slots__ = (
        "primaryie",
        "secondaryie",
    )

    _fields = __slots__

    def __init__(
        self,
        primaryie: str,
        secondaryie: str,
    ):  # pragma: no cover
        self.primaryie = primaryie
        self.secondaryie = secondaryie


class Node_docListType(ListNode["Node_docListItemType"]):
    __slots__ = (
        "type",
        "start",
    )

    _fields = __slots__

    def __init__(
        self,
        __children,
        type: DoxOlType | None = None,
        start: int | None = None,
    ):  # pragma: no cover
        super().__init__(__children)
        self.type = type
        self.start = start


class Node_docSimpleSectType(Node):
    __slots__ = (
        "kind",
        "title",
        "para",
    )

    _fields = __slots__

    def __init__(
        self,
        kind: DoxSimpleSectKind,
        title: Node_docTitleType | None = None,
        para: Iterable[Node_docParaType] = (),
    ):  # pragma: no cover
        self.kind = kind
        self.title = title
        self.para = para if _GLOBAL_type(para) is _GLOBAL_list else _GLOBAL_list(para)


class ListItem_docVariableListType(NamedTuple):
    varlistentry: Node_docVarListEntryType
    listitem: Node_docListItemType


class Node_docVariableListType(ListNode["ListItem_docVariableListType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docTableType(Node):
    __slots__ = (
        "rows",
        "cols",
        "width",
        "caption",
        "row",
    )

    _fields = __slots__

    def __init__(
        self,
        rows: int,
        cols: int,
        width: str | None = None,
        caption: Node_docCaptionType | None = None,
        row: Iterable[Node_docRowType] = (),
    ):  # pragma: no cover
        self.rows = rows
        self.cols = cols
        self.width = width
        self.caption = caption
        self.row = row if _GLOBAL_type(row) is _GLOBAL_list else _GLOBAL_list(row)


class Node_docTocListType(ListNode["Node_docTocItemType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docLanguageType(ListNode["Node_docParaType"]):
    __slots__ = ("langid",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        langid: str,
    ):  # pragma: no cover
        super().__init__(__children)
        self.langid = langid


class Node_docParamListType(ListNode["Node_docParamListItem"]):
    __slots__ = ("kind",)

    _fields = __slots__

    def __init__(
        self,
        __children,
        kind: DoxParamListKind,
    ):  # pragma: no cover
        super().__init__(__children)
        self.kind = kind


class Node_docXRefSectType(Node):
    __slots__ = (
        "id",
        "xreftitle",
        "xrefdescription",
    )

    _fields = __slots__

    def __init__(
        self,
        id: str,
        xrefdescription: Node_descriptionType,
        xreftitle: Iterable[str] = (),
    ):  # pragma: no cover
        self.id = id
        self.xreftitle = (
            xreftitle if _GLOBAL_type(xreftitle) is _GLOBAL_list else _GLOBAL_list(xreftitle)
        )
        self.xrefdescription = xrefdescription


class Node_docCopyType(Node):
    __slots__ = (
        "link",
        "para",
        "sec1",
        "internal",
    )

    _fields = __slots__

    def __init__(
        self,
        link: str,
        para: Iterable[Node_docParaType] = (),
        sec1: Iterable[Node_docSect1Type] = (),
        internal: Node_docInternalType | None = None,
    ):  # pragma: no cover
        self.link = link
        self.para = para if _GLOBAL_type(para) is _GLOBAL_list else _GLOBAL_list(para)
        self.sec1 = sec1 if _GLOBAL_type(sec1) is _GLOBAL_list else _GLOBAL_list(sec1)
        self.internal = internal


class Node_docDetailsType(Node):
    __slots__ = (
        "summary",
        "para",
    )

    _fields = __slots__

    def __init__(
        self,
        summary: Node_docSummaryType | None = None,
        para: Iterable[Node_docParaType] = (),
    ):  # pragma: no cover
        self.summary = summary
        self.para = para if _GLOBAL_type(para) is _GLOBAL_list else _GLOBAL_list(para)


class Node_docBlockQuoteType(ListNode["Node_docParaType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docParBlockType(ListNode["Node_docParaType"]):
    __slots__ = ()

    _fields = __slots__

    def __init__(
        self,
        __children,
    ):  # pragma: no cover
        super().__init__(__children)


class Node_docVarListEntryType(Node):
    __slots__ = ("term",)

    _fields = __slots__

    def __init__(
        self,
        term: Node_docTitleType,
    ):  # pragma: no cover
        self.term = term


class DoxCompoundKind(enum.Enum):
    class_ = "class"
    struct = "struct"
    union = "union"
    interface = "interface"
    protocol = "protocol"
    category = "category"
    exception = "exception"
    service = "service"
    singleton = "singleton"
    module = "module"
    type = "type"
    file = "file"
    namespace = "namespace"
    group = "group"
    page = "page"
    example = "example"
    dir = "dir"
    concept = "concept"


class CompoundKind(enum.Enum):
    class_ = "class"
    struct = "struct"
    union = "union"
    interface = "interface"
    protocol = "protocol"
    category = "category"
    exception = "exception"
    module = "module"
    type = "type"
    file = "file"
    namespace = "namespace"
    group = "group"
    page = "page"
    example = "example"
    dir = "dir"
    concept = "concept"


class DoxLanguage(enum.Enum):
    Unknown = "Unknown"
    IDL = "IDL"
    Java = "Java"
    CSharp = "C#"
    D = "D"
    PHP = "PHP"
    Objective_C = "Objective-C"
    CPlusPlus = "C++"
    JavaScript = "JavaScript"
    Python = "Python"
    Fortran = "Fortran"
    VHDL = "VHDL"
    XML = "XML"
    SQL = "SQL"
    Markdown = "Markdown"
    Slice = "Slice"
    Lex = "Lex"


class DoxProtectionKind(enum.Enum):
    public = "public"
    protected = "protected"
    private = "private"
    package = "package"


class DoxRefQualifierKind(enum.Enum):
    lvalue = "lvalue"
    rvalue = "rvalue"


class DoxVirtualKind(enum.Enum):
    non_virtual = "non-virtual"
    virtual = "virtual"
    pure_virtual = "pure-virtual"


class DoxSectionKind(enum.Enum):
    user_defined = "user-defined"
    public_type = "public-type"
    public_func = "public-func"
    public_attrib = "public-attrib"
    public_slot = "public-slot"
    signal = "signal"
    dcop_func = "dcop-func"
    property = "property"
    event = "event"
    public_static_func = "public-static-func"
    public_static_attrib = "public-static-attrib"
    protected_type = "protected-type"
    protected_func = "protected-func"
    protected_attrib = "protected-attrib"
    protected_slot = "protected-slot"
    protected_static_func = "protected-static-func"
    protected_static_attrib = "protected-static-attrib"
    package_type = "package-type"
    package_func = "package-func"
    package_attrib = "package-attrib"
    package_static_func = "package-static-func"
    package_static_attrib = "package-static-attrib"
    private_type = "private-type"
    private_func = "private-func"
    private_attrib = "private-attrib"
    private_slot = "private-slot"
    private_static_func = "private-static-func"
    private_static_attrib = "private-static-attrib"
    friend = "friend"
    related = "related"
    define = "define"
    prototype = "prototype"
    typedef = "typedef"
    enum = "enum"
    func = "func"
    var = "var"


class DoxHighlightClass(enum.Enum):
    comment = "comment"
    normal = "normal"
    preprocessor = "preprocessor"
    keyword = "keyword"
    keywordtype = "keywordtype"
    keywordflow = "keywordflow"
    stringliteral = "stringliteral"
    xmlcdata = "xmlcdata"
    charliteral = "charliteral"
    vhdlkeyword = "vhdlkeyword"
    vhdllogic = "vhdllogic"
    vhdlchar = "vhdlchar"
    vhdldigit = "vhdldigit"


class DoxSimpleSectKind(enum.Enum):
    see = "see"
    return_ = "return"
    author = "author"
    authors = "authors"
    version = "version"
    since = "since"
    date = "date"
    note = "note"
    warning = "warning"
    pre = "pre"
    post = "post"
    copyright = "copyright"
    invariant = "invariant"
    remark = "remark"
    attention = "attention"
    important = "important"
    par = "par"
    rcs = "rcs"


class DoxImageKind(enum.Enum):
    html = "html"
    latex = "latex"
    docbook = "docbook"
    rtf = "rtf"
    xml = "xml"


class DoxPlantumlEngine(enum.Enum):
    uml = "uml"
    bpm = "bpm"
    wire = "wire"
    dot = "dot"
    ditaa = "ditaa"
    salt = "salt"
    math = "math"
    latex = "latex"
    gantt = "gantt"
    mindmap = "mindmap"
    wbs = "wbs"
    yaml = "yaml"
    creole = "creole"
    json = "json"
    flow = "flow"
    board = "board"
    git = "git"
    hcl = "hcl"
    regex = "regex"
    ebnf = "ebnf"
    files = "files"


class DoxParamListKind(enum.Enum):
    param = "param"
    retval = "retval"
    exception = "exception"
    templateparam = "templateparam"


class DoxParamDir(enum.Enum):
    in_ = "in"
    out = "out"
    inout = "inout"


class DoxAccessor(enum.Enum):
    retain = "retain"
    copy = "copy"
    assign = "assign"
    weak = "weak"
    strong = "strong"
    unretained = "unretained"


class DoxAlign(enum.Enum):
    left = "left"
    right = "right"
    center = "center"


class DoxVerticalAlign(enum.Enum):
    bottom = "bottom"
    top = "top"
    middle = "middle"


class DoxGraphRelation(enum.Enum):
    include = "include"
    usage = "usage"
    template_instance = "template-instance"
    public_inheritance = "public-inheritance"
    protected_inheritance = "protected-inheritance"
    private_inheritance = "private-inheritance"
    type_constraint = "type-constraint"


class DoxCheck(enum.Enum):
    checked = "checked"
    unchecked = "unchecked"


DoxOlType = Literal["1", "a", "A", "i", "I"]


def parse_str(data: str, /):
    return _parse(data, expat.XMLParserType.Parse)


def parse_file(file, /):
    return _parse(file, expat.XMLParserType.ParseFile)


if TYPE_CHECKING:
    _ChildStartCallback = Callable[["_ParseState", Any, Iterable[tuple[str, str]]], None]
    _FinishCallback = Callable[["_ParseState"], None]
    _TextCallback = Callable[["_ParseState", str], None]
    _Setter = Callable[[Any], None]

    _T_covar = TypeVar("_T_covar", covariant=True)
    _U_covar = TypeVar("_U_covar", covariant=True)


_GLOBAL_type = type
_GLOBAL_list = list


class _ParseCallbacks:
    __slots__ = "value", "setter", "cs_call", "f_call", "t_call"

    value: Any
    """The value corresponding the currently visited XML element."""

    setter: _Setter | None
    """A callback given by the parent element to consume the value.
    
    This may be None if no action is needed.
    """

    cs_call: dict[str, _ChildStartCallback] | None
    """A mapping of element names to callbacks for a children of the current
    element.
    
    This may be None if no child elements are allowed.
    """

    f_call: _FinishCallback | None
    """A callback for when the current element is closed.
    
    This may be None if no action is needed.
    """

    t_call: _TextCallback | None
    """A callback for text contained directly inside the current element.
    
    This may be None if text is not allowed. If None, whitespace is ignored.
    """

    def __init__(self, value=None, setter=None, cs_call=None, f_call=None, t_call=None):
        self.value = value
        self.setter = setter
        self.cs_call = cs_call
        self.f_call = f_call
        self.t_call = t_call


class _ParseState:
    def __init__(self, parser: expat.XMLParserType, /):
        self.parser = parser
        self.parse_callbacks: list[_ParseCallbacks] = []

        # While this is greater than zero all XML content is ignored.
        #
        # This starts at zero. When an unexpected element start is encountered,
        # a warning is issued (via PyErr_WarnFormat) and this is set to 1. Any
        # subsequent element-starts increment this and element-ends decrement
        # this until this is zero again, and normal parsing resumes.
        self.ignore_level: int = 0

    def start_element(self, name: str, attrs: dict[str, str], /) -> None:
        if self.ignore_level:
            self.ignore_level += 1
            return

        cb = self.parse_callbacks[-1]

        if cb.cs_call is not None:
            handler = cb.cs_call.get(name)
            if handler is not None:
                handler(self, cb.value, attrs.items())
                return

        self.set_parse_warning(f'unexpected element "{name}"')

        self.ignore_level = 1

    def end_element(self, unused, /) -> None:
        if self.ignore_level:
            self.ignore_level -= 1
            return

        cb = self.parse_callbacks[-1]

        if cb.f_call is not None:
            cb.f_call(self)

        if cb.setter is not None:
            cb.setter(cb.value)
        self.parse_callbacks.pop()

    def character_data(self, s: str, /) -> None:
        if self.ignore_level:
            return

        cb = self.parse_callbacks[-1]

        if cb.t_call is not None:
            cb.t_call(self, s)
        elif s and not s.isspace():
            self.set_parse_warning("unexpected character data")

    def raise_parse_error(self, msg, /) -> NoReturn:
        raise ParseError(msg, self.parser.CurrentLineNumber)

    def set_parse_warning(self, msg, /) -> None:
        warnings.warn(ParseWarning(f"Warning on line {self.parser.CurrentLineNumber}: {msg}"))


def _node_list_common_text(state: _ParseState, data: str, /):
    value = state.parse_callbacks[-1].value

    if value and type(value[-1]) is str:
        value[-1] += data
    else:
        value.append(data)


def _push_tuple_item(state: _ParseState, tuple_i: int, tag_names: Sequence[str], cls, obj, /):
    if tuple_i == 0:
        if len(obj):
            tuple_size = len(tag_names)
            if len(obj[-1]) < tuple_size:
                state.raise_parse_error(
                    f'"{tag_names[0]}" element can only come after "{tag_names[tuple_size - 1]}" element or be the first in its group',
                )

            obj[-1] = cls._make(obj[-1])

        # tuples are immutable so a list is used while collecting the values
        new_tuple: list[Any] = []
        obj.append(new_tuple)

        return new_tuple.append

    if not obj or len(obj[-1]) < tuple_i:
        state.raise_parse_error(
            f'"{tag_names[tuple_i]}" element can only come after "{tag_names[tuple_i - 1]}" element'
        )

    return obj[-1].append


def _check_complete_tuple(state: _ParseState, tag_names: Sequence[str], cls, obj, /):
    if obj:
        last = obj[-1]

        if len(last) != len(tag_names):
            state.raise_parse_error(
                f'"{tag_names[len(last)]}" element must come after "{tag_names[len(last) - 1]}" element'
            )

        obj[-1] = cls._make(last)


def _warn_unexpected_attribute(state: _ParseState, name: str, /):
    state.set_parse_warning(f'unexpected attribute "{name}"')


def _raise_missing_attribute_error(state: _ParseState, name: str, /):
    state.raise_parse_error(f'missing "{name}" attribute')


def _raise_duplicate_element_error(state: _ParseState, name: str, /):
    state.raise_parse_error(f'"{name}" cannot appear more than once in this context')


def _raise_missing_element_error(state: _ParseState, parent: Any, name: str, /):
    state.raise_parse_error(f'"{parent}" missing "{name}" child')


def _raise_empty_list_element_error(state: _ParseState, name: str, /):
    state.raise_parse_error(f'at least one "{name}" child is required')


def _raise_invalid_int_error(state: _ParseState, value: str, /):
    state.raise_parse_error(f'"{value}" is not a valid integer')


def _raise_invalid_enum_error(state: _ParseState, value: str, /):
    state.raise_parse_error(f'"{value}" is not one of the allowed enumeration values')


def _raise_invalid_char_enum_error(state: _ParseState, c: str, allowed: str, /):
    state.raise_parse_error(
        f'"{c}" is not one of the allowed character values; must be one of "{allowed}"'
    )


def _parse_DoxBool_attribute(state: _ParseState, name: str, value: str, /) -> bool:
    if value == "yes":
        return True
    if value == "no":
        return False

    state.raise_parse_error(f'"{name}" must be "yes" or "no"')


def _node_string_text(state: _ParseState, data: str) -> None:
    state.parse_callbacks[-1].value += data


def _node_start_string(state: _ParseState, setter: _Setter, attr: Iterable[tuple[str, str]], /):
    for name, _ in attr:
        _warn_unexpected_attribute(state, name)

    state.parse_callbacks.append(_ParseCallbacks("", setter, None, None, _node_string_text))


def _node_start_empty(state: _ParseState, setter: _Setter, attr: Iterable[tuple[str, str]], /):
    for name, _ in attr:
        _warn_unexpected_attribute(state, name)

    setter(None)
    state.parse_callbacks.append(_ParseCallbacks())


def _node_start_spType(state: _ParseState, attr: Iterable[tuple[str, str]], /):
    c = " "

    for name, value in attr:
        if name != "value":
            _warn_unexpected_attribute(state, name)

        try:
            c_i = int(value, 10)
        except ValueError:
            state.raise_parse_error('"value" must be a valid integer')
        if 0 > c_i > 127:
            state.raise_parse_error('"value" must be between 0 and 127')

        c = chr(c_i)

    state.parse_callbacks.append(_ParseCallbacks())
    return c


def _node_start_const_char(state: _ParseState, attr: Iterable[tuple[str, str]], /):
    for name, _ in attr:
        _warn_unexpected_attribute(state, name)

    state.parse_callbacks.append(_ParseCallbacks())


def _union_codepoint_element(c):
    def inner(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /) -> None:
        if obj and type(obj[-1]) is str:
            obj[-1] += c
        else:
            obj.append(c)

        _node_start_const_char(state, attr)

    return inner


_cur_list: dict[str, Callable]


def _add_to_list(name):
    def inner(f):
        global _cur_list
        _cur_list[name] = f

    return inner


_node_class_attr__DoxygenTypeIndex = _cur_list = {}


@_add_to_list("version")
def _a__DoxygenTypeIndex__version(state: _ParseState, obj, value: str, /):
    obj.version = value


def _node_class_attr_end__DoxygenTypeIndex(state: _ParseState, obj, /):
    if not hasattr(obj, "version"):
        _raise_missing_attribute_error(state, "version")


_node_class_child__DoxygenTypeIndex = _cur_list = {}


@_add_to_list("compound")
def _e__DoxygenTypeIndex__compound(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__CompoundType(state, obj.compound.append, attr)


def _node_class_start__DoxygenTypeIndex(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_DoxygenTypeIndex.__new__(Node_DoxygenTypeIndex)

    n.compound = []
    for name, value in attr:
        handler = _node_class_attr__DoxygenTypeIndex.get(name)

        if handler is not None:
            handler(state, n, value)
    _node_class_attr_end__DoxygenTypeIndex(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__DoxygenTypeIndex,
            None,
            None,
        )
    )


_node_class_attr__CompoundType = _cur_list = {}


@_add_to_list("refid")
def _a__CompoundType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("kind")
def _a__CompoundType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = CompoundKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__CompoundType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")


_node_class_child__CompoundType = _cur_list = {}


@_add_to_list("name")
def _e__CompoundType__name(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "name"):
        _raise_duplicate_element_error(state, "name")

    _node_start_string(state, functools.partial(setattr, obj, "name"), attr)


@_add_to_list("member")
def _e__CompoundType__member(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__MemberType(state, obj.member.append, attr)


def _node_class_start__CompoundType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_CompoundType.__new__(Node_CompoundType)

    n.member = []
    for name, value in attr:
        handler = _node_class_attr__CompoundType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__CompoundType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__CompoundType,
            _node_class_finish__CompoundType,
            None,
        )
    )


def _node_class_finish_fields__CompoundType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "name"):
        _raise_missing_element_error(state, obj, "name")


def _node_class_finish__CompoundType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__CompoundType(state, n)


_node_class_attr__DoxygenType = _cur_list = {}


@_add_to_list("version")
def _a__DoxygenType__version(state: _ParseState, obj, value: str, /):
    obj.version = value


def _node_class_attr_end__DoxygenType(state: _ParseState, obj, /):
    if not hasattr(obj, "version"):
        _raise_missing_attribute_error(state, "version")


_node_class_child__DoxygenType = _cur_list = {}


@_add_to_list("compounddef")
def _e__DoxygenType__compounddef(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__compounddefType(state, obj.compounddef.append, attr)


def _node_class_start__DoxygenType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_DoxygenType.__new__(Node_DoxygenType)

    n.compounddef = []
    for name, value in attr:
        handler = _node_class_attr__DoxygenType.get(name)

        if handler is not None:
            handler(state, n, value)
    _node_class_attr_end__DoxygenType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__DoxygenType,
            None,
            None,
        )
    )


_node_class_attr__compounddefType = _cur_list = {}


@_add_to_list("abstract")
def _a__compounddefType__abstract(state: _ParseState, obj, value: str, /):
    obj.abstract = _parse_DoxBool_attribute(state, "abstract", value)


@_add_to_list("final")
def _a__compounddefType__final(state: _ParseState, obj, value: str, /):
    obj.final = _parse_DoxBool_attribute(state, "final", value)


@_add_to_list("id")
def _a__compounddefType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


@_add_to_list("inline")
def _a__compounddefType__inline(state: _ParseState, obj, value: str, /):
    obj.inline = _parse_DoxBool_attribute(state, "inline", value)


@_add_to_list("kind")
def _a__compounddefType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = DoxCompoundKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("language")
def _a__compounddefType__language(state: _ParseState, obj, value: str, /):
    try:
        obj.language = DoxLanguage(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("prot")
def _a__compounddefType__prot(state: _ParseState, obj, value: str, /):
    try:
        obj.prot = DoxProtectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("sealed")
def _a__compounddefType__sealed(state: _ParseState, obj, value: str, /):
    obj.sealed = _parse_DoxBool_attribute(state, "sealed", value)


def _node_class_attr_end__compounddefType(state: _ParseState, obj, /):
    if not hasattr(obj, "abstract"):
        obj.abstract = None
    if not hasattr(obj, "final"):
        obj.final = None
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")
    if not hasattr(obj, "inline"):
        obj.inline = None
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")
    if not hasattr(obj, "language"):
        obj.language = None
    if not hasattr(obj, "prot"):
        obj.prot = None
    if not hasattr(obj, "sealed"):
        obj.sealed = None


_node_class_child__compounddefType = _cur_list = {}


@_add_to_list("basecompoundref")
def _e__compounddefType__basecompoundref(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__compoundRefType(state, obj.basecompoundref.append, attr)


@_add_to_list("briefdescription")
def _e__compounddefType__briefdescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "briefdescription"):
        _raise_duplicate_element_error(state, "briefdescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "briefdescription"), attr
    )


@_add_to_list("collaborationgraph")
def _e__compounddefType__collaborationgraph(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "collaborationgraph"):
        _raise_duplicate_element_error(state, "collaborationgraph")

    _node_class_start__graphType(state, functools.partial(setattr, obj, "collaborationgraph"), attr)


@_add_to_list("compoundname")
def _e__compounddefType__compoundname(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "compoundname"):
        _raise_duplicate_element_error(state, "compoundname")

    _node_start_string(state, functools.partial(setattr, obj, "compoundname"), attr)


@_add_to_list("derivedcompoundref")
def _e__compounddefType__derivedcompoundref(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__compoundRefType(state, obj.derivedcompoundref.append, attr)


@_add_to_list("detaileddescription")
def _e__compounddefType__detaileddescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "detaileddescription"):
        _raise_duplicate_element_error(state, "detaileddescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "detaileddescription"), attr
    )


@_add_to_list("exports")
def _e__compounddefType__exports(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "exports"):
        _raise_duplicate_element_error(state, "exports")

    _node_class_start__exportsType(state, functools.partial(setattr, obj, "exports"), attr)


@_add_to_list("incdepgraph")
def _e__compounddefType__incdepgraph(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "incdepgraph"):
        _raise_duplicate_element_error(state, "incdepgraph")

    _node_class_start__graphType(state, functools.partial(setattr, obj, "incdepgraph"), attr)


@_add_to_list("includedby")
def _e__compounddefType__includedby(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__incType(state, obj.includedby.append, attr)


@_add_to_list("includes")
def _e__compounddefType__includes(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__incType(state, obj.includes.append, attr)


@_add_to_list("inheritancegraph")
def _e__compounddefType__inheritancegraph(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "inheritancegraph"):
        _raise_duplicate_element_error(state, "inheritancegraph")

    _node_class_start__graphType(state, functools.partial(setattr, obj, "inheritancegraph"), attr)


@_add_to_list("initializer")
def _e__compounddefType__initializer(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "initializer"):
        _raise_duplicate_element_error(state, "initializer")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "initializer"), attr)


@_add_to_list("innerclass")
def _e__compounddefType__innerclass(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innerclass.append, attr)


@_add_to_list("innerconcept")
def _e__compounddefType__innerconcept(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innerconcept.append, attr)


@_add_to_list("innerdir")
def _e__compounddefType__innerdir(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innerdir.append, attr)


@_add_to_list("innerfile")
def _e__compounddefType__innerfile(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innerfile.append, attr)


@_add_to_list("innergroup")
def _e__compounddefType__innergroup(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innergroup.append, attr)


@_add_to_list("innermodule")
def _e__compounddefType__innermodule(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innermodule.append, attr)


@_add_to_list("innernamespace")
def _e__compounddefType__innernamespace(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__refType(state, obj.innernamespace.append, attr)


@_add_to_list("innerpage")
def _e__compounddefType__innerpage(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refType(state, obj.innerpage.append, attr)


@_add_to_list("invincdepgraph")
def _e__compounddefType__invincdepgraph(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "invincdepgraph"):
        _raise_duplicate_element_error(state, "invincdepgraph")

    _node_class_start__graphType(state, functools.partial(setattr, obj, "invincdepgraph"), attr)


@_add_to_list("listofallmembers")
def _e__compounddefType__listofallmembers(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "listofallmembers"):
        _raise_duplicate_element_error(state, "listofallmembers")

    _node_class_start__listofallmembersType(
        state, functools.partial(setattr, obj, "listofallmembers"), attr
    )


@_add_to_list("location")
def _e__compounddefType__location(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "location"):
        _raise_duplicate_element_error(state, "location")

    _node_class_start__locationType(state, functools.partial(setattr, obj, "location"), attr)


@_add_to_list("programlisting")
def _e__compounddefType__programlisting(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "programlisting"):
        _raise_duplicate_element_error(state, "programlisting")

    _node_class_start__listingType(state, functools.partial(setattr, obj, "programlisting"), attr)


@_add_to_list("qualifier")
def _e__compounddefType__qualifier(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, obj.qualifier.append, attr)


@_add_to_list("requiresclause")
def _e__compounddefType__requiresclause(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "requiresclause"):
        _raise_duplicate_element_error(state, "requiresclause")

    _node_class_start__linkedTextType(
        state, functools.partial(setattr, obj, "requiresclause"), attr
    )


@_add_to_list("sectiondef")
def _e__compounddefType__sectiondef(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__sectiondefType(state, obj.sectiondef.append, attr)


@_add_to_list("tableofcontents")
def _e__compounddefType__tableofcontents(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "tableofcontents"):
        _raise_duplicate_element_error(state, "tableofcontents")

    _node_class_start__tableofcontentsType(
        state, functools.partial(setattr, obj, "tableofcontents"), attr
    )


@_add_to_list("templateparamlist")
def _e__compounddefType__templateparamlist(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "templateparamlist"):
        _raise_duplicate_element_error(state, "templateparamlist")

    _node_class_start__templateparamlistType(
        state, functools.partial(setattr, obj, "templateparamlist"), attr
    )


@_add_to_list("title")
def _e__compounddefType__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_start_string(state, functools.partial(setattr, obj, "title"), attr)


def _node_class_start__compounddefType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_compounddefType.__new__(Node_compounddefType)

    n.basecompoundref = []
    n.derivedcompoundref = []
    n.includedby = []
    n.includes = []
    n.innerclass = []
    n.innerconcept = []
    n.innerdir = []
    n.innerfile = []
    n.innergroup = []
    n.innermodule = []
    n.innernamespace = []
    n.innerpage = []
    n.qualifier = []
    n.sectiondef = []
    for name, value in attr:
        handler = _node_class_attr__compounddefType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__compounddefType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__compounddefType,
            _node_class_finish__compounddefType,
            None,
        )
    )


def _node_class_finish_fields__compounddefType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "briefdescription"):
        obj.briefdescription = None
    if not hasattr(obj, "collaborationgraph"):
        obj.collaborationgraph = None
    if not hasattr(obj, "compoundname"):
        _raise_missing_element_error(state, obj, "compoundname")
    if not hasattr(obj, "detaileddescription"):
        obj.detaileddescription = None
    if not hasattr(obj, "exports"):
        obj.exports = None
    if not hasattr(obj, "incdepgraph"):
        obj.incdepgraph = None
    if not hasattr(obj, "inheritancegraph"):
        obj.inheritancegraph = None
    if not hasattr(obj, "initializer"):
        obj.initializer = None
    if not hasattr(obj, "invincdepgraph"):
        obj.invincdepgraph = None
    if not hasattr(obj, "listofallmembers"):
        obj.listofallmembers = None
    if not hasattr(obj, "location"):
        obj.location = None
    if not hasattr(obj, "programlisting"):
        obj.programlisting = None
    if not hasattr(obj, "requiresclause"):
        obj.requiresclause = None
    if not hasattr(obj, "tableofcontents"):
        obj.tableofcontents = None
    if not hasattr(obj, "templateparamlist"):
        obj.templateparamlist = None
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__compounddefType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__compounddefType(state, n)


_node_class_child__graphType = _cur_list = {}


@_add_to_list("node")
def _e__graphType__node(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__nodeType(state, obj.node.append, attr)


def _node_class_start__graphType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_graphType.__new__(Node_graphType)

    n.node = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__graphType,
            _node_class_finish__graphType,
            None,
        )
    )


def _node_class_finish_fields__graphType(state: _ParseState, obj, /) -> None:
    if len(obj.node) < 1:
        _raise_empty_list_element_error(state, "node")


def _node_class_finish__graphType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__graphType(state, n)


_node_class_child__templateparamlistType = _cur_list = {}


@_add_to_list("param")
def _e__templateparamlistType__param(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__paramType(state, obj.param.append, attr)


def _node_class_start__templateparamlistType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_templateparamlistType.__new__(Node_templateparamlistType)

    n.param = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__templateparamlistType,
            None,
            None,
        )
    )


_node_class_attr__sectiondefType = _cur_list = {}


@_add_to_list("kind")
def _a__sectiondefType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = DoxSectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__sectiondefType(state: _ParseState, obj, /):
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")


_node_class_child__sectiondefType = _cur_list = {}


@_add_to_list("description")
def _e__sectiondefType__description(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "description"):
        _raise_duplicate_element_error(state, "description")

    _node_class_start__descriptionType(state, functools.partial(setattr, obj, "description"), attr)


@_add_to_list("header")
def _e__sectiondefType__header(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "header"):
        _raise_duplicate_element_error(state, "header")

    _node_start_string(state, functools.partial(setattr, obj, "header"), attr)


@_add_to_list("member")
def _e__sectiondefType__member(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__MemberType(state, obj.member.append, attr)


@_add_to_list("memberdef")
def _e__sectiondefType__memberdef(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__memberdefType(state, obj.memberdef.append, attr)


def _node_class_start__sectiondefType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_sectiondefType.__new__(Node_sectiondefType)

    n.member = []
    n.memberdef = []
    for name, value in attr:
        handler = _node_class_attr__sectiondefType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__sectiondefType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__sectiondefType,
            _node_class_finish__sectiondefType,
            None,
        )
    )


def _node_class_finish_fields__sectiondefType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "description"):
        obj.description = None
    if not hasattr(obj, "header"):
        obj.header = None


def _node_class_finish__sectiondefType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__sectiondefType(state, n)


_node_class_child__tableofcontentsType = _cur_list = {}


@_add_to_list("tocsect")
def _e__tableofcontentsType__tocsect(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__tableofcontentsKindType(state, obj.tocsect.append, attr)


@_add_to_list("tableofcontents")
def _e__tableofcontentsType__tableofcontents(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__tableofcontentsType(state, obj.tableofcontents.append, attr)


def _node_class_start__tableofcontentsType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_tableofcontentsType.__new__(Node_tableofcontentsType)

    n.tocsect = []
    n.tableofcontents = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__tableofcontentsType,
            _node_class_finish__tableofcontentsType,
            None,
        )
    )


def _node_class_finish_fields__tableofcontentsType(state: _ParseState, obj, /) -> None:
    if len(obj.tocsect) < 1:
        _raise_empty_list_element_error(state, "tocsect")


def _node_class_finish__tableofcontentsType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__tableofcontentsType(state, n)


_node_class_child__linkedTextType = _cur_list = {}


@_add_to_list("ref")
def _e__linkedTextType__ref(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refTextType(state, (lambda x: obj.append(TaggedValue("ref", x))), attr)


def _node_class_start__linkedTextType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_linkedTextType.__new__(Node_linkedTextType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__linkedTextType,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__descriptionType = _cur_list = {}


@_add_to_list("title")
def _e__descriptionType__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_start_string(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("internal")
def _e__descriptionType__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalType(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


@_add_to_list("para")
def _e__descriptionType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect1")
def _e__descriptionType__sect1(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect1Type(state, (lambda x: obj.append(TaggedValue("sect1", x))), attr)


@_add_to_list("sect2")
def _e__descriptionType__sect2(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect2Type(state, (lambda x: obj.append(TaggedValue("sect2", x))), attr)


@_add_to_list("sect3")
def _e__descriptionType__sect3(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect3Type(state, (lambda x: obj.append(TaggedValue("sect3", x))), attr)


def _node_class_start__descriptionType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_descriptionType.__new__(Node_descriptionType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__descriptionType,
            _node_class_finish__descriptionType,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__descriptionType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__descriptionType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__descriptionType(state, n)


_node_class_child__exportsType = _cur_list = {}


@_add_to_list("export")
def _e__exportsType__export(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__exportType(state, obj.export.append, attr)


def _node_class_start__exportsType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_exportsType.__new__(Node_exportsType)

    n.export = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__exportsType,
            _node_class_finish__exportsType,
            None,
        )
    )


def _node_class_finish_fields__exportsType(state: _ParseState, obj, /) -> None:
    if len(obj.export) < 1:
        _raise_empty_list_element_error(state, "export")


def _node_class_finish__exportsType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__exportsType(state, n)


_node_class_attr__listingType = _cur_list = {}


@_add_to_list("filename")
def _a__listingType__filename(state: _ParseState, obj, value: str, /):
    obj.filename = value


def _node_class_attr_end__listingType(state: _ParseState, obj, /):
    if not hasattr(obj, "filename"):
        obj.filename = None


_node_class_child__listingType = _cur_list = {}


@_add_to_list("codeline")
def _e__listingType__codeline(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__codelineType(state, obj.codeline.append, attr)


def _node_class_start__listingType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_listingType.__new__(Node_listingType)

    n.codeline = []
    for name, value in attr:
        handler = _node_class_attr__listingType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__listingType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__listingType,
            None,
            None,
        )
    )


_node_class_attr__locationType = _cur_list = {}


@_add_to_list("bodyend")
def _a__locationType__bodyend(state: _ParseState, obj, value: str, /):
    try:
        obj.bodyend = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("bodyfile")
def _a__locationType__bodyfile(state: _ParseState, obj, value: str, /):
    obj.bodyfile = value


@_add_to_list("bodystart")
def _a__locationType__bodystart(state: _ParseState, obj, value: str, /):
    try:
        obj.bodystart = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("column")
def _a__locationType__column(state: _ParseState, obj, value: str, /):
    try:
        obj.column = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("declcolumn")
def _a__locationType__declcolumn(state: _ParseState, obj, value: str, /):
    try:
        obj.declcolumn = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("declfile")
def _a__locationType__declfile(state: _ParseState, obj, value: str, /):
    obj.declfile = value


@_add_to_list("declline")
def _a__locationType__declline(state: _ParseState, obj, value: str, /):
    try:
        obj.declline = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("file")
def _a__locationType__file(state: _ParseState, obj, value: str, /):
    obj.file = value


@_add_to_list("line")
def _a__locationType__line(state: _ParseState, obj, value: str, /):
    try:
        obj.line = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


def _node_class_attr_end__locationType(state: _ParseState, obj, /):
    if not hasattr(obj, "bodyend"):
        obj.bodyend = None
    if not hasattr(obj, "bodyfile"):
        obj.bodyfile = None
    if not hasattr(obj, "bodystart"):
        obj.bodystart = None
    if not hasattr(obj, "column"):
        obj.column = None
    if not hasattr(obj, "declcolumn"):
        obj.declcolumn = None
    if not hasattr(obj, "declfile"):
        obj.declfile = None
    if not hasattr(obj, "declline"):
        obj.declline = None
    if not hasattr(obj, "file"):
        _raise_missing_attribute_error(state, "file")
    if not hasattr(obj, "line"):
        obj.line = None


def _node_class_start__locationType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_locationType.__new__(Node_locationType)

    for name, value in attr:
        handler = _node_class_attr__locationType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__locationType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            None,
        )
    )


_node_class_child__listofallmembersType = _cur_list = {}


@_add_to_list("member")
def _e__listofallmembersType__member(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__memberRefType(state, obj.member.append, attr)


def _node_class_start__listofallmembersType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_listofallmembersType.__new__(Node_listofallmembersType)

    n.member = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__listofallmembersType,
            None,
            None,
        )
    )


_node_class_attr__memberRefType = _cur_list = {}


@_add_to_list("ambiguityscope")
def _a__memberRefType__ambiguityscope(state: _ParseState, obj, value: str, /):
    obj.ambiguityscope = value


@_add_to_list("prot")
def _a__memberRefType__prot(state: _ParseState, obj, value: str, /):
    try:
        obj.prot = DoxProtectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("refid")
def _a__memberRefType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("virt")
def _a__memberRefType__virt(state: _ParseState, obj, value: str, /):
    try:
        obj.virt = DoxVirtualKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__memberRefType(state: _ParseState, obj, /):
    if not hasattr(obj, "ambiguityscope"):
        obj.ambiguityscope = None
    if not hasattr(obj, "prot"):
        _raise_missing_attribute_error(state, "prot")
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "virt"):
        _raise_missing_attribute_error(state, "virt")


_node_class_child__memberRefType = _cur_list = {}


@_add_to_list("name")
def _e__memberRefType__name(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "name"):
        _raise_duplicate_element_error(state, "name")

    _node_start_string(state, functools.partial(setattr, obj, "name"), attr)


@_add_to_list("scope")
def _e__memberRefType__scope(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "scope"):
        _raise_duplicate_element_error(state, "scope")

    _node_start_string(state, functools.partial(setattr, obj, "scope"), attr)


def _node_class_start__memberRefType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_memberRefType.__new__(Node_memberRefType)

    for name, value in attr:
        handler = _node_class_attr__memberRefType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__memberRefType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__memberRefType,
            _node_class_finish__memberRefType,
            None,
        )
    )


def _node_class_finish_fields__memberRefType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "name"):
        _raise_missing_element_error(state, obj, "name")
    if not hasattr(obj, "scope"):
        _raise_missing_element_error(state, obj, "scope")


def _node_class_finish__memberRefType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__memberRefType(state, n)


_node_class_attr__memberdefType = _cur_list = {}


@_add_to_list("accessor")
def _a__memberdefType__accessor(state: _ParseState, obj, value: str, /):
    try:
        obj.accessor = DoxAccessor(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("add")
def _a__memberdefType__add(state: _ParseState, obj, value: str, /):
    obj.add = _parse_DoxBool_attribute(state, "add", value)


@_add_to_list("attribute")
def _a__memberdefType__attribute(state: _ParseState, obj, value: str, /):
    obj.attribute = _parse_DoxBool_attribute(state, "attribute", value)


@_add_to_list("bound")
def _a__memberdefType__bound(state: _ParseState, obj, value: str, /):
    obj.bound = _parse_DoxBool_attribute(state, "bound", value)


@_add_to_list("const")
def _a__memberdefType__const(state: _ParseState, obj, value: str, /):
    obj.const = _parse_DoxBool_attribute(state, "const", value)


@_add_to_list("constexpr")
def _a__memberdefType__constexpr(state: _ParseState, obj, value: str, /):
    obj.constexpr = _parse_DoxBool_attribute(state, "constexpr", value)


@_add_to_list("consteval")
def _a__memberdefType__consteval(state: _ParseState, obj, value: str, /):
    obj.consteval = _parse_DoxBool_attribute(state, "consteval", value)


@_add_to_list("constinit")
def _a__memberdefType__constinit(state: _ParseState, obj, value: str, /):
    obj.constinit = _parse_DoxBool_attribute(state, "constinit", value)


@_add_to_list("constrained")
def _a__memberdefType__constrained(state: _ParseState, obj, value: str, /):
    obj.constrained = _parse_DoxBool_attribute(state, "constrained", value)


@_add_to_list("explicit")
def _a__memberdefType__explicit(state: _ParseState, obj, value: str, /):
    obj.explicit = _parse_DoxBool_attribute(state, "explicit", value)


@_add_to_list("extern")
def _a__memberdefType__extern(state: _ParseState, obj, value: str, /):
    obj.extern = _parse_DoxBool_attribute(state, "extern", value)


@_add_to_list("final")
def _a__memberdefType__final(state: _ParseState, obj, value: str, /):
    obj.final = _parse_DoxBool_attribute(state, "final", value)


@_add_to_list("gettable")
def _a__memberdefType__gettable(state: _ParseState, obj, value: str, /):
    obj.gettable = _parse_DoxBool_attribute(state, "gettable", value)


@_add_to_list("id")
def _a__memberdefType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


@_add_to_list("initonly")
def _a__memberdefType__initonly(state: _ParseState, obj, value: str, /):
    obj.initonly = _parse_DoxBool_attribute(state, "initonly", value)


@_add_to_list("inline")
def _a__memberdefType__inline(state: _ParseState, obj, value: str, /):
    obj.inline = _parse_DoxBool_attribute(state, "inline", value)


@_add_to_list("kind")
def _a__memberdefType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = DoxMemberKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("maybeambiguous")
def _a__memberdefType__maybeambiguous(state: _ParseState, obj, value: str, /):
    obj.maybeambiguous = _parse_DoxBool_attribute(state, "maybeambiguous", value)


@_add_to_list("maybedefault")
def _a__memberdefType__maybedefault(state: _ParseState, obj, value: str, /):
    obj.maybedefault = _parse_DoxBool_attribute(state, "maybedefault", value)


@_add_to_list("maybevoid")
def _a__memberdefType__maybevoid(state: _ParseState, obj, value: str, /):
    obj.maybevoid = _parse_DoxBool_attribute(state, "maybevoid", value)


@_add_to_list("mutable")
def _a__memberdefType__mutable(state: _ParseState, obj, value: str, /):
    obj.mutable = _parse_DoxBool_attribute(state, "mutable", value)


@_add_to_list("new")
def _a__memberdefType__new(state: _ParseState, obj, value: str, /):
    obj.new = _parse_DoxBool_attribute(state, "new", value)


@_add_to_list("nodiscard")
def _a__memberdefType__nodiscard(state: _ParseState, obj, value: str, /):
    obj.nodiscard = _parse_DoxBool_attribute(state, "nodiscard", value)


@_add_to_list("noexcept")
def _a__memberdefType__noexcept(state: _ParseState, obj, value: str, /):
    obj.noexcept = _parse_DoxBool_attribute(state, "noexcept", value)


@_add_to_list("noexceptexpression")
def _a__memberdefType__noexceptexpression(state: _ParseState, obj, value: str, /):
    obj.noexceptexpression = value


@_add_to_list("optional")
def _a__memberdefType__optional(state: _ParseState, obj, value: str, /):
    obj.optional = _parse_DoxBool_attribute(state, "optional", value)


@_add_to_list("privategettable")
def _a__memberdefType__privategettable(state: _ParseState, obj, value: str, /):
    obj.privategettable = _parse_DoxBool_attribute(state, "privategettable", value)


@_add_to_list("privatesettable")
def _a__memberdefType__privatesettable(state: _ParseState, obj, value: str, /):
    obj.privatesettable = _parse_DoxBool_attribute(state, "privatesettable", value)


@_add_to_list("property")
def _a__memberdefType__property(state: _ParseState, obj, value: str, /):
    obj.property = _parse_DoxBool_attribute(state, "property", value)


@_add_to_list("prot")
def _a__memberdefType__prot(state: _ParseState, obj, value: str, /):
    try:
        obj.prot = DoxProtectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("protectedgettable")
def _a__memberdefType__protectedgettable(state: _ParseState, obj, value: str, /):
    obj.protectedgettable = _parse_DoxBool_attribute(state, "protectedgettable", value)


@_add_to_list("protectedsettable")
def _a__memberdefType__protectedsettable(state: _ParseState, obj, value: str, /):
    obj.protectedsettable = _parse_DoxBool_attribute(state, "protectedsettable", value)


@_add_to_list("raise")
def _a__memberdefType__raise(state: _ParseState, obj, value: str, /):
    obj.raise_ = _parse_DoxBool_attribute(state, "raise", value)


@_add_to_list("readable")
def _a__memberdefType__readable(state: _ParseState, obj, value: str, /):
    obj.readable = _parse_DoxBool_attribute(state, "readable", value)


@_add_to_list("readonly")
def _a__memberdefType__readonly(state: _ParseState, obj, value: str, /):
    obj.readonly = _parse_DoxBool_attribute(state, "readonly", value)


@_add_to_list("refqual")
def _a__memberdefType__refqual(state: _ParseState, obj, value: str, /):
    try:
        obj.refqual = DoxRefQualifierKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("removable")
def _a__memberdefType__removable(state: _ParseState, obj, value: str, /):
    obj.removable = _parse_DoxBool_attribute(state, "removable", value)


@_add_to_list("remove")
def _a__memberdefType__remove(state: _ParseState, obj, value: str, /):
    obj.remove = _parse_DoxBool_attribute(state, "remove", value)


@_add_to_list("required")
def _a__memberdefType__required(state: _ParseState, obj, value: str, /):
    obj.required = _parse_DoxBool_attribute(state, "required", value)


@_add_to_list("sealed")
def _a__memberdefType__sealed(state: _ParseState, obj, value: str, /):
    obj.sealed = _parse_DoxBool_attribute(state, "sealed", value)


@_add_to_list("settable")
def _a__memberdefType__settable(state: _ParseState, obj, value: str, /):
    obj.settable = _parse_DoxBool_attribute(state, "settable", value)


@_add_to_list("static")
def _a__memberdefType__static(state: _ParseState, obj, value: str, /):
    obj.static = _parse_DoxBool_attribute(state, "static", value)


@_add_to_list("strong")
def _a__memberdefType__strong(state: _ParseState, obj, value: str, /):
    obj.strong = _parse_DoxBool_attribute(state, "strong", value)


@_add_to_list("transient")
def _a__memberdefType__transient(state: _ParseState, obj, value: str, /):
    obj.transient = _parse_DoxBool_attribute(state, "transient", value)


@_add_to_list("virt")
def _a__memberdefType__virt(state: _ParseState, obj, value: str, /):
    try:
        obj.virt = DoxVirtualKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("volatile")
def _a__memberdefType__volatile(state: _ParseState, obj, value: str, /):
    obj.volatile = _parse_DoxBool_attribute(state, "volatile", value)


@_add_to_list("writable")
def _a__memberdefType__writable(state: _ParseState, obj, value: str, /):
    obj.writable = _parse_DoxBool_attribute(state, "writable", value)


def _node_class_attr_end__memberdefType(state: _ParseState, obj, /):
    if not hasattr(obj, "accessor"):
        obj.accessor = None
    if not hasattr(obj, "add"):
        obj.add = None
    if not hasattr(obj, "attribute"):
        obj.attribute = None
    if not hasattr(obj, "bound"):
        obj.bound = None
    if not hasattr(obj, "const"):
        obj.const = None
    if not hasattr(obj, "constexpr"):
        obj.constexpr = None
    if not hasattr(obj, "consteval"):
        obj.consteval = None
    if not hasattr(obj, "constinit"):
        obj.constinit = None
    if not hasattr(obj, "constrained"):
        obj.constrained = None
    if not hasattr(obj, "explicit"):
        obj.explicit = None
    if not hasattr(obj, "extern"):
        obj.extern = None
    if not hasattr(obj, "final"):
        obj.final = None
    if not hasattr(obj, "gettable"):
        obj.gettable = None
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")
    if not hasattr(obj, "initonly"):
        obj.initonly = None
    if not hasattr(obj, "inline"):
        obj.inline = None
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")
    if not hasattr(obj, "maybeambiguous"):
        obj.maybeambiguous = None
    if not hasattr(obj, "maybedefault"):
        obj.maybedefault = None
    if not hasattr(obj, "maybevoid"):
        obj.maybevoid = None
    if not hasattr(obj, "mutable"):
        obj.mutable = None
    if not hasattr(obj, "new"):
        obj.new = None
    if not hasattr(obj, "nodiscard"):
        obj.nodiscard = None
    if not hasattr(obj, "noexcept"):
        obj.noexcept = None
    if not hasattr(obj, "noexceptexpression"):
        obj.noexceptexpression = None
    if not hasattr(obj, "optional"):
        obj.optional = None
    if not hasattr(obj, "privategettable"):
        obj.privategettable = None
    if not hasattr(obj, "privatesettable"):
        obj.privatesettable = None
    if not hasattr(obj, "property"):
        obj.property = None
    if not hasattr(obj, "prot"):
        _raise_missing_attribute_error(state, "prot")
    if not hasattr(obj, "protectedgettable"):
        obj.protectedgettable = None
    if not hasattr(obj, "protectedsettable"):
        obj.protectedsettable = None
    if not hasattr(obj, "raise_"):
        obj.raise_ = None
    if not hasattr(obj, "readable"):
        obj.readable = None
    if not hasattr(obj, "readonly"):
        obj.readonly = None
    if not hasattr(obj, "refqual"):
        obj.refqual = None
    if not hasattr(obj, "removable"):
        obj.removable = None
    if not hasattr(obj, "remove"):
        obj.remove = None
    if not hasattr(obj, "required"):
        obj.required = None
    if not hasattr(obj, "sealed"):
        obj.sealed = None
    if not hasattr(obj, "settable"):
        obj.settable = None
    if not hasattr(obj, "static"):
        _raise_missing_attribute_error(state, "static")
    if not hasattr(obj, "strong"):
        obj.strong = None
    if not hasattr(obj, "transient"):
        obj.transient = None
    if not hasattr(obj, "virt"):
        obj.virt = None
    if not hasattr(obj, "volatile"):
        obj.volatile = None
    if not hasattr(obj, "writable"):
        obj.writable = None


_node_class_child__memberdefType = _cur_list = {}


@_add_to_list("argsstring")
def _e__memberdefType__argsstring(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "argsstring"):
        _raise_duplicate_element_error(state, "argsstring")

    _node_start_string(state, functools.partial(setattr, obj, "argsstring"), attr)


@_add_to_list("bitfield")
def _e__memberdefType__bitfield(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "bitfield"):
        _raise_duplicate_element_error(state, "bitfield")

    _node_start_string(state, functools.partial(setattr, obj, "bitfield"), attr)


@_add_to_list("briefdescription")
def _e__memberdefType__briefdescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "briefdescription"):
        _raise_duplicate_element_error(state, "briefdescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "briefdescription"), attr
    )


@_add_to_list("definition")
def _e__memberdefType__definition(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "definition"):
        _raise_duplicate_element_error(state, "definition")

    _node_start_string(state, functools.partial(setattr, obj, "definition"), attr)


@_add_to_list("detaileddescription")
def _e__memberdefType__detaileddescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "detaileddescription"):
        _raise_duplicate_element_error(state, "detaileddescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "detaileddescription"), attr
    )


@_add_to_list("enumvalue")
def _e__memberdefType__enumvalue(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__enumvalueType(state, obj.enumvalue.append, attr)


@_add_to_list("exceptions")
def _e__memberdefType__exceptions(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "exceptions"):
        _raise_duplicate_element_error(state, "exceptions")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "exceptions"), attr)


@_add_to_list("inbodydescription")
def _e__memberdefType__inbodydescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "inbodydescription"):
        _raise_duplicate_element_error(state, "inbodydescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "inbodydescription"), attr
    )


@_add_to_list("initializer")
def _e__memberdefType__initializer(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "initializer"):
        _raise_duplicate_element_error(state, "initializer")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "initializer"), attr)


@_add_to_list("location")
def _e__memberdefType__location(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "location"):
        _raise_duplicate_element_error(state, "location")

    _node_class_start__locationType(state, functools.partial(setattr, obj, "location"), attr)


@_add_to_list("name")
def _e__memberdefType__name(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "name"):
        _raise_duplicate_element_error(state, "name")

    _node_start_string(state, functools.partial(setattr, obj, "name"), attr)


@_add_to_list("param")
def _e__memberdefType__param(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__paramType(state, obj.param.append, attr)


@_add_to_list("qualifiedname")
def _e__memberdefType__qualifiedname(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "qualifiedname"):
        _raise_duplicate_element_error(state, "qualifiedname")

    _node_start_string(state, functools.partial(setattr, obj, "qualifiedname"), attr)


@_add_to_list("qualifier")
def _e__memberdefType__qualifier(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, obj.qualifier.append, attr)


@_add_to_list("read")
def _e__memberdefType__read(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "read"):
        _raise_duplicate_element_error(state, "read")

    _node_start_string(state, functools.partial(setattr, obj, "read"), attr)


@_add_to_list("referencedby")
def _e__memberdefType__referencedby(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__referenceType(state, obj.referencedby.append, attr)


@_add_to_list("references")
def _e__memberdefType__references(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__referenceType(state, obj.references.append, attr)


@_add_to_list("reimplementedby")
def _e__memberdefType__reimplementedby(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__reimplementType(state, obj.reimplementedby.append, attr)


@_add_to_list("reimplements")
def _e__memberdefType__reimplements(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__reimplementType(state, obj.reimplements.append, attr)


@_add_to_list("requiresclause")
def _e__memberdefType__requiresclause(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "requiresclause"):
        _raise_duplicate_element_error(state, "requiresclause")

    _node_class_start__linkedTextType(
        state, functools.partial(setattr, obj, "requiresclause"), attr
    )


@_add_to_list("templateparamlist")
def _e__memberdefType__templateparamlist(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "templateparamlist"):
        _raise_duplicate_element_error(state, "templateparamlist")

    _node_class_start__templateparamlistType(
        state, functools.partial(setattr, obj, "templateparamlist"), attr
    )


@_add_to_list("type")
def _e__memberdefType__type(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "type"):
        _raise_duplicate_element_error(state, "type")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "type"), attr)


@_add_to_list("write")
def _e__memberdefType__write(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "write"):
        _raise_duplicate_element_error(state, "write")

    _node_start_string(state, functools.partial(setattr, obj, "write"), attr)


def _node_class_start__memberdefType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_memberdefType.__new__(Node_memberdefType)

    n.enumvalue = []
    n.param = []
    n.qualifier = []
    n.referencedby = []
    n.references = []
    n.reimplementedby = []
    n.reimplements = []
    for name, value in attr:
        handler = _node_class_attr__memberdefType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__memberdefType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__memberdefType,
            _node_class_finish__memberdefType,
            None,
        )
    )


def _node_class_finish_fields__memberdefType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "argsstring"):
        obj.argsstring = None
    if not hasattr(obj, "bitfield"):
        obj.bitfield = None
    if not hasattr(obj, "briefdescription"):
        obj.briefdescription = None
    if not hasattr(obj, "definition"):
        obj.definition = None
    if not hasattr(obj, "detaileddescription"):
        obj.detaileddescription = None
    if not hasattr(obj, "exceptions"):
        obj.exceptions = None
    if not hasattr(obj, "inbodydescription"):
        obj.inbodydescription = None
    if not hasattr(obj, "initializer"):
        obj.initializer = None
    if not hasattr(obj, "location"):
        _raise_missing_element_error(state, obj, "location")
    if not hasattr(obj, "name"):
        _raise_missing_element_error(state, obj, "name")
    if not hasattr(obj, "qualifiedname"):
        obj.qualifiedname = None
    if not hasattr(obj, "read"):
        obj.read = None
    if not hasattr(obj, "requiresclause"):
        obj.requiresclause = None
    if not hasattr(obj, "templateparamlist"):
        obj.templateparamlist = None
    if not hasattr(obj, "type"):
        obj.type = None
    if not hasattr(obj, "write"):
        obj.write = None


def _node_class_finish__memberdefType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__memberdefType(state, n)


_node_class_attr__MemberType = _cur_list = {}


@_add_to_list("kind")
def _a__MemberType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = MemberKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("refid")
def _a__MemberType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


def _node_class_attr_end__MemberType(state: _ParseState, obj, /):
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")


_node_class_child__MemberType = _cur_list = {}


@_add_to_list("name")
def _e__MemberType__name(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "name"):
        _raise_duplicate_element_error(state, "name")

    _node_start_string(state, functools.partial(setattr, obj, "name"), attr)


def _node_class_start__MemberType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_MemberType.__new__(Node_MemberType)

    for name, value in attr:
        handler = _node_class_attr__MemberType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__MemberType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__MemberType,
            _node_class_finish__MemberType,
            None,
        )
    )


def _node_class_finish_fields__MemberType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "name"):
        _raise_missing_element_error(state, obj, "name")


def _node_class_finish__MemberType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__MemberType(state, n)


_node_class_child__paramType = _cur_list = {}


@_add_to_list("array")
def _e__paramType__array(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "array"):
        _raise_duplicate_element_error(state, "array")

    _node_start_string(state, functools.partial(setattr, obj, "array"), attr)


@_add_to_list("attributes")
def _e__paramType__attributes(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "attributes"):
        _raise_duplicate_element_error(state, "attributes")

    _node_start_string(state, functools.partial(setattr, obj, "attributes"), attr)


@_add_to_list("briefdescription")
def _e__paramType__briefdescription(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "briefdescription"):
        _raise_duplicate_element_error(state, "briefdescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "briefdescription"), attr
    )


@_add_to_list("declname")
def _e__paramType__declname(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "declname"):
        _raise_duplicate_element_error(state, "declname")

    _node_start_string(state, functools.partial(setattr, obj, "declname"), attr)


@_add_to_list("defname")
def _e__paramType__defname(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "defname"):
        _raise_duplicate_element_error(state, "defname")

    _node_start_string(state, functools.partial(setattr, obj, "defname"), attr)


@_add_to_list("defval")
def _e__paramType__defval(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "defval"):
        _raise_duplicate_element_error(state, "defval")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "defval"), attr)


@_add_to_list("type")
def _e__paramType__type(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "type"):
        _raise_duplicate_element_error(state, "type")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "type"), attr)


@_add_to_list("typeconstraint")
def _e__paramType__typeconstraint(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "typeconstraint"):
        _raise_duplicate_element_error(state, "typeconstraint")

    _node_class_start__linkedTextType(
        state, functools.partial(setattr, obj, "typeconstraint"), attr
    )


def _node_class_start__paramType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_paramType.__new__(Node_paramType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__paramType,
            _node_class_finish__paramType,
            None,
        )
    )


def _node_class_finish_fields__paramType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "array"):
        obj.array = None
    if not hasattr(obj, "attributes"):
        obj.attributes = None
    if not hasattr(obj, "briefdescription"):
        obj.briefdescription = None
    if not hasattr(obj, "declname"):
        obj.declname = None
    if not hasattr(obj, "defname"):
        obj.defname = None
    if not hasattr(obj, "defval"):
        obj.defval = None
    if not hasattr(obj, "type"):
        obj.type = None
    if not hasattr(obj, "typeconstraint"):
        obj.typeconstraint = None


def _node_class_finish__paramType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__paramType(state, n)


_node_class_attr__enumvalueType = _cur_list = {}


@_add_to_list("id")
def _a__enumvalueType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


@_add_to_list("prot")
def _a__enumvalueType__prot(state: _ParseState, obj, value: str, /):
    try:
        obj.prot = DoxProtectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__enumvalueType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")
    if not hasattr(obj, "prot"):
        _raise_missing_attribute_error(state, "prot")


_node_class_child__enumvalueType = _cur_list = {}


@_add_to_list("briefdescription")
def _e__enumvalueType__briefdescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "briefdescription"):
        _raise_duplicate_element_error(state, "briefdescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "briefdescription"), attr
    )


@_add_to_list("detaileddescription")
def _e__enumvalueType__detaileddescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "detaileddescription"):
        _raise_duplicate_element_error(state, "detaileddescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "detaileddescription"), attr
    )


@_add_to_list("initializer")
def _e__enumvalueType__initializer(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "initializer"):
        _raise_duplicate_element_error(state, "initializer")

    _node_class_start__linkedTextType(state, functools.partial(setattr, obj, "initializer"), attr)


@_add_to_list("name")
def _e__enumvalueType__name(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "name"):
        _raise_duplicate_element_error(state, "name")

    _node_start_string(state, functools.partial(setattr, obj, "name"), attr)


def _node_class_start__enumvalueType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_enumvalueType.__new__(Node_enumvalueType)

    for name, value in attr:
        handler = _node_class_attr__enumvalueType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__enumvalueType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__enumvalueType,
            _node_class_finish__enumvalueType,
            None,
        )
    )


def _node_class_finish_fields__enumvalueType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "briefdescription"):
        obj.briefdescription = None
    if not hasattr(obj, "detaileddescription"):
        obj.detaileddescription = None
    if not hasattr(obj, "initializer"):
        obj.initializer = None
    if not hasattr(obj, "name"):
        _raise_missing_element_error(state, obj, "name")


def _node_class_finish__enumvalueType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__enumvalueType(state, n)


_node_class_attr__referenceType = _cur_list = {}


@_add_to_list("compoundref")
def _a__referenceType__compoundref(state: _ParseState, obj, value: str, /):
    obj.compoundref = value


@_add_to_list("endline")
def _a__referenceType__endline(state: _ParseState, obj, value: str, /):
    try:
        obj.endline = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("refid")
def _a__referenceType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("startline")
def _a__referenceType__startline(state: _ParseState, obj, value: str, /):
    try:
        obj.startline = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


def _node_class_attr_end__referenceType(state: _ParseState, obj, /):
    if not hasattr(obj, "compoundref"):
        obj.compoundref = None
    if not hasattr(obj, "endline"):
        obj.endline = None
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "startline"):
        obj.startline = None


def _node_class_start__referenceType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_referenceType.__new__(Node_referenceType)

    for name, value in attr:
        handler = _node_class_attr__referenceType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__referenceType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docInternalType = _cur_list = {}


@_add_to_list("para")
def _e__docInternalType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect1")
def _e__docInternalType__sect1(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect1Type(state, (lambda x: obj.append(TaggedValue("sect1", x))), attr)


def _node_class_start__docInternalType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalType.__new__(Node_docInternalType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docSect1Type = _cur_list = {}


@_add_to_list("id")
def _a__docSect1Type__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docSect1Type(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docSect1Type = _cur_list = {}


@_add_to_list("title")
def _e__docSect1Type__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSect1Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect2")
def _e__docSect1Type__sect2(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect2Type(state, (lambda x: obj.append(TaggedValue("sect2", x))), attr)


@_add_to_list("sect3")
def _e__docSect1Type__sect3(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect3Type(state, (lambda x: obj.append(TaggedValue("sect3", x))), attr)


@_add_to_list("internal")
def _e__docSect1Type__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalS1Type(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


def _node_class_start__docSect1Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSect1Type.__new__(Node_docSect1Type)

    for name, value in attr:
        handler = _node_class_attr__docSect1Type.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSect1Type(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSect1Type,
            _node_class_finish__docSect1Type,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__docSect1Type(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__docSect1Type(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSect1Type(state, n)


_node_class_attr__nodeType = _cur_list = {}


@_add_to_list("id")
def _a__nodeType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__nodeType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__nodeType = _cur_list = {}


@_add_to_list("childnode")
def _e__nodeType__childnode(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__childnodeType(state, obj.childnode.append, attr)


@_add_to_list("label")
def _e__nodeType__label(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "label"):
        _raise_duplicate_element_error(state, "label")

    _node_start_string(state, functools.partial(setattr, obj, "label"), attr)


@_add_to_list("link")
def _e__nodeType__link(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "link"):
        _raise_duplicate_element_error(state, "link")

    _node_class_start__linkType(state, functools.partial(setattr, obj, "link"), attr)


def _node_class_start__nodeType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_nodeType.__new__(Node_nodeType)

    n.childnode = []
    for name, value in attr:
        handler = _node_class_attr__nodeType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__nodeType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__nodeType,
            _node_class_finish__nodeType,
            None,
        )
    )


def _node_class_finish_fields__nodeType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "label"):
        _raise_missing_element_error(state, obj, "label")
    if not hasattr(obj, "link"):
        obj.link = None


def _node_class_finish__nodeType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__nodeType(state, n)


_node_class_attr__linkType = _cur_list = {}


@_add_to_list("external")
def _a__linkType__external(state: _ParseState, obj, value: str, /):
    obj.external = value


@_add_to_list("refid")
def _a__linkType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


def _node_class_attr_end__linkType(state: _ParseState, obj, /):
    if not hasattr(obj, "external"):
        obj.external = None
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")


def _node_class_start__linkType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_linkType.__new__(Node_linkType)

    for name, value in attr:
        handler = _node_class_attr__linkType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__linkType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            None,
        )
    )


_node_class_attr__childnodeType = _cur_list = {}


@_add_to_list("refid")
def _a__childnodeType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("relation")
def _a__childnodeType__relation(state: _ParseState, obj, value: str, /):
    try:
        obj.relation = DoxGraphRelation(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__childnodeType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "relation"):
        _raise_missing_attribute_error(state, "relation")


_node_class_child__childnodeType = _cur_list = {}


@_add_to_list("edgelabel")
def _e__childnodeType__edgelabel(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, obj.edgelabel.append, attr)


def _node_class_start__childnodeType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_childnodeType.__new__(Node_childnodeType)

    n.edgelabel = []
    for name, value in attr:
        handler = _node_class_attr__childnodeType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__childnodeType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__childnodeType,
            None,
            None,
        )
    )


_node_class_attr__codelineType = _cur_list = {}


@_add_to_list("external")
def _a__codelineType__external(state: _ParseState, obj, value: str, /):
    obj.external = _parse_DoxBool_attribute(state, "external", value)


@_add_to_list("lineno")
def _a__codelineType__lineno(state: _ParseState, obj, value: str, /):
    try:
        obj.lineno = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("refid")
def _a__codelineType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("refkind")
def _a__codelineType__refkind(state: _ParseState, obj, value: str, /):
    try:
        obj.refkind = DoxRefKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__codelineType(state: _ParseState, obj, /):
    if not hasattr(obj, "external"):
        obj.external = None
    if not hasattr(obj, "lineno"):
        obj.lineno = None
    if not hasattr(obj, "refid"):
        obj.refid = None
    if not hasattr(obj, "refkind"):
        obj.refkind = None


_node_class_child__codelineType = _cur_list = {}


@_add_to_list("highlight")
def _e__codelineType__highlight(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__highlightType(state, obj.highlight.append, attr)


def _node_class_start__codelineType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_codelineType.__new__(Node_codelineType)

    n.highlight = []
    for name, value in attr:
        handler = _node_class_attr__codelineType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__codelineType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__codelineType,
            None,
            None,
        )
    )


_node_class_attr__highlightType = _cur_list = {}


@_add_to_list("class")
def _a__highlightType__class(state: _ParseState, obj, value: str, /):
    try:
        obj.class_ = DoxHighlightClass(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__highlightType(state: _ParseState, obj, /):
    if not hasattr(obj, "class_"):
        _raise_missing_attribute_error(state, "class")


_node_class_child__highlightType = _cur_list = {}


@_add_to_list("sp")
def _e__highlightType__sp(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    c = _node_start_spType(state, attr)
    if obj and type(obj[-1]) is str:
        obj[-1] += c
    else:
        obj.append(c)


@_add_to_list("ref")
def _e__highlightType__ref(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refTextType(state, (lambda x: obj.append(TaggedValue("ref", x))), attr)


def _node_class_start__highlightType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_highlightType.__new__(Node_highlightType)

    for name, value in attr:
        handler = _node_class_attr__highlightType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__highlightType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__highlightType,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docInternalS1Type = _cur_list = {}


@_add_to_list("para")
def _e__docInternalS1Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect2")
def _e__docInternalS1Type__sect2(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect2Type(state, (lambda x: obj.append(TaggedValue("sect2", x))), attr)


def _node_class_start__docInternalS1Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalS1Type.__new__(Node_docInternalS1Type)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalS1Type,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docSect2Type = _cur_list = {}


@_add_to_list("id")
def _a__docSect2Type__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docSect2Type(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docSect2Type = _cur_list = {}


@_add_to_list("title")
def _e__docSect2Type__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSect2Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect3")
def _e__docSect2Type__sect3(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect3Type(state, (lambda x: obj.append(TaggedValue("sect3", x))), attr)


@_add_to_list("internal")
def _e__docSect2Type__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalS2Type(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


def _node_class_start__docSect2Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSect2Type.__new__(Node_docSect2Type)

    for name, value in attr:
        handler = _node_class_attr__docSect2Type.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSect2Type(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSect2Type,
            _node_class_finish__docSect2Type,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__docSect2Type(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__docSect2Type(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSect2Type(state, n)


_node_class_attr__docSect3Type = _cur_list = {}


@_add_to_list("id")
def _a__docSect3Type__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docSect3Type(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docSect3Type = _cur_list = {}


@_add_to_list("title")
def _e__docSect3Type__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSect3Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect4")
def _e__docSect3Type__sect4(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect4Type(state, (lambda x: obj.append(TaggedValue("sect4", x))), attr)


@_add_to_list("internal")
def _e__docSect3Type__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalS3Type(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


def _node_class_start__docSect3Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSect3Type.__new__(Node_docSect3Type)

    for name, value in attr:
        handler = _node_class_attr__docSect3Type.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSect3Type(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSect3Type,
            _node_class_finish__docSect3Type,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__docSect3Type(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__docSect3Type(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSect3Type(state, n)


_node_class_child__docInternalS2Type = _cur_list = {}


@_add_to_list("para")
def _e__docInternalS2Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect3")
def _e__docInternalS2Type__sect3(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect3Type(state, (lambda x: obj.append(TaggedValue("sect3", x))), attr)


def _node_class_start__docInternalS2Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalS2Type.__new__(Node_docInternalS2Type)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalS2Type,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docSect4Type = _cur_list = {}


@_add_to_list("id")
def _a__docSect4Type__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docSect4Type(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docSect4Type = _cur_list = {}


@_add_to_list("title")
def _e__docSect4Type__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSect4Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect5")
def _e__docSect4Type__sect5(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect5Type(state, (lambda x: obj.append(TaggedValue("sect5", x))), attr)


@_add_to_list("internal")
def _e__docSect4Type__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalS4Type(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


def _node_class_start__docSect4Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSect4Type.__new__(Node_docSect4Type)

    for name, value in attr:
        handler = _node_class_attr__docSect4Type.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSect4Type(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSect4Type,
            _node_class_finish__docSect4Type,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__docSect4Type(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__docSect4Type(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSect4Type(state, n)


_node_class_attr__docSect5Type = _cur_list = {}


@_add_to_list("id")
def _a__docSect5Type__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docSect5Type(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docSect5Type = _cur_list = {}


@_add_to_list("title")
def _e__docSect5Type__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSect5Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect6")
def _e__docSect5Type__sect6(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect6Type(state, (lambda x: obj.append(TaggedValue("sect6", x))), attr)


@_add_to_list("internal")
def _e__docSect5Type__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalS5Type(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


def _node_class_start__docSect5Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSect5Type.__new__(Node_docSect5Type)

    for name, value in attr:
        handler = _node_class_attr__docSect5Type.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSect5Type(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSect5Type,
            _node_class_finish__docSect5Type,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__docSect5Type(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__docSect5Type(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSect5Type(state, n)


_node_class_attr__docSect6Type = _cur_list = {}


@_add_to_list("id")
def _a__docSect6Type__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docSect6Type(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docSect6Type = _cur_list = {}


@_add_to_list("title")
def _e__docSect6Type__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSect6Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("internal")
def _e__docSect6Type__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docInternalS6Type(
        state, (lambda x: obj.append(TaggedValue("internal", x))), attr
    )


def _node_class_start__docSect6Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSect6Type.__new__(Node_docSect6Type)

    for name, value in attr:
        handler = _node_class_attr__docSect6Type.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSect6Type(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSect6Type,
            _node_class_finish__docSect6Type,
            _node_list_common_text,
        )
    )


def _node_class_finish_fields__docSect6Type(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None


def _node_class_finish__docSect6Type(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSect6Type(state, n)


_node_class_child__docInternalS3Type = _cur_list = {}


@_add_to_list("para")
def _e__docInternalS3Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect3")
def _e__docInternalS3Type__sect3(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect4Type(state, (lambda x: obj.append(TaggedValue("sect3", x))), attr)


def _node_class_start__docInternalS3Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalS3Type.__new__(Node_docInternalS3Type)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalS3Type,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docInternalS4Type = _cur_list = {}


@_add_to_list("para")
def _e__docInternalS4Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect5")
def _e__docInternalS4Type__sect5(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect4Type(state, (lambda x: obj.append(TaggedValue("sect5", x))), attr)


def _node_class_start__docInternalS4Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalS4Type.__new__(Node_docInternalS4Type)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalS4Type,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docInternalS5Type = _cur_list = {}


@_add_to_list("para")
def _e__docInternalS5Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


@_add_to_list("sect6")
def _e__docInternalS5Type__sect6(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect4Type(state, (lambda x: obj.append(TaggedValue("sect6", x))), attr)


def _node_class_start__docInternalS5Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalS5Type.__new__(Node_docInternalS5Type)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalS5Type,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docInternalS6Type = _cur_list = {}


@_add_to_list("para")
def _e__docInternalS6Type__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, (lambda x: obj.append(TaggedValue("para", x))), attr)


def _node_class_start__docInternalS6Type(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docInternalS6Type.__new__(Node_docInternalS6Type)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docInternalS6Type,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docListItemType = _cur_list = {}


@_add_to_list("override")
def _a__docListItemType__override(state: _ParseState, obj, value: str, /):
    try:
        obj.override = DoxCheck(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("value")
def _a__docListItemType__value(state: _ParseState, obj, value: str, /):
    try:
        obj.value = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


def _node_class_attr_end__docListItemType(state: _ParseState, obj, /):
    if not hasattr(obj, "override"):
        obj.override = None
    if not hasattr(obj, "value"):
        obj.value = None


_node_class_child__docListItemType = _cur_list = {}


@_add_to_list("para")
def _e__docListItemType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.append, attr)


def _node_class_start__docListItemType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docListItemType.__new__(Node_docListItemType)

    for name, value in attr:
        handler = _node_class_attr__docListItemType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docListItemType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docListItemType,
            None,
            None,
        )
    )


_node_class_attr__docCaptionType = _cur_list = {}


@_add_to_list("id")
def _a__docCaptionType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docCaptionType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


def _node_class_start__docCaptionType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docCaptionType.__new__(Node_docCaptionType)

    for name, value in attr:
        handler = _node_class_attr__docCaptionType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docCaptionType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docRowType = _cur_list = {}


@_add_to_list("entry")
def _e__docRowType__entry(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docEntryType(state, obj.entry.append, attr)


def _node_class_start__docRowType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docRowType.__new__(Node_docRowType)

    n.entry = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docRowType,
            None,
            None,
        )
    )


_node_class_attr__docEntryType = _cur_list = {}


@_add_to_list("align")
def _a__docEntryType__align(state: _ParseState, obj, value: str, /):
    try:
        obj.align = DoxAlign(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("class")
def _a__docEntryType__class(state: _ParseState, obj, value: str, /):
    obj.class_ = value


@_add_to_list("colspan")
def _a__docEntryType__colspan(state: _ParseState, obj, value: str, /):
    try:
        obj.colspan = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("rowspan")
def _a__docEntryType__rowspan(state: _ParseState, obj, value: str, /):
    try:
        obj.rowspan = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("thead")
def _a__docEntryType__thead(state: _ParseState, obj, value: str, /):
    obj.thead = _parse_DoxBool_attribute(state, "thead", value)


@_add_to_list("valign")
def _a__docEntryType__valign(state: _ParseState, obj, value: str, /):
    try:
        obj.valign = DoxVerticalAlign(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("width")
def _a__docEntryType__width(state: _ParseState, obj, value: str, /):
    obj.width = value


def _node_class_attr_end__docEntryType(state: _ParseState, obj, /):
    if not hasattr(obj, "align"):
        obj.align = None
    if not hasattr(obj, "class_"):
        obj.class_ = None
    if not hasattr(obj, "colspan"):
        obj.colspan = None
    if not hasattr(obj, "rowspan"):
        obj.rowspan = None
    if not hasattr(obj, "thead"):
        _raise_missing_attribute_error(state, "thead")
    if not hasattr(obj, "valign"):
        obj.valign = None
    if not hasattr(obj, "width"):
        obj.width = None


_node_class_child__docEntryType = _cur_list = {}


@_add_to_list("para")
def _e__docEntryType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.para.append, attr)


def _node_class_start__docEntryType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docEntryType.__new__(Node_docEntryType)

    n.para = []
    for name, value in attr:
        handler = _node_class_attr__docEntryType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docEntryType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docEntryType,
            None,
            None,
        )
    )


_node_class_attr__docTocItemType = _cur_list = {}


@_add_to_list("id")
def _a__docTocItemType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docTocItemType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


def _node_class_start__docTocItemType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docTocItemType.__new__(Node_docTocItemType)

    for name, value in attr:
        handler = _node_class_attr__docTocItemType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docTocItemType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docParamListItem = _cur_list = {}


@_add_to_list("parameterdescription")
def _e__docParamListItem__parameterdescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "parameterdescription"):
        _raise_duplicate_element_error(state, "parameterdescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "parameterdescription"), attr
    )


@_add_to_list("parameternamelist")
def _e__docParamListItem__parameternamelist(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__docParamNameList(state, obj.parameternamelist.append, attr)


def _node_class_start__docParamListItem(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParamListItem.__new__(Node_docParamListItem)

    n.parameternamelist = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParamListItem,
            _node_class_finish__docParamListItem,
            None,
        )
    )


def _node_class_finish_fields__docParamListItem(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "parameterdescription"):
        _raise_missing_element_error(state, obj, "parameterdescription")


def _node_class_finish__docParamListItem(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docParamListItem(state, n)


_node_class_child__docParamNameList = _cur_list = {}


@_add_to_list("parametername")
def _e__docParamNameList__parametername(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__docParamName(state, obj.parametername.append, attr)


@_add_to_list("parametertype")
def _e__docParamNameList__parametertype(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__docParamType(state, obj.parametertype.append, attr)


def _node_class_start__docParamNameList(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParamNameList.__new__(Node_docParamNameList)

    n.parametername = []
    n.parametertype = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParamNameList,
            None,
            None,
        )
    )


_node_class_child__docParamType = _cur_list = {}


@_add_to_list("ref")
def _e__docParamType__ref(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refTextType(state, (lambda x: obj.append(TaggedValue("ref", x))), attr)


def _node_class_start__docParamType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParamType.__new__(Node_docParamType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParamType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docParamName = _cur_list = {}


@_add_to_list("direction")
def _a__docParamName__direction(state: _ParseState, obj, value: str, /):
    try:
        obj.direction = DoxParamDir(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__docParamName(state: _ParseState, obj, /):
    if not hasattr(obj, "direction"):
        obj.direction = None


_node_class_child__docParamName = _cur_list = {}


@_add_to_list("ref")
def _e__docParamName__ref(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__refTextType(state, (lambda x: obj.append(TaggedValue("ref", x))), attr)


def _node_class_start__docParamName(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParamName.__new__(Node_docParamName)

    for name, value in attr:
        handler = _node_class_attr__docParamName.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docParamName(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParamName,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__tableofcontentsKindType = _cur_list = {}


@_add_to_list("name")
def _e__tableofcontentsKindType__name(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "name"):
        _raise_duplicate_element_error(state, "name")

    _node_start_string(state, functools.partial(setattr, obj, "name"), attr)


@_add_to_list("reference")
def _e__tableofcontentsKindType__reference(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "reference"):
        _raise_duplicate_element_error(state, "reference")

    _node_start_string(state, functools.partial(setattr, obj, "reference"), attr)


@_add_to_list("tableofcontents")
def _e__tableofcontentsKindType__tableofcontents(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__tableofcontentsType(state, obj.tableofcontents.append, attr)


def _node_class_start__tableofcontentsKindType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_tableofcontentsKindType.__new__(Node_tableofcontentsKindType)

    n.tableofcontents = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__tableofcontentsKindType,
            _node_class_finish__tableofcontentsKindType,
            None,
        )
    )


def _node_class_finish_fields__tableofcontentsKindType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "name"):
        _raise_missing_element_error(state, obj, "name")
    if not hasattr(obj, "reference"):
        _raise_missing_element_error(state, obj, "reference")


def _node_class_finish__tableofcontentsKindType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__tableofcontentsKindType(state, n)


_node_class_attr__incType = _cur_list = {}


@_add_to_list("refid")
def _a__incType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("local")
def _a__incType__local(state: _ParseState, obj, value: str, /):
    obj.local = _parse_DoxBool_attribute(state, "local", value)


def _node_class_attr_end__incType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        obj.refid = None
    if not hasattr(obj, "local"):
        _raise_missing_attribute_error(state, "local")


def _node_class_start__incType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_incType.__new__(Node_incType)

    for name, value in attr:
        handler = _node_class_attr__incType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__incType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__compoundRefType = _cur_list = {}


@_add_to_list("refid")
def _a__compoundRefType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("prot")
def _a__compoundRefType__prot(state: _ParseState, obj, value: str, /):
    try:
        obj.prot = DoxProtectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("virt")
def _a__compoundRefType__virt(state: _ParseState, obj, value: str, /):
    try:
        obj.virt = DoxVirtualKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__compoundRefType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        obj.refid = None
    if not hasattr(obj, "prot"):
        _raise_missing_attribute_error(state, "prot")
    if not hasattr(obj, "virt"):
        _raise_missing_attribute_error(state, "virt")


def _node_class_start__compoundRefType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_compoundRefType.__new__(Node_compoundRefType)

    for name, value in attr:
        handler = _node_class_attr__compoundRefType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__compoundRefType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__refType = _cur_list = {}


@_add_to_list("refid")
def _a__refType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("prot")
def _a__refType__prot(state: _ParseState, obj, value: str, /):
    try:
        obj.prot = DoxProtectionKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("inline")
def _a__refType__inline(state: _ParseState, obj, value: str, /):
    obj.inline = _parse_DoxBool_attribute(state, "inline", value)


def _node_class_attr_end__refType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "prot"):
        obj.prot = None
    if not hasattr(obj, "inline"):
        obj.inline = None


def _node_class_start__refType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_refType.__new__(Node_refType)

    for name, value in attr:
        handler = _node_class_attr__refType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__refType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__exportType = _cur_list = {}


@_add_to_list("refid")
def _a__exportType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


def _node_class_attr_end__exportType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        obj.refid = None


def _node_class_start__exportType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_exportType.__new__(Node_exportType)

    for name, value in attr:
        handler = _node_class_attr__exportType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__exportType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__refTextType = _cur_list = {}


@_add_to_list("refid")
def _a__refTextType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("kindref")
def _a__refTextType__kindref(state: _ParseState, obj, value: str, /):
    try:
        obj.kindref = DoxRefKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("external")
def _a__refTextType__external(state: _ParseState, obj, value: str, /):
    obj.external = value


@_add_to_list("tooltip")
def _a__refTextType__tooltip(state: _ParseState, obj, value: str, /):
    obj.tooltip = value


def _node_class_attr_end__refTextType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "kindref"):
        _raise_missing_attribute_error(state, "kindref")
    if not hasattr(obj, "external"):
        obj.external = None
    if not hasattr(obj, "tooltip"):
        obj.tooltip = None


def _node_class_start__refTextType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_refTextType.__new__(Node_refTextType)

    for name, value in attr:
        handler = _node_class_attr__refTextType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__refTextType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__reimplementType = _cur_list = {}


@_add_to_list("refid")
def _a__reimplementType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


def _node_class_attr_end__reimplementType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")


def _node_class_start__reimplementType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_reimplementType.__new__(Node_reimplementType)

    for name, value in attr:
        handler = _node_class_attr__reimplementType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__reimplementType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docTitleCmdGroup = _cur_list = {}


@_add_to_list("ulink")
def _e__docTitleCmdGroup__ulink(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docURLLink(state, (lambda x: obj.append(TaggedValue("ulink", x))), attr)


@_add_to_list("bold")
def _e__docTitleCmdGroup__bold(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("bold", x))), attr)


@_add_to_list("s")
def _e__docTitleCmdGroup__s(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("s", x))), attr)


@_add_to_list("strike")
def _e__docTitleCmdGroup__strike(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("strike", x))), attr)


@_add_to_list("underline")
def _e__docTitleCmdGroup__underline(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(
        state, (lambda x: obj.append(TaggedValue("underline", x))), attr
    )


@_add_to_list("emphasis")
def _e__docTitleCmdGroup__emphasis(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(
        state, (lambda x: obj.append(TaggedValue("emphasis", x))), attr
    )


@_add_to_list("computeroutput")
def _e__docTitleCmdGroup__computeroutput(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__docMarkupType(
        state, (lambda x: obj.append(TaggedValue("computeroutput", x))), attr
    )


@_add_to_list("subscript")
def _e__docTitleCmdGroup__subscript(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(
        state, (lambda x: obj.append(TaggedValue("subscript", x))), attr
    )


@_add_to_list("superscript")
def _e__docTitleCmdGroup__superscript(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(
        state, (lambda x: obj.append(TaggedValue("superscript", x))), attr
    )


@_add_to_list("center")
def _e__docTitleCmdGroup__center(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("center", x))), attr)


@_add_to_list("small")
def _e__docTitleCmdGroup__small(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("small", x))), attr)


@_add_to_list("cite")
def _e__docTitleCmdGroup__cite(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("cite", x))), attr)


@_add_to_list("del")
def _e__docTitleCmdGroup__del(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("del", x))), attr)


@_add_to_list("ins")
def _e__docTitleCmdGroup__ins(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(state, (lambda x: obj.append(TaggedValue("ins", x))), attr)


@_add_to_list("htmlonly")
def _e__docTitleCmdGroup__htmlonly(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docHtmlOnlyType(
        state, (lambda x: obj.append(TaggedValue("htmlonly", x))), attr
    )


@_add_to_list("manonly")
def _e__docTitleCmdGroup__manonly(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("manonly", x))), attr)


@_add_to_list("xmlonly")
def _e__docTitleCmdGroup__xmlonly(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("xmlonly", x))), attr)


@_add_to_list("rtfonly")
def _e__docTitleCmdGroup__rtfonly(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("rtfonly", x))), attr)


@_add_to_list("latexonly")
def _e__docTitleCmdGroup__latexonly(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("latexonly", x))), attr)


@_add_to_list("docbookonly")
def _e__docTitleCmdGroup__docbookonly(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("docbookonly", x))), attr)


@_add_to_list("image")
def _e__docTitleCmdGroup__image(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docImageType(state, (lambda x: obj.append(TaggedValue("image", x))), attr)


@_add_to_list("dot")
def _e__docTitleCmdGroup__dot(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docDotMscType(state, (lambda x: obj.append(TaggedValue("dot", x))), attr)


@_add_to_list("msc")
def _e__docTitleCmdGroup__msc(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docDotMscType(state, (lambda x: obj.append(TaggedValue("msc", x))), attr)


@_add_to_list("plantuml")
def _e__docTitleCmdGroup__plantuml(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docPlantumlType(
        state, (lambda x: obj.append(TaggedValue("plantuml", x))), attr
    )


@_add_to_list("anchor")
def _e__docTitleCmdGroup__anchor(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docAnchorType(state, (lambda x: obj.append(TaggedValue("anchor", x))), attr)


@_add_to_list("formula")
def _e__docTitleCmdGroup__formula(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docFormulaType(
        state, (lambda x: obj.append(TaggedValue("formula", x))), attr
    )


@_add_to_list("ref")
def _e__docTitleCmdGroup__ref(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docRefTextType(state, (lambda x: obj.append(TaggedValue("ref", x))), attr)


@_add_to_list("emoji")
def _e__docTitleCmdGroup__emoji(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docEmojiType(state, (lambda x: obj.append(TaggedValue("emoji", x))), attr)


_add_to_list("linebreak")(_union_codepoint_element(chr(10)))
_add_to_list("nonbreakablespace")(_union_codepoint_element(chr(160)))
_add_to_list("iexcl")(_union_codepoint_element(chr(161)))
_add_to_list("cent")(_union_codepoint_element(chr(162)))
_add_to_list("pound")(_union_codepoint_element(chr(163)))
_add_to_list("curren")(_union_codepoint_element(chr(164)))
_add_to_list("yen")(_union_codepoint_element(chr(165)))
_add_to_list("brvbar")(_union_codepoint_element(chr(166)))
_add_to_list("sect")(_union_codepoint_element(chr(167)))
_add_to_list("umlaut")(_union_codepoint_element(chr(168)))
_add_to_list("copy")(_union_codepoint_element(chr(169)))
_add_to_list("ordf")(_union_codepoint_element(chr(170)))
_add_to_list("laquo")(_union_codepoint_element(chr(171)))
_add_to_list("not")(_union_codepoint_element(chr(172)))
_add_to_list("shy")(_union_codepoint_element(chr(173)))
_add_to_list("registered")(_union_codepoint_element(chr(174)))
_add_to_list("macr")(_union_codepoint_element(chr(175)))
_add_to_list("deg")(_union_codepoint_element(chr(176)))
_add_to_list("plusmn")(_union_codepoint_element(chr(177)))
_add_to_list("sup2")(_union_codepoint_element(chr(178)))
_add_to_list("sup3")(_union_codepoint_element(chr(179)))
_add_to_list("acute")(_union_codepoint_element(chr(180)))
_add_to_list("micro")(_union_codepoint_element(chr(181)))
_add_to_list("para")(_union_codepoint_element(chr(182)))
_add_to_list("middot")(_union_codepoint_element(chr(183)))
_add_to_list("cedil")(_union_codepoint_element(chr(184)))
_add_to_list("sup1")(_union_codepoint_element(chr(185)))
_add_to_list("ordm")(_union_codepoint_element(chr(186)))
_add_to_list("raquo")(_union_codepoint_element(chr(187)))
_add_to_list("frac14")(_union_codepoint_element(chr(188)))
_add_to_list("frac12")(_union_codepoint_element(chr(189)))
_add_to_list("frac34")(_union_codepoint_element(chr(190)))
_add_to_list("iquest")(_union_codepoint_element(chr(191)))
_add_to_list("Agrave")(_union_codepoint_element(chr(192)))
_add_to_list("Aacute")(_union_codepoint_element(chr(193)))
_add_to_list("Acirc")(_union_codepoint_element(chr(194)))
_add_to_list("Atilde")(_union_codepoint_element(chr(195)))
_add_to_list("Aumlaut")(_union_codepoint_element(chr(196)))
_add_to_list("Aring")(_union_codepoint_element(chr(197)))
_add_to_list("AElig")(_union_codepoint_element(chr(198)))
_add_to_list("Ccedil")(_union_codepoint_element(chr(199)))
_add_to_list("Egrave")(_union_codepoint_element(chr(200)))
_add_to_list("Eacute")(_union_codepoint_element(chr(201)))
_add_to_list("Ecirc")(_union_codepoint_element(chr(202)))
_add_to_list("Eumlaut")(_union_codepoint_element(chr(203)))
_add_to_list("Igrave")(_union_codepoint_element(chr(204)))
_add_to_list("Iacute")(_union_codepoint_element(chr(205)))
_add_to_list("Icirc")(_union_codepoint_element(chr(206)))
_add_to_list("Iumlaut")(_union_codepoint_element(chr(207)))
_add_to_list("ETH")(_union_codepoint_element(chr(208)))
_add_to_list("Ntilde")(_union_codepoint_element(chr(209)))
_add_to_list("Ograve")(_union_codepoint_element(chr(210)))
_add_to_list("Oacute")(_union_codepoint_element(chr(211)))
_add_to_list("Ocirc")(_union_codepoint_element(chr(212)))
_add_to_list("Otilde")(_union_codepoint_element(chr(213)))
_add_to_list("Oumlaut")(_union_codepoint_element(chr(214)))
_add_to_list("times")(_union_codepoint_element(chr(215)))
_add_to_list("Oslash")(_union_codepoint_element(chr(216)))
_add_to_list("Ugrave")(_union_codepoint_element(chr(217)))
_add_to_list("Uacute")(_union_codepoint_element(chr(218)))
_add_to_list("Ucirc")(_union_codepoint_element(chr(219)))
_add_to_list("Uumlaut")(_union_codepoint_element(chr(220)))
_add_to_list("Yacute")(_union_codepoint_element(chr(221)))
_add_to_list("THORN")(_union_codepoint_element(chr(222)))
_add_to_list("szlig")(_union_codepoint_element(chr(223)))
_add_to_list("agrave")(_union_codepoint_element(chr(224)))
_add_to_list("aacute")(_union_codepoint_element(chr(225)))
_add_to_list("acirc")(_union_codepoint_element(chr(226)))
_add_to_list("atilde")(_union_codepoint_element(chr(227)))
_add_to_list("aumlaut")(_union_codepoint_element(chr(228)))
_add_to_list("aring")(_union_codepoint_element(chr(229)))
_add_to_list("aelig")(_union_codepoint_element(chr(230)))
_add_to_list("ccedil")(_union_codepoint_element(chr(231)))
_add_to_list("egrave")(_union_codepoint_element(chr(232)))
_add_to_list("eacute")(_union_codepoint_element(chr(233)))
_add_to_list("ecirc")(_union_codepoint_element(chr(234)))
_add_to_list("eumlaut")(_union_codepoint_element(chr(235)))
_add_to_list("igrave")(_union_codepoint_element(chr(236)))
_add_to_list("iacute")(_union_codepoint_element(chr(237)))
_add_to_list("icirc")(_union_codepoint_element(chr(238)))
_add_to_list("iumlaut")(_union_codepoint_element(chr(239)))
_add_to_list("eth")(_union_codepoint_element(chr(240)))
_add_to_list("ntilde")(_union_codepoint_element(chr(241)))
_add_to_list("ograve")(_union_codepoint_element(chr(242)))
_add_to_list("oacute")(_union_codepoint_element(chr(243)))
_add_to_list("ocirc")(_union_codepoint_element(chr(244)))
_add_to_list("otilde")(_union_codepoint_element(chr(245)))
_add_to_list("oumlaut")(_union_codepoint_element(chr(246)))
_add_to_list("divide")(_union_codepoint_element(chr(247)))
_add_to_list("oslash")(_union_codepoint_element(chr(248)))
_add_to_list("ugrave")(_union_codepoint_element(chr(249)))
_add_to_list("uacute")(_union_codepoint_element(chr(250)))
_add_to_list("ucirc")(_union_codepoint_element(chr(251)))
_add_to_list("uumlaut")(_union_codepoint_element(chr(252)))
_add_to_list("yacute")(_union_codepoint_element(chr(253)))
_add_to_list("thorn")(_union_codepoint_element(chr(254)))
_add_to_list("yumlaut")(_union_codepoint_element(chr(255)))
_add_to_list("fnof")(_union_codepoint_element(chr(402)))
_add_to_list("Alpha")(_union_codepoint_element(chr(913)))
_add_to_list("Beta")(_union_codepoint_element(chr(914)))
_add_to_list("Gamma")(_union_codepoint_element(chr(915)))
_add_to_list("Delta")(_union_codepoint_element(chr(916)))
_add_to_list("Epsilon")(_union_codepoint_element(chr(917)))
_add_to_list("Zeta")(_union_codepoint_element(chr(918)))
_add_to_list("Eta")(_union_codepoint_element(chr(919)))
_add_to_list("Theta")(_union_codepoint_element(chr(920)))
_add_to_list("Iota")(_union_codepoint_element(chr(921)))
_add_to_list("Kappa")(_union_codepoint_element(chr(922)))
_add_to_list("Lambda")(_union_codepoint_element(chr(923)))
_add_to_list("Mu")(_union_codepoint_element(chr(924)))
_add_to_list("Nu")(_union_codepoint_element(chr(925)))
_add_to_list("Xi")(_union_codepoint_element(chr(926)))
_add_to_list("Omicron")(_union_codepoint_element(chr(927)))
_add_to_list("Pi")(_union_codepoint_element(chr(928)))
_add_to_list("Rho")(_union_codepoint_element(chr(929)))
_add_to_list("Sigma")(_union_codepoint_element(chr(931)))
_add_to_list("Tau")(_union_codepoint_element(chr(932)))
_add_to_list("Upsilon")(_union_codepoint_element(chr(933)))
_add_to_list("Phi")(_union_codepoint_element(chr(934)))
_add_to_list("Chi")(_union_codepoint_element(chr(935)))
_add_to_list("Psi")(_union_codepoint_element(chr(936)))
_add_to_list("Omega")(_union_codepoint_element(chr(937)))
_add_to_list("alpha")(_union_codepoint_element(chr(945)))
_add_to_list("beta")(_union_codepoint_element(chr(946)))
_add_to_list("gamma")(_union_codepoint_element(chr(947)))
_add_to_list("delta")(_union_codepoint_element(chr(948)))
_add_to_list("epsilon")(_union_codepoint_element(chr(949)))
_add_to_list("zeta")(_union_codepoint_element(chr(950)))
_add_to_list("eta")(_union_codepoint_element(chr(951)))
_add_to_list("theta")(_union_codepoint_element(chr(952)))
_add_to_list("iota")(_union_codepoint_element(chr(953)))
_add_to_list("kappa")(_union_codepoint_element(chr(954)))
_add_to_list("lambda")(_union_codepoint_element(chr(955)))
_add_to_list("mu")(_union_codepoint_element(chr(956)))
_add_to_list("nu")(_union_codepoint_element(chr(957)))
_add_to_list("xi")(_union_codepoint_element(chr(958)))
_add_to_list("omicron")(_union_codepoint_element(chr(959)))
_add_to_list("pi")(_union_codepoint_element(chr(960)))
_add_to_list("rho")(_union_codepoint_element(chr(961)))
_add_to_list("sigmaf")(_union_codepoint_element(chr(962)))
_add_to_list("sigma")(_union_codepoint_element(chr(963)))
_add_to_list("tau")(_union_codepoint_element(chr(964)))
_add_to_list("upsilon")(_union_codepoint_element(chr(965)))
_add_to_list("phi")(_union_codepoint_element(chr(966)))
_add_to_list("chi")(_union_codepoint_element(chr(967)))
_add_to_list("psi")(_union_codepoint_element(chr(968)))
_add_to_list("omega")(_union_codepoint_element(chr(969)))
_add_to_list("thetasym")(_union_codepoint_element(chr(977)))
_add_to_list("upsih")(_union_codepoint_element(chr(978)))
_add_to_list("piv")(_union_codepoint_element(chr(982)))
_add_to_list("bull")(_union_codepoint_element(chr(8226)))
_add_to_list("hellip")(_union_codepoint_element(chr(8230)))
_add_to_list("prime")(_union_codepoint_element(chr(8242)))
_add_to_list("Prime")(_union_codepoint_element(chr(8243)))
_add_to_list("oline")(_union_codepoint_element(chr(8254)))
_add_to_list("frasl")(_union_codepoint_element(chr(8260)))
_add_to_list("weierp")(_union_codepoint_element(chr(8472)))
_add_to_list("imaginary")(_union_codepoint_element(chr(8465)))
_add_to_list("real")(_union_codepoint_element(chr(8476)))
_add_to_list("trademark")(_union_codepoint_element(chr(8482)))
_add_to_list("alefsym")(_union_codepoint_element(chr(8501)))
_add_to_list("larr")(_union_codepoint_element(chr(8592)))
_add_to_list("uarr")(_union_codepoint_element(chr(8593)))
_add_to_list("rarr")(_union_codepoint_element(chr(8594)))
_add_to_list("darr")(_union_codepoint_element(chr(8595)))
_add_to_list("harr")(_union_codepoint_element(chr(8596)))
_add_to_list("crarr")(_union_codepoint_element(chr(8629)))
_add_to_list("lArr")(_union_codepoint_element(chr(8656)))
_add_to_list("uArr")(_union_codepoint_element(chr(8657)))
_add_to_list("rArr")(_union_codepoint_element(chr(8658)))
_add_to_list("dArr")(_union_codepoint_element(chr(8659)))
_add_to_list("hArr")(_union_codepoint_element(chr(8660)))
_add_to_list("forall")(_union_codepoint_element(chr(8704)))
_add_to_list("part")(_union_codepoint_element(chr(8706)))
_add_to_list("exist")(_union_codepoint_element(chr(8707)))
_add_to_list("empty")(_union_codepoint_element(chr(8709)))
_add_to_list("nabla")(_union_codepoint_element(chr(8711)))
_add_to_list("isin")(_union_codepoint_element(chr(8712)))
_add_to_list("notin")(_union_codepoint_element(chr(8713)))
_add_to_list("ni")(_union_codepoint_element(chr(8715)))
_add_to_list("prod")(_union_codepoint_element(chr(8719)))
_add_to_list("sum")(_union_codepoint_element(chr(8721)))
_add_to_list("minus")(_union_codepoint_element(chr(8722)))
_add_to_list("lowast")(_union_codepoint_element(chr(8727)))
_add_to_list("radic")(_union_codepoint_element(chr(8730)))
_add_to_list("prop")(_union_codepoint_element(chr(8733)))
_add_to_list("infin")(_union_codepoint_element(chr(8734)))
_add_to_list("ang")(_union_codepoint_element(chr(8736)))
_add_to_list("and")(_union_codepoint_element(chr(8743)))
_add_to_list("or")(_union_codepoint_element(chr(8744)))
_add_to_list("cap")(_union_codepoint_element(chr(8745)))
_add_to_list("cup")(_union_codepoint_element(chr(8746)))
_add_to_list("int")(_union_codepoint_element(chr(8747)))
_add_to_list("there4")(_union_codepoint_element(chr(8756)))
_add_to_list("sim")(_union_codepoint_element(chr(8764)))
_add_to_list("cong")(_union_codepoint_element(chr(8773)))
_add_to_list("asymp")(_union_codepoint_element(chr(8776)))
_add_to_list("ne")(_union_codepoint_element(chr(8800)))
_add_to_list("equiv")(_union_codepoint_element(chr(8801)))
_add_to_list("le")(_union_codepoint_element(chr(8804)))
_add_to_list("ge")(_union_codepoint_element(chr(8805)))
_add_to_list("sub")(_union_codepoint_element(chr(8834)))
_add_to_list("sup")(_union_codepoint_element(chr(8835)))
_add_to_list("nsub")(_union_codepoint_element(chr(8836)))
_add_to_list("sube")(_union_codepoint_element(chr(8838)))
_add_to_list("supe")(_union_codepoint_element(chr(8839)))
_add_to_list("oplus")(_union_codepoint_element(chr(8853)))
_add_to_list("otimes")(_union_codepoint_element(chr(8855)))
_add_to_list("perp")(_union_codepoint_element(chr(8869)))
_add_to_list("sdot")(_union_codepoint_element(chr(8901)))
_add_to_list("lceil")(_union_codepoint_element(chr(8968)))
_add_to_list("rceil")(_union_codepoint_element(chr(8969)))
_add_to_list("lfloor")(_union_codepoint_element(chr(8970)))
_add_to_list("rfloor")(_union_codepoint_element(chr(8971)))
_add_to_list("lang")(_union_codepoint_element(chr(10216)))
_add_to_list("rang")(_union_codepoint_element(chr(10217)))
_add_to_list("loz")(_union_codepoint_element(chr(9674)))
_add_to_list("spades")(_union_codepoint_element(chr(9824)))
_add_to_list("clubs")(_union_codepoint_element(chr(9827)))
_add_to_list("hearts")(_union_codepoint_element(chr(9829)))
_add_to_list("diams")(_union_codepoint_element(chr(9830)))
_add_to_list("OElig")(_union_codepoint_element(chr(338)))
_add_to_list("oelig")(_union_codepoint_element(chr(339)))
_add_to_list("Scaron")(_union_codepoint_element(chr(352)))
_add_to_list("scaron")(_union_codepoint_element(chr(353)))
_add_to_list("Yumlaut")(_union_codepoint_element(chr(376)))
_add_to_list("circ")(_union_codepoint_element(chr(710)))
_add_to_list("tilde")(_union_codepoint_element(chr(732)))
_add_to_list("ensp")(_union_codepoint_element(chr(8194)))
_add_to_list("emsp")(_union_codepoint_element(chr(8195)))
_add_to_list("thinsp")(_union_codepoint_element(chr(8201)))
_add_to_list("zwnj")(_union_codepoint_element(chr(8204)))
_add_to_list("zwj")(_union_codepoint_element(chr(8205)))
_add_to_list("lrm")(_union_codepoint_element(chr(8206)))
_add_to_list("rlm")(_union_codepoint_element(chr(8207)))
_add_to_list("ndash")(_union_codepoint_element(chr(8211)))
_add_to_list("mdash")(_union_codepoint_element(chr(8212)))
_add_to_list("lsquo")(_union_codepoint_element(chr(8216)))
_add_to_list("rsquo")(_union_codepoint_element(chr(8217)))
_add_to_list("sbquo")(_union_codepoint_element(chr(8218)))
_add_to_list("ldquo")(_union_codepoint_element(chr(8220)))
_add_to_list("rdquo")(_union_codepoint_element(chr(8221)))
_add_to_list("bdquo")(_union_codepoint_element(chr(8222)))
_add_to_list("dagger")(_union_codepoint_element(chr(8224)))
_add_to_list("Dagger")(_union_codepoint_element(chr(8225)))
_add_to_list("permil")(_union_codepoint_element(chr(8240)))
_add_to_list("lsaquo")(_union_codepoint_element(chr(8249)))
_add_to_list("rsaquo")(_union_codepoint_element(chr(8250)))
_add_to_list("euro")(_union_codepoint_element(chr(8364)))
_add_to_list("tm")(_union_codepoint_element(chr(8482)))
_node_class_child__docCmdGroup = _cur_list = {}
_node_class_child__docCmdGroup.update(_node_class_child__docTitleCmdGroup)


@_add_to_list("hruler")
def _e__docCmdGroup__hruler(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_empty(state, (lambda x: obj.append(TaggedValue("hruler", x))), attr)


@_add_to_list("preformatted")
def _e__docCmdGroup__preformatted(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docMarkupType(
        state, (lambda x: obj.append(TaggedValue("preformatted", x))), attr
    )


@_add_to_list("programlisting")
def _e__docCmdGroup__programlisting(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__listingType(
        state, (lambda x: obj.append(TaggedValue("programlisting", x))), attr
    )


@_add_to_list("verbatim")
def _e__docCmdGroup__verbatim(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("verbatim", x))), attr)


@_add_to_list("javadocliteral")
def _e__docCmdGroup__javadocliteral(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("javadocliteral", x))), attr)


@_add_to_list("javadoccode")
def _e__docCmdGroup__javadoccode(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, (lambda x: obj.append(TaggedValue("javadoccode", x))), attr)


@_add_to_list("indexentry")
def _e__docCmdGroup__indexentry(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docIndexEntryType(
        state, (lambda x: obj.append(TaggedValue("indexentry", x))), attr
    )


@_add_to_list("orderedlist")
def _e__docCmdGroup__orderedlist(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docListType(
        state, (lambda x: obj.append(TaggedValue("orderedlist", x))), attr
    )


@_add_to_list("itemizedlist")
def _e__docCmdGroup__itemizedlist(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docListType(
        state, (lambda x: obj.append(TaggedValue("itemizedlist", x))), attr
    )


@_add_to_list("simplesect")
def _e__docCmdGroup__simplesect(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSimpleSectType(
        state, (lambda x: obj.append(TaggedValue("simplesect", x))), attr
    )


@_add_to_list("title")
def _e__docCmdGroup__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docTitleType(state, (lambda x: obj.append(TaggedValue("title", x))), attr)


@_add_to_list("variablelist")
def _e__docCmdGroup__variablelist(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docVariableListType(
        state, (lambda x: obj.append(TaggedValue("variablelist", x))), attr
    )


@_add_to_list("table")
def _e__docCmdGroup__table(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docTableType(state, (lambda x: obj.append(TaggedValue("table", x))), attr)


@_add_to_list("heading")
def _e__docCmdGroup__heading(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docHeadingType(
        state, (lambda x: obj.append(TaggedValue("heading", x))), attr
    )


@_add_to_list("dotfile")
def _e__docCmdGroup__dotfile(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docImageFileType(
        state, (lambda x: obj.append(TaggedValue("dotfile", x))), attr
    )


@_add_to_list("mscfile")
def _e__docCmdGroup__mscfile(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docImageFileType(
        state, (lambda x: obj.append(TaggedValue("mscfile", x))), attr
    )


@_add_to_list("diafile")
def _e__docCmdGroup__diafile(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docImageFileType(
        state, (lambda x: obj.append(TaggedValue("diafile", x))), attr
    )


@_add_to_list("plantumlfile")
def _e__docCmdGroup__plantumlfile(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docImageFileType(
        state, (lambda x: obj.append(TaggedValue("plantumlfile", x))), attr
    )


@_add_to_list("toclist")
def _e__docCmdGroup__toclist(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docTocListType(
        state, (lambda x: obj.append(TaggedValue("toclist", x))), attr
    )


@_add_to_list("language")
def _e__docCmdGroup__language(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docLanguageType(
        state, (lambda x: obj.append(TaggedValue("language", x))), attr
    )


@_add_to_list("parameterlist")
def _e__docCmdGroup__parameterlist(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParamListType(
        state, (lambda x: obj.append(TaggedValue("parameterlist", x))), attr
    )


@_add_to_list("xrefsect")
def _e__docCmdGroup__xrefsect(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docXRefSectType(
        state, (lambda x: obj.append(TaggedValue("xrefsect", x))), attr
    )


@_add_to_list("copydoc")
def _e__docCmdGroup__copydoc(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docCopyType(state, (lambda x: obj.append(TaggedValue("copydoc", x))), attr)


@_add_to_list("details")
def _e__docCmdGroup__details(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docDetailsType(
        state, (lambda x: obj.append(TaggedValue("details", x))), attr
    )


@_add_to_list("blockquote")
def _e__docCmdGroup__blockquote(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docBlockQuoteType(
        state, (lambda x: obj.append(TaggedValue("blockquote", x))), attr
    )


@_add_to_list("parblock")
def _e__docCmdGroup__parblock(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParBlockType(
        state, (lambda x: obj.append(TaggedValue("parblock", x))), attr
    )


_node_class_child__docParaType = _cur_list = {}
_node_class_child__docParaType.update(_node_class_child__docCmdGroup)


def _node_class_start__docParaType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParaType.__new__(Node_docParaType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParaType,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docMarkupType = _cur_list = {}
_node_class_child__docMarkupType.update(_node_class_child__docCmdGroup)


def _node_class_start__docMarkupType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docMarkupType.__new__(Node_docMarkupType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docMarkupType,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docTitleType = _cur_list = {}
_node_class_child__docTitleType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docTitleType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docTitleType.__new__(Node_docTitleType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docTitleType,
            None,
            _node_list_common_text,
        )
    )


_node_class_child__docSummaryType = _cur_list = {}
_node_class_child__docSummaryType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docSummaryType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSummaryType.__new__(Node_docSummaryType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSummaryType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docURLLink = _cur_list = {}


@_add_to_list("url")
def _a__docURLLink__url(state: _ParseState, obj, value: str, /):
    obj.url = value


def _node_class_attr_end__docURLLink(state: _ParseState, obj, /):
    if not hasattr(obj, "url"):
        _raise_missing_attribute_error(state, "url")


_node_class_child__docURLLink = _cur_list = {}
_node_class_child__docURLLink.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docURLLink(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docURLLink.__new__(Node_docURLLink)

    for name, value in attr:
        handler = _node_class_attr__docURLLink.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docURLLink(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docURLLink,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docHtmlOnlyType = _cur_list = {}


@_add_to_list("block")
def _a__docHtmlOnlyType__block(state: _ParseState, obj, value: str, /):
    obj.block = value


def _node_class_attr_end__docHtmlOnlyType(state: _ParseState, obj, /):
    if not hasattr(obj, "block"):
        obj.block = None


def _node_class_start__docHtmlOnlyType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docHtmlOnlyType.__new__(Node_docHtmlOnlyType)

    for name, value in attr:
        handler = _node_class_attr__docHtmlOnlyType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docHtmlOnlyType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docImageType = _cur_list = {}


@_add_to_list("type")
def _a__docImageType__type(state: _ParseState, obj, value: str, /):
    try:
        obj.type = DoxImageKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


@_add_to_list("name")
def _a__docImageType__name(state: _ParseState, obj, value: str, /):
    obj.name = value


@_add_to_list("width")
def _a__docImageType__width(state: _ParseState, obj, value: str, /):
    obj.width = value


@_add_to_list("height")
def _a__docImageType__height(state: _ParseState, obj, value: str, /):
    obj.height = value


@_add_to_list("alt")
def _a__docImageType__alt(state: _ParseState, obj, value: str, /):
    obj.alt = value


@_add_to_list("inline")
def _a__docImageType__inline(state: _ParseState, obj, value: str, /):
    obj.inline = _parse_DoxBool_attribute(state, "inline", value)


@_add_to_list("caption")
def _a__docImageType__caption(state: _ParseState, obj, value: str, /):
    obj.caption = value


def _node_class_attr_end__docImageType(state: _ParseState, obj, /):
    if not hasattr(obj, "type"):
        obj.type = None
    if not hasattr(obj, "name"):
        obj.name = None
    if not hasattr(obj, "width"):
        obj.width = None
    if not hasattr(obj, "height"):
        obj.height = None
    if not hasattr(obj, "alt"):
        obj.alt = None
    if not hasattr(obj, "inline"):
        obj.inline = None
    if not hasattr(obj, "caption"):
        obj.caption = None


_node_class_child__docImageType = _cur_list = {}
_node_class_child__docImageType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docImageType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docImageType.__new__(Node_docImageType)

    for name, value in attr:
        handler = _node_class_attr__docImageType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docImageType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docImageType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docDotMscType = _cur_list = {}


@_add_to_list("name")
def _a__docDotMscType__name(state: _ParseState, obj, value: str, /):
    obj.name = value


@_add_to_list("width")
def _a__docDotMscType__width(state: _ParseState, obj, value: str, /):
    obj.width = value


@_add_to_list("height")
def _a__docDotMscType__height(state: _ParseState, obj, value: str, /):
    obj.height = value


@_add_to_list("caption")
def _a__docDotMscType__caption(state: _ParseState, obj, value: str, /):
    obj.caption = value


def _node_class_attr_end__docDotMscType(state: _ParseState, obj, /):
    if not hasattr(obj, "name"):
        obj.name = None
    if not hasattr(obj, "width"):
        obj.width = None
    if not hasattr(obj, "height"):
        obj.height = None
    if not hasattr(obj, "caption"):
        obj.caption = None


_node_class_child__docDotMscType = _cur_list = {}
_node_class_child__docDotMscType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docDotMscType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docDotMscType.__new__(Node_docDotMscType)

    for name, value in attr:
        handler = _node_class_attr__docDotMscType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docDotMscType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docDotMscType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docPlantumlType = _cur_list = {}


@_add_to_list("name")
def _a__docPlantumlType__name(state: _ParseState, obj, value: str, /):
    obj.name = value


@_add_to_list("width")
def _a__docPlantumlType__width(state: _ParseState, obj, value: str, /):
    obj.width = value


@_add_to_list("height")
def _a__docPlantumlType__height(state: _ParseState, obj, value: str, /):
    obj.height = value


@_add_to_list("caption")
def _a__docPlantumlType__caption(state: _ParseState, obj, value: str, /):
    obj.caption = value


@_add_to_list("engine")
def _a__docPlantumlType__engine(state: _ParseState, obj, value: str, /):
    try:
        obj.engine = DoxPlantumlEngine(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__docPlantumlType(state: _ParseState, obj, /):
    if not hasattr(obj, "name"):
        obj.name = None
    if not hasattr(obj, "width"):
        obj.width = None
    if not hasattr(obj, "height"):
        obj.height = None
    if not hasattr(obj, "caption"):
        obj.caption = None
    if not hasattr(obj, "engine"):
        obj.engine = None


_node_class_child__docPlantumlType = _cur_list = {}
_node_class_child__docPlantumlType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docPlantumlType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docPlantumlType.__new__(Node_docPlantumlType)

    for name, value in attr:
        handler = _node_class_attr__docPlantumlType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docPlantumlType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docPlantumlType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docRefTextType = _cur_list = {}


@_add_to_list("refid")
def _a__docRefTextType__refid(state: _ParseState, obj, value: str, /):
    obj.refid = value


@_add_to_list("kindref")
def _a__docRefTextType__kindref(state: _ParseState, obj, value: str, /):
    obj.kindref = value


@_add_to_list("external")
def _a__docRefTextType__external(state: _ParseState, obj, value: str, /):
    obj.external = value


def _node_class_attr_end__docRefTextType(state: _ParseState, obj, /):
    if not hasattr(obj, "refid"):
        _raise_missing_attribute_error(state, "refid")
    if not hasattr(obj, "kindref"):
        _raise_missing_attribute_error(state, "kindref")
    if not hasattr(obj, "external"):
        obj.external = None


_node_class_child__docRefTextType = _cur_list = {}
_node_class_child__docRefTextType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docRefTextType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docRefTextType.__new__(Node_docRefTextType)

    for name, value in attr:
        handler = _node_class_attr__docRefTextType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docRefTextType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docRefTextType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docHeadingType = _cur_list = {}


@_add_to_list("level")
def _a__docHeadingType__level(state: _ParseState, obj, value: str, /):
    try:
        obj.level = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


def _node_class_attr_end__docHeadingType(state: _ParseState, obj, /):
    if not hasattr(obj, "level"):
        _raise_missing_attribute_error(state, "level")


_node_class_child__docHeadingType = _cur_list = {}
_node_class_child__docHeadingType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docHeadingType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docHeadingType.__new__(Node_docHeadingType)

    for name, value in attr:
        handler = _node_class_attr__docHeadingType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docHeadingType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docHeadingType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docImageFileType = _cur_list = {}


@_add_to_list("name")
def _a__docImageFileType__name(state: _ParseState, obj, value: str, /):
    obj.name = value


@_add_to_list("width")
def _a__docImageFileType__width(state: _ParseState, obj, value: str, /):
    obj.width = value


@_add_to_list("height")
def _a__docImageFileType__height(state: _ParseState, obj, value: str, /):
    obj.height = value


def _node_class_attr_end__docImageFileType(state: _ParseState, obj, /):
    if not hasattr(obj, "name"):
        obj.name = None
    if not hasattr(obj, "width"):
        obj.width = None
    if not hasattr(obj, "height"):
        obj.height = None


_node_class_child__docImageFileType = _cur_list = {}
_node_class_child__docImageFileType.update(_node_class_child__docTitleCmdGroup)


def _node_class_start__docImageFileType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docImageFileType.__new__(Node_docImageFileType)

    for name, value in attr:
        handler = _node_class_attr__docImageFileType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docImageFileType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docImageFileType,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docAnchorType = _cur_list = {}


@_add_to_list("id")
def _a__docAnchorType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docAnchorType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


def _node_class_start__docAnchorType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docAnchorType.__new__(Node_docAnchorType)

    for name, value in attr:
        handler = _node_class_attr__docAnchorType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docAnchorType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docFormulaType = _cur_list = {}


@_add_to_list("id")
def _a__docFormulaType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docFormulaType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


def _node_class_start__docFormulaType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docFormulaType.__new__(Node_docFormulaType)

    for name, value in attr:
        handler = _node_class_attr__docFormulaType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docFormulaType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            _node_list_common_text,
        )
    )


_node_class_attr__docEmojiType = _cur_list = {}


@_add_to_list("name")
def _a__docEmojiType__name(state: _ParseState, obj, value: str, /):
    obj.name = value


@_add_to_list("unicode")
def _a__docEmojiType__unicode(state: _ParseState, obj, value: str, /):
    obj.unicode = value


def _node_class_attr_end__docEmojiType(state: _ParseState, obj, /):
    if not hasattr(obj, "name"):
        _raise_missing_attribute_error(state, "name")
    if not hasattr(obj, "unicode"):
        _raise_missing_attribute_error(state, "unicode")


def _node_class_start__docEmojiType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docEmojiType.__new__(Node_docEmojiType)

    for name, value in attr:
        handler = _node_class_attr__docEmojiType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docEmojiType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            None,
            None,
            None,
        )
    )


_node_class_child__docIndexEntryType = _cur_list = {}


@_add_to_list("primaryie")
def _e__docIndexEntryType__primaryie(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "primaryie"):
        _raise_duplicate_element_error(state, "primaryie")

    _node_start_string(state, functools.partial(setattr, obj, "primaryie"), attr)


@_add_to_list("secondaryie")
def _e__docIndexEntryType__secondaryie(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "secondaryie"):
        _raise_duplicate_element_error(state, "secondaryie")

    _node_start_string(state, functools.partial(setattr, obj, "secondaryie"), attr)


def _node_class_start__docIndexEntryType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docIndexEntryType.__new__(Node_docIndexEntryType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docIndexEntryType,
            _node_class_finish__docIndexEntryType,
            None,
        )
    )


def _node_class_finish_fields__docIndexEntryType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "primaryie"):
        _raise_missing_element_error(state, obj, "primaryie")
    if not hasattr(obj, "secondaryie"):
        _raise_missing_element_error(state, obj, "secondaryie")


def _node_class_finish__docIndexEntryType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docIndexEntryType(state, n)


_node_class_attr__docListType = _cur_list = {}


@_add_to_list("type")
def _a__docListType__type(state: _ParseState, obj, value: str, /):
    obj.type = _parse__DoxOlType(state, value)


@_add_to_list("start")
def _a__docListType__start(state: _ParseState, obj, value: str, /):
    try:
        obj.start = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


def _node_class_attr_end__docListType(state: _ParseState, obj, /):
    if not hasattr(obj, "type"):
        obj.type = None
    if not hasattr(obj, "start"):
        obj.start = None


_node_class_child__docListType = _cur_list = {}


@_add_to_list("listitem")
def _e__docListType__listitem(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docListItemType(state, obj.append, attr)


def _node_class_start__docListType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docListType.__new__(Node_docListType)

    for name, value in attr:
        handler = _node_class_attr__docListType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docListType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docListType,
            None,
            None,
        )
    )


_node_class_attr__docSimpleSectType = _cur_list = {}


@_add_to_list("kind")
def _a__docSimpleSectType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = DoxSimpleSectKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__docSimpleSectType(state: _ParseState, obj, /):
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")


_node_class_child__docSimpleSectType = _cur_list = {}


@_add_to_list("title")
def _e__docSimpleSectType__title(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "title"):
        _raise_duplicate_element_error(state, "title")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "title"), attr)


@_add_to_list("para")
def _e__docSimpleSectType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.para.append, attr)


def _node_class_start__docSimpleSectType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docSimpleSectType.__new__(Node_docSimpleSectType)

    n.para = []
    for name, value in attr:
        handler = _node_class_attr__docSimpleSectType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docSimpleSectType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docSimpleSectType,
            _node_class_finish__docSimpleSectType,
            None,
        )
    )


def _node_class_finish_fields__docSimpleSectType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "title"):
        obj.title = None
    if len(obj.para) < 1:
        _raise_empty_list_element_error(state, "para")


def _node_class_finish__docSimpleSectType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docSimpleSectType(state, n)


_node_class_child__docVariableListType = _cur_list = {}


@_add_to_list("varlistentry")
def _e__docVariableListType__varlistentry(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__docVarListEntryType(
        state,
        _push_tuple_item(
            state, 0, _tuple_item_tag_names__docVariableListType, ListItem_docVariableListType, obj
        ),
        attr,
    )


@_add_to_list("listitem")
def _e__docVariableListType__listitem(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docListItemType(
        state,
        _push_tuple_item(
            state, 1, _tuple_item_tag_names__docVariableListType, ListItem_docVariableListType, obj
        ),
        attr,
    )


_tuple_item_tag_names__docVariableListType = ListItem_docVariableListType._fields


def _node_class_start__docVariableListType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docVariableListType.__new__(Node_docVariableListType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docVariableListType,
            _node_class_finish__docVariableListType,
            None,
        )
    )


def _node_class_finish__docVariableListType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _check_complete_tuple(
        state, _tuple_item_tag_names__docVariableListType, ListItem_docVariableListType, n
    )


_node_class_attr__docTableType = _cur_list = {}


@_add_to_list("rows")
def _a__docTableType__rows(state: _ParseState, obj, value: str, /):
    try:
        obj.rows = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("cols")
def _a__docTableType__cols(state: _ParseState, obj, value: str, /):
    try:
        obj.cols = int(value, 10)
    except ValueError:
        _raise_invalid_int_error(state, value)


@_add_to_list("width")
def _a__docTableType__width(state: _ParseState, obj, value: str, /):
    obj.width = value


def _node_class_attr_end__docTableType(state: _ParseState, obj, /):
    if not hasattr(obj, "rows"):
        _raise_missing_attribute_error(state, "rows")
    if not hasattr(obj, "cols"):
        _raise_missing_attribute_error(state, "cols")
    if not hasattr(obj, "width"):
        obj.width = None


_node_class_child__docTableType = _cur_list = {}


@_add_to_list("caption")
def _e__docTableType__caption(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "caption"):
        _raise_duplicate_element_error(state, "caption")

    _node_class_start__docCaptionType(state, functools.partial(setattr, obj, "caption"), attr)


@_add_to_list("row")
def _e__docTableType__row(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docRowType(state, obj.row.append, attr)


def _node_class_start__docTableType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docTableType.__new__(Node_docTableType)

    n.row = []
    for name, value in attr:
        handler = _node_class_attr__docTableType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docTableType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docTableType,
            _node_class_finish__docTableType,
            None,
        )
    )


def _node_class_finish_fields__docTableType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "caption"):
        obj.caption = None


def _node_class_finish__docTableType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docTableType(state, n)


_node_class_child__docTocListType = _cur_list = {}


@_add_to_list("tocitem")
def _e__docTocListType__tocitem(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docTocItemType(state, obj.append, attr)


def _node_class_start__docTocListType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docTocListType.__new__(Node_docTocListType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docTocListType,
            None,
            None,
        )
    )


_node_class_attr__docLanguageType = _cur_list = {}


@_add_to_list("langid")
def _a__docLanguageType__langid(state: _ParseState, obj, value: str, /):
    obj.langid = value


def _node_class_attr_end__docLanguageType(state: _ParseState, obj, /):
    if not hasattr(obj, "langid"):
        _raise_missing_attribute_error(state, "langid")


_node_class_child__docLanguageType = _cur_list = {}


@_add_to_list("para")
def _e__docLanguageType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.append, attr)


def _node_class_start__docLanguageType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docLanguageType.__new__(Node_docLanguageType)

    for name, value in attr:
        handler = _node_class_attr__docLanguageType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docLanguageType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docLanguageType,
            None,
            None,
        )
    )


_node_class_attr__docParamListType = _cur_list = {}


@_add_to_list("kind")
def _a__docParamListType__kind(state: _ParseState, obj, value: str, /):
    try:
        obj.kind = DoxParamListKind(value.strip())
    except ValueError:
        _raise_invalid_enum_error(state, value)


def _node_class_attr_end__docParamListType(state: _ParseState, obj, /):
    if not hasattr(obj, "kind"):
        _raise_missing_attribute_error(state, "kind")


_node_class_child__docParamListType = _cur_list = {}


@_add_to_list("parameteritem")
def _e__docParamListType__parameteritem(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    _node_class_start__docParamListItem(state, obj.append, attr)


def _node_class_start__docParamListType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParamListType.__new__(Node_docParamListType)

    for name, value in attr:
        handler = _node_class_attr__docParamListType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docParamListType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParamListType,
            None,
            None,
        )
    )


_node_class_attr__docXRefSectType = _cur_list = {}


@_add_to_list("id")
def _a__docXRefSectType__id(state: _ParseState, obj, value: str, /):
    obj.id = value


def _node_class_attr_end__docXRefSectType(state: _ParseState, obj, /):
    if not hasattr(obj, "id"):
        _raise_missing_attribute_error(state, "id")


_node_class_child__docXRefSectType = _cur_list = {}


@_add_to_list("xreftitle")
def _e__docXRefSectType__xreftitle(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_start_string(state, obj.xreftitle.append, attr)


@_add_to_list("xrefdescription")
def _e__docXRefSectType__xrefdescription(
    state: _ParseState, obj, attr: Iterable[tuple[str, str]], /
):
    if hasattr(obj, "xrefdescription"):
        _raise_duplicate_element_error(state, "xrefdescription")

    _node_class_start__descriptionType(
        state, functools.partial(setattr, obj, "xrefdescription"), attr
    )


def _node_class_start__docXRefSectType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docXRefSectType.__new__(Node_docXRefSectType)

    n.xreftitle = []
    for name, value in attr:
        handler = _node_class_attr__docXRefSectType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docXRefSectType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docXRefSectType,
            _node_class_finish__docXRefSectType,
            None,
        )
    )


def _node_class_finish_fields__docXRefSectType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "xrefdescription"):
        _raise_missing_element_error(state, obj, "xrefdescription")


def _node_class_finish__docXRefSectType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docXRefSectType(state, n)


_node_class_attr__docCopyType = _cur_list = {}


@_add_to_list("link")
def _a__docCopyType__link(state: _ParseState, obj, value: str, /):
    obj.link = value


def _node_class_attr_end__docCopyType(state: _ParseState, obj, /):
    if not hasattr(obj, "link"):
        _raise_missing_attribute_error(state, "link")


_node_class_child__docCopyType = _cur_list = {}


@_add_to_list("para")
def _e__docCopyType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.para.append, attr)


@_add_to_list("sec1")
def _e__docCopyType__sec1(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docSect1Type(state, obj.sec1.append, attr)


@_add_to_list("internal")
def _e__docCopyType__internal(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "internal"):
        _raise_duplicate_element_error(state, "internal")

    _node_class_start__docInternalType(state, functools.partial(setattr, obj, "internal"), attr)


def _node_class_start__docCopyType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docCopyType.__new__(Node_docCopyType)

    n.para = []
    n.sec1 = []
    for name, value in attr:
        handler = _node_class_attr__docCopyType.get(name)

        if handler is not None:
            handler(state, n, value)
        else:
            _warn_unexpected_attribute(state, name)
    _node_class_attr_end__docCopyType(state, n)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docCopyType,
            _node_class_finish__docCopyType,
            None,
        )
    )


def _node_class_finish_fields__docCopyType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "internal"):
        obj.internal = None


def _node_class_finish__docCopyType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docCopyType(state, n)


_node_class_child__docDetailsType = _cur_list = {}


@_add_to_list("summary")
def _e__docDetailsType__summary(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "summary"):
        _raise_duplicate_element_error(state, "summary")

    _node_class_start__docSummaryType(state, functools.partial(setattr, obj, "summary"), attr)


@_add_to_list("para")
def _e__docDetailsType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.para.append, attr)


def _node_class_start__docDetailsType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docDetailsType.__new__(Node_docDetailsType)

    n.para = []
    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docDetailsType,
            _node_class_finish__docDetailsType,
            None,
        )
    )


def _node_class_finish_fields__docDetailsType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "summary"):
        obj.summary = None


def _node_class_finish__docDetailsType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docDetailsType(state, n)


_node_class_child__docBlockQuoteType = _cur_list = {}


@_add_to_list("para")
def _e__docBlockQuoteType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.append, attr)


def _node_class_start__docBlockQuoteType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docBlockQuoteType.__new__(Node_docBlockQuoteType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docBlockQuoteType,
            None,
            None,
        )
    )


_node_class_child__docParBlockType = _cur_list = {}


@_add_to_list("para")
def _e__docParBlockType__para(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    _node_class_start__docParaType(state, obj.append, attr)


def _node_class_start__docParBlockType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docParBlockType.__new__(Node_docParBlockType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docParBlockType,
            None,
            None,
        )
    )


_node_class_child__docVarListEntryType = _cur_list = {}


@_add_to_list("term")
def _e__docVarListEntryType__term(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    if hasattr(obj, "term"):
        _raise_duplicate_element_error(state, "term")

    _node_class_start__docTitleType(state, functools.partial(setattr, obj, "term"), attr)


def _node_class_start__docVarListEntryType(
    state: _ParseState, setter: Callable, attr: Iterable[tuple[str, str]], /
):
    n = Node_docVarListEntryType.__new__(Node_docVarListEntryType)

    for name, value in attr:
        _warn_unexpected_attribute(state, name)
    state.parse_callbacks.append(
        _ParseCallbacks(
            n,
            setter,
            _node_class_child__docVarListEntryType,
            _node_class_finish__docVarListEntryType,
            None,
        )
    )


def _node_class_finish_fields__docVarListEntryType(state: _ParseState, obj, /) -> None:
    if not hasattr(obj, "term"):
        _raise_missing_element_error(state, obj, "term")


def _node_class_finish__docVarListEntryType(state: _ParseState, /):
    n = state.parse_callbacks[-1].value
    _node_class_finish_fields__docVarListEntryType(state, n)


def _parse__DoxOlType(state: _ParseState, data, /):
    data = data.strip()
    if len(data) != 1:
        state.raise_parse_error("value must be a single character")

    if data not in "1aAiI":
        _raise_invalid_char_enum_error(state, data, "1aAiI")

    return data


_top_level_handlers = _cur_list = {}


@_add_to_list("doxygen")
def _(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    cb = state.parse_callbacks[-1]
    if obj is not None:
        state.raise_parse_error("cannot have more than one root element")

    def setter(x):
        cb.value = TaggedValue("doxygen", x)

    return _node_class_start__DoxygenType(state, setter, attr)


@_add_to_list("doxygenindex")
def _(state: _ParseState, obj, attr: Iterable[tuple[str, str]], /):
    cb = state.parse_callbacks[-1]
    if obj is not None:
        state.raise_parse_error("cannot have more than one root element")

    def setter(x):
        cb.value = TaggedValue("doxygenindex", x)

    return _node_class_start__DoxygenTypeIndex(state, setter, attr)


def _parse(obj, meth, /):
    p = expat.ParserCreate()
    state = _ParseState(p)

    state.parse_callbacks.append(_ParseCallbacks(None, None, _top_level_handlers))

    p.StartElementHandler = state.start_element
    p.EndElementHandler = state.end_element
    p.CharacterDataHandler = state.character_data

    try:
        meth(p, obj)
    except expat.ExpatError as e:
        raise ParseError(expat.errors.messages[e.code], e.lineno)
    finally:
        # break reference cycle for faster garbage collection
        p.StartElementHandler = None
        p.EndElementHandler = None
        p.CharacterDataHandler = None

    value = state.parse_callbacks[0].value
    if value is None:
        raise ParseError("document without a recognized root element", None)

    return value
