# Copyright (c) d-SEAMS developers
"""1.x Lua globals, later Lua names, and the 2.14 call that replaces each.

``V1_GLOBALS`` is every function ``src/main.cpp`` registers with
``lua.set_function`` in v1.0.0 and v1.0.1; the two tags register the
same set. 1.x had no ``dseams`` table and no ``require("dseams")``.

``LUA_RELEASES`` are the yodaStruct releases that first registered a
later Lua name. No seams-core 2.x tag registers Lua names.

Each row of ``PAIRS`` holds the exact v1.0.1 globals (empty when the
analysis is new in 2.x) and the later names keyed by the release that
first shipped them. The migration page names every string here. The
tutorial test runs each kind on one input and checks the two 2.14 sides
agree.
"""

V1_TAG = "v1.0.1"

V1_GLOBALS = frozenset(
    (
        "averageQ6",
        "bondNetworkByIndex",
        "bulkTopologicalNetworkCriterion",
        "calcRDF",
        "chill_cij",
        "chill_iceType",
        "chillPlus_cij",
        "chillPlus_iceType",
        "clusterAnalysis",
        "getHbondNetwork",
        "getPrimitiveRings",
        "modifyChill",
        "neighborList",
        "percentage_Ice",
        "prismAnalysis",
        "readFrame",
        "readFrameOnlyOne",
        "recenterCluster",
        "ringAnalysis",
        "writeDump",
        "writeHistogram",
    )
)

LUA_RELEASES = ("2.0.0", "2.2.0", "2.4.1", "2.7.0", "2.8.0", "2.9.0")

PAIRS = (
    {
        "v1": ("readFrame", "readFrameOnlyOne"),
        "lua": {
            "2.0.0": ("readLammpsTrjO",),
            "2.2.0": ("dseams.read",),
            "2.7.0": ("readLammpsTrj",),
        },
        "cli": "seams read",
        "c_api": "seams_frame_read",
        "kind": "read",
    },
    {
        "v1": ("chillPlus_cij", "chillPlus_iceType"),
        "lua": {
            "2.0.0": ("getCorrelPlus", "getIceTypePlus"),
            "2.2.0": ("dseams.chill_plus",),
        },
        "cli": "seams chill-plus",
        "c_api": "seams_frame_chill_plus",
        "kind": "chill_plus",
    },
    {
        "v1": ("chill_cij", "chill_iceType"),
        "lua": {
            "2.0.0": ("getCorrel", "getIceType"),
            "2.2.0": ("dseams.chill",),
        },
        "cli": "seams chill",
        "c_api": "seams_frame_chill",
        "kind": "chill",
    },
    {
        "v1": ("getPrimitiveRings",),
        "lua": {"2.0.0": ("ringNetwork",)},
        "cli": "seams rings",
        "c_api": "seams_frame_rings",
        "kind": "rings",
    },
    {
        "v1": ("averageQ6",),
        "lua": {"2.0.0": ("steinhardtQl",)},
        "cli": "seams steinhardt",
        "c_api": "seams_frame_steinhardt_graph",
        "kind": "steinhardt",
    },
    {
        "v1": ("bulkTopologicalNetworkCriterion",),
        "lua": {
            "2.0.0": ("seededCageAffiliation",),
            "2.2.0": ("dseams.cages",),
        },
        "cli": "seams cages",
        "c_api": "seams_frame_cages",
        "kind": "cages",
    },
    {
        "v1": (),
        "lua": {"2.9.0": ("findBySignature", "dseams.cages_by_signature")},
        "cli": "seams cages --signature",
        "c_api": "seams_frame_signature_edge",
        "kind": "signature",
    },
    {
        "v1": ("neighborList", "bondNetworkByIndex"),
        "lua": {
            "2.0.0": ("neighListO",),
            "2.2.0": ("dseams.neighbors",),
        },
        "cli": "seams_frame_neighbours",
        "c_api": "seams_neighbours_build",
        "kind": "neighbors",
    },
    {
        "v1": (),
        "lua": {
            "2.0.0": ("kNearestNeighbourList",),
            "2.2.0": ("dseams.knn",),
        },
        "cli": "seams_frame_neighbours",
        "c_api": "seams_neighbours_build",
        "kind": "knn",
    },
    {
        "v1": ("calcRDF",),
        "lua": {
            "2.4.1": ("calcRDF3D",),
            "2.7.0": ("dseams.rdf",),
        },
        "cli": "seams rdf",
        "c_api": "seams_frame_rdf",
        "kind": "rdf",
    },
    {
        "v1": ("prismAnalysis",),
        "lua": {},
        "cli": "seams prisms",
        "c_api": "seams_frame_prisms",
        "kind": "prisms",
    },
    {
        "v1": ("ringAnalysis",),
        "lua": {},
        "cli": "seams_frame_monolayer",
        "c_api": "seams_frame_monolayer",
        "kind": "monolayer",
    },
    {
        "v1": ("getHbondNetwork",),
        "lua": {"2.7.0": ("dseams.hbonds",)},
        "cli": "seams hbonds",
        "c_api": "seams hbonds",
        "kind": "hbonds",
    },
    {
        "v1": (),
        "lua": {"2.4.1": ("calcCN", "dseams.cn")},
        "cli": "seams cn",
        "c_api": "seams cn",
        "kind": "cn",
    },
    {
        "v1": (),
        "lua": {"2.7.0": ("densityByType", "dseams.density")},
        "cli": "seams density-z",
        "c_api": "seams density-z",
        "kind": "density",
    },
    {
        "v1": (),
        "lua": {"2.8.0": ("topologyFingerprint", "dseams.fingerprint")},
        "cli": "seams fingerprint",
        "c_api": "seams fingerprint",
        "kind": "fingerprint",
    },
    {
        "v1": (),
        "lua": {"2.7.0": ("contactPairs", "dseams.pairs")},
        "cli": "seams pairs",
        "c_api": "seams pairs",
        "kind": "pairs",
    },
    {
        "v1": (),
        "lua": {"2.7.0": ("domainStats", "dseams.domain")},
        "cli": "seams domains",
        "c_api": "seams domains",
        "kind": "domains",
    },
)


def first_release(pair):
    """The yodaStruct release of the row's earliest later Lua name, or None."""
    for release in LUA_RELEASES:
        if release in pair["lua"]:
            return release
    return None


def v1_cell(pair):
    """The 1.x column text for a row."""
    if pair["v1"]:
        return " / ".join("~{}~".format(name) for name in pair["v1"])
    return "new in 2.x (yodaStruct {})".format(first_release(pair))


def lua_cell(pair, release):
    """The text under one yodaStruct release column for a row."""
    return " / ".join("~{}~".format(name) for name in pair["lua"].get(release, ()))


def table_rows(pair):
    """The cells of one migration table row, in column order."""
    return (
        v1_cell(pair),
        *(lua_cell(pair, release) for release in LUA_RELEASES),
        "~{}~".format(pair["cli"]),
        "~{}~".format(pair["c_api"]),
    )
