# Copyright (c) d-SEAMS developers
"""1.x entry points and the 2.14 call that replaces each one.

The migration page names every string in this list. The tutorial test
runs each kind on one input and checks the two 2.14 sides agree.
"""

PAIRS = (
    {
        "old": "dseams.read",
        "also": ("readLammpsTrjO", "readLammpsTrj"),
        "cli": "seams read",
        "c_api": "seams_frame_read",
        "kind": "read",
    },
    {
        "old": "dseams.chill_plus",
        "also": ("getIceTypePlus", "getCorrelPlus"),
        "cli": "seams chill-plus",
        "c_api": "seams_frame_chill_plus",
        "kind": "chill_plus",
    },
    {
        "old": "dseams.chill",
        "also": ("getIceType", "getCorrel"),
        "cli": "seams chill",
        "c_api": "seams_frame_chill",
        "kind": "chill",
    },
    {
        "old": "ringNetwork",
        "also": ("bulkPrimitiveRings",),
        "cli": "seams rings",
        "c_api": "seams_frame_rings",
        "kind": "rings",
    },
    {
        "old": "steinhardtQl",
        "also": ("bulk.bondOrderParameters",),
        "cli": "seams steinhardt",
        "c_api": "seams_frame_steinhardt_graph",
        "kind": "steinhardt",
    },
    {
        "old": "dseams.cages",
        "also": ("seededCageAffiliation", "bulk.topologicalNetworkCriterion"),
        "cli": "seams cages",
        "c_api": "seams_frame_cages",
        "kind": "cages",
    },
    {
        "old": "dseams.cages_by_signature",
        "also": ("findBySignature",),
        "cli": "seams cages --signature",
        "c_api": "seams_frame_signature_edge",
        "kind": "signature",
    },
    {
        "old": "dseams.neighbors",
        "also": ("neighListO",),
        "cli": "seams_frame_neighbours",
        "c_api": "seams_neighbours_build",
        "kind": "neighbors",
    },
    {
        "old": "kNearestNeighbourList",
        "also": ("dseams.knn",),
        "cli": "seams_frame_neighbours",
        "c_api": "seams_neighbours_build",
        "kind": "knn",
    },
    {
        "old": "dseams.rdf",
        "also": ("calcRDF3D",),
        "cli": "seams rdf",
        "c_api": "seams_frame_rdf",
        "kind": "rdf",
    },
    {
        "old": "topoOneDim",
        "also": ("dseams prisms",),
        "cli": "seams prisms",
        "c_api": "seams_frame_prisms",
        "kind": "prisms",
    },
    {
        "old": "topoTwoDim",
        "also": ("shapeMatching",),
        "cli": "seams_frame_monolayer",
        "c_api": "seams_frame_monolayer",
        "kind": "monolayer",
    },
    {
        "old": "dseams.hbonds",
        "also": ("getHbondNetwork",),
        "cli": "seams hbonds",
        "c_api": "seams hbonds",
        "kind": "hbonds",
    },
    {
        "old": "dseams.cn",
        "also": ("calcCN",),
        "cli": "seams cn",
        "c_api": "seams cn",
        "kind": "cn",
    },
    {
        "old": "dseams.density",
        "also": ("densityByType",),
        "cli": "seams density-z",
        "c_api": "seams density-z",
        "kind": "density",
    },
    {
        "old": "dseams.fingerprint",
        "also": ("topologyFingerprint",),
        "cli": "seams fingerprint",
        "c_api": "seams fingerprint",
        "kind": "fingerprint",
    },
    {
        "old": "dseams.pairs",
        "also": ("contactPairs",),
        "cli": "seams pairs",
        "c_api": "seams pairs",
        "kind": "pairs",
    },
    {
        "old": "dseams.domain",
        "also": ("domainStats",),
        "cli": "seams domains",
        "c_api": "seams domains",
        "kind": "domains",
    },
)
