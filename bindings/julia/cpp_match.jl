# Copyright (c) d-SEAMS developers
# Compare the seams.h ccall wrappers with the C++ engine census.

function _oracle_rows()
    oracle = get(ENV, "SEAMS_CPP_ORACLE", "")
    root = get(ENV, "SEAMS_SOURCE_ROOT", get(ENV, "DSEAMS_ROOT", ""))
    @assert !isempty(oracle)
    @assert !isempty(root)
    text = read(`$oracle $root`, String)
    rows = Dict{String,Vector{Float64}}()
    for line in split(strip(text), '\n')
        parts = split(line)
        rows[parts[1] * " " * parts[2]] = parse.(Float64, parts[3:end])
    end
    return root, rows
end

function _close(got, want)
    scale = max(abs(got), abs(want), 1.0)
    @assert abs(got - want) <= 1e-8 * scale
end

function _mean(values)
    isempty(values) && return 0.0
    sum(values) / length(values)
end

function check_seams_against_cpp()
    root, want = _oracle_rows()
    ih = seams_system_read(System(), joinpath(root, "input/traj/ice1h.lammpstrj"), 1, -1)
    @test seams_system_count(ih) == want["ih n"][1]
    labels, counts = seams_chill(ih)
    @test all(Float64.(counts) .== want["ih chill0"])
    q3, _ = seams_steinhardt(ih, 3)
    q4, _ = seams_steinhardt(ih, 4)
    q6, bar = seams_steinhardt(ih, 6; average=true)
    q8, _ = seams_steinhardt(ih, 8)
    q12, _ = seams_steinhardt(ih, 12)
    _close(_mean(q3), want["ih q3"][1])
    _close(_mean(q4), want["ih q4"][1])
    _close(_mean(q6), want["ih q6"][1])
    _close(_mean(bar), want["ih q6bar"][1])
    _close(_mean(q8), want["ih q8"][1])
    _close(_mean(q12), want["ih q12"][1])
    vert = seams_rings_vertices(ih)
    @test length(vert) == want["ih vert"][1]
    @test (isempty(vert) ? -1 : Int(vert[1])) == want["ih vert"][2]
    @test all(Float64.(collect(seams_cages(ih))) .== want["ih cages"])
    @test seams_neighbours_knn(ih) == want["ih knn"][1]
    r, g, volume = seams_rdf(ih)
    @test length(r) == want["ih rdf"][1]
    _close(volume, want["ih rdf"][2])
    _close(isempty(g) ? 0.0 : g[1], want["ih rdf"][3])
    cn = seams_running_cn(ih)
    _close(isempty(cn) ? 0.0 : cn[end], want["ih rcn"][1])
    _close(seams_coordination(ih), want["ih cn"][1])
    @test seams_hbonds(ih) == want["ih hbond"][1]
    rho = seams_density_axis(ih)
    @test length(rho) == want["ih rho"][1]
    _close(isempty(rho) ? 0.0 : rho[1], want["ih rho"][2])
    key_len, ring_sum = seams_fingerprint(ih)
    @test key_len == want["ih fp"][1]
    @test ring_sum == want["ih fp"][2]
    @test seams_site_table(ih, "1=waterO")[1] == want["ih site"][1]
    @test all(Float64.(collect(seams_contact_pairs(ih, "1=cationHead"))) .== want["ih pairs"])
    n_subset, largest, percolation = seams_domains(ih, "1=waterO")
    @test n_subset == want["ih domain"][1]
    @test largest == want["ih domain"][2]
    _close(percolation, want["ih domain"][3])
    @test seams_prisms(ih)[7] == want["ih prism"][1]
    @test seams_monolayer(ih)[7] == want["ih sheet"][1]

    guest = seams_system_read(System(), joinpath(root, "input/traj/genice_sI.lammpstrj"), 1, 1)
    @test seams_system_count(guest) == want["si n"][1]
    @test all(Float64.(collect(seams_signature_cages(guest, "512"))) .== want["si sig"])

    xyz = seams_system_read(System(), joinpath(root, "input/xyz/configuration.xyz"), 1, -1)
    @test seams_system_count(xyz) == want["xyz n"][1]
    _close(seams_coordination(xyz; rmax=6.0), want["xyz cn"][1])

    con = seams_system_read(System(), joinpath(root, "input/con/tiny_multi_cuh2.con"), 1, -1)
    @test seams_system_count(con) == want["con n"][1]
    @test seams_hbonds(con; oxygen_type=0) == want["con hbond"][1]

    tri = seams_system_read(System(), joinpath(root, "input/traj/genice_sH.lammpstrj"), 1, 1)
    @test seams_system_count(tri) == want["tri n"][1]
    q6, _ = seams_steinhardt(tri, 6)
    _close(_mean(q6), want["tri q6"][1])
end
