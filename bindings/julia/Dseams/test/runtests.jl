using Test
using Dseams

# The counts and the relative tolerance are the R front checks in
# bindings/r/dseams/tests/testthat/test-ice.R.

const root = get(ENV, "DSEAMS_ROOT", "")
@test !isempty(root)

const tip = joinpath(root, "input", "traj", "exampleTraj.lammpstrj")
const cubic = joinpath(root, "input", "traj", "mW_cubic.lammpstrj")
const nucleation = joinpath(root, "input", "traj", "nucleation_f1.lammpstrj")

@testset "CHILL+ on the mixed frame" begin
    frame = read_frame(tip; type = 2)
    @test nop(frame) == 250
    ice = chill_plus(frame)
    @test get(ice.counts, "clathrate", 0) == 0
    @test ice.counts["interClathrate"] == 9
    @test ice.counts["water"] == 241
    cutoff = chill_plus(frame; graph = "cutoff")
    @test ice.labels == cutoff.labels
    mutual = chill_plus(frame; graph = "mutual")
    @test mutual.counts["clathrate"] == 1
    @test mutual.counts["interClathrate"] == 7
    @test mutual.counts["water"] == 242
    rings = primitive_rings(frame; max_depth = 8)
    @test rings[6] == 224
    sig = signature_cages(frame, "512")
    @test sig.cages >= 0 && sig.atoms >= 0
    labels = symmetry(frame)
    @test length(labels) == nop(frame)
    q3 = steinhardt(frame, 3)
    finite = filter(isfinite, q3.ql)
    got = sum(finite) / length(finite)
    # Orthorhombic frames use the cutoff shell. 0.389362 is the four-nearest mean.
    @test isapprox(got, 0.374897; rtol = 1e-4)
end

@testset "seeded cages match with offload" begin
    @test isfile(nucleation)
    frame = read_frame(nucleation; type = 1)
    delete!(ENV, "SEAMS_OFFLOAD")
    host = cages(frame; graph = "seeded")
    ENV["SEAMS_OFFLOAD"] = "1"
    off = cages(frame; graph = "seeded")
    delete!(ENV, "SEAMS_OFFLOAD")
    @test off.cage == host.cage
    @test off.n_ih == host.n_ih
    @test off.n_ic == host.n_ic
    @test off.n_both == host.n_both
    @test off.n_water == host.n_water
    @test host.n_ih + host.n_ic + host.n_both + host.n_water == length(host.cage)
    @test host.n_both > 0
end

@testset "cubic ice is 4096 cubic and 4096 double-diamond" begin
    frame = read_frame(cubic; type = 1)
    ice = chill_plus(frame)
    @test ice.counts["cubic"] == 4096
    found = cages(frame)
    @test found.n_ic == 4096
end
