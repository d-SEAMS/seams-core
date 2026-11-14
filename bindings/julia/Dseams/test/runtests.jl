using Test
using Dseams

# The counts and the relative tolerance are the R front checks in
# bindings/r/dseams/tests/testthat/test-ice.R.

const tip = Dseams.example_path("exampleTraj.lammpstrj")
const nucleation = Dseams.example_path("nucleation_f1.lammpstrj")
const root = get(ENV, "DSEAMS_ROOT", "")

@testset "library and examples" begin
    @test !isempty(Dseams.library_path())
    @test VersionNumber(Dseams.engine_version()) >= v"2.14.0"
    @test_throws ErrorException Dseams.example_path("missing.lammpstrj")
    @test_throws ErrorException read_frame(joinpath(mktempdir(), "missing.lammpstrj"))
end

@testset "ice Ih" begin
    frame = read_frame(Dseams.ice_ih_path(); type = -1)
    @test nop(frame) == 128
    @test nframes(Dseams.ice_ih_path()) == 1
    @test chill_plus(frame; graph = "cutoff").counts["hexagonal"] == 128
    @test primitive_rings(frame; graph = "cutoff", max_depth = 6)[6] == 256
    @test cages(frame; graph = "cutoff").n_ih == 128
end

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
    @test neighbour_pairs(frame) > 0
end

@testset "seeded cages match with offload" begin
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
    order = q4_q6(frame)
    @test length(order.q6) == nop(frame)
    tau = jump_tau90(frame, frame, 1.0)
    @test isnan(tau) || isfinite(tau)
end

# A seams-core checkout has the larger trajectories, and the originals of
# the copies in data/.
if !isempty(root)
    @testset "checkout data" begin
        for name in readdir(Dseams.DATA)
            original = joinpath(root, "input", "traj", name)
            @test read(original) == read(joinpath(Dseams.DATA, name))
        end
        cubic = joinpath(root, "input", "traj", "mW_cubic.lammpstrj")
        frame = read_frame(cubic; type = 1)
        @test chill_plus(frame).counts["cubic"] == 4096
        @test cages(frame).n_ic == 4096
        @test nframes(cubic) >= 1
    end
end
