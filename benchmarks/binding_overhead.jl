# Copyright (c) d-SEAMS developers
# Times the Julia ccall wrappers. Commands are "threads op path".

using Seams

const sink = Ref(0.0)

function retune(threads::Integer)
    ENV["OMP_NUM_THREADS"] = string(threads)
    ENV["PLUMED_NUM_THREADS"] = string(threads)
    set_num_threads(threads)
end

frames = Dict{String,System}()
loaded = Dict{String,Bool}()

function system_for(path::AbstractString)
    get!(frames, path) do
        System()
    end
end

function ensure_frame(path::AbstractString)
    sys = system_for(path)
    if !get(loaded, path, false)
        seams_system_read(sys, path, 1, -1)
        loaded[path] = true
    end
    return sys
end

for line in eachline(stdin)
    isempty(line) && continue
    startswith(line, "quit") && break
    threads_s, op, path = split(line)
    retune(parse(Int, threads_s))
    t0 = time_ns()
    if op == "read"
        sys = system_for(path)
        seams_system_read(sys, path, 1, -1)
        loaded[path] = true
        sink[] += seams_system_count(sys)
    elseif op == "knn"
        sys = ensure_frame(path)
        sink[] += seams_neighbours_knn(sys; k=4, cutoff=3.5, mutual=true)
    elseif op == "steinhardt"
        sys = ensure_frame(path)
        ql, bar = seams_steinhardt(sys, 6; cutoff=3.5, average=true)
        sink[] += ql[1] + bar[1]
    elseif op == "chill"
        sys = ensure_frame(path)
        labels, counts = seams_chill(sys; cutoff=3.5, type=0, graph=0, rule=0)
        sink[] += counts[2]
        sink[] += labels[1]
    elseif op == "rings"
        sys = ensure_frame(path)
        atoms = seams_rings_vertices(sys; cutoff=3.5, max_depth=6, index=0)
        sink[] += length(atoms)
    else
        println(stderr, "error ", op)
        exit(1)
    end
    ms = (time_ns() - t0) / 1e6
    println(stdout, ms)
    flush(stdout)
end
println(stderr, "sink ", sink[])
