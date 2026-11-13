# Science tutorials in Julia. Same numbers as science_tutorials.py.
# Usage: julia science_tutorials.jl ROOT [brine_t0 brine_tend]

using Libdl

const LIB = get(ENV, "SEAMS_LIB", "libseams")

function chill_file(path)
    sys = ccall((:seams_system_create, LIB), Ptr{Cvoid}, ())
    code = ccall((:seams_system_read, LIB), Cint, (Ptr{Cvoid}, Cstring, Cint, Cint), sys, path, 1, -1)
    code == 0 || error("read $path")
    n = ccall((:seams_system_count, LIB), Cint, (Ptr{Cvoid},), sys)
    labels = Vector{Cint}(undef, n)
    counts = Vector{Cint}(undef, 9)
    code = ccall((:seams_chill_plus, LIB), Cint,
                 (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
                 sys, 3.5, 0, 0, 0, labels, counts)
    ccall((:seams_system_destroy, LIB), Cvoid, (Ptr{Cvoid},), sys)
    code == 0 || error("chill+ $path")
    Int.(labels), Int.(counts)
end

function cages_of(path, spec)
    sys = ccall((:seams_system_create, LIB), Ptr{Cvoid}, ())
    ccall((:seams_system_read, LIB), Cint, (Ptr{Cvoid}, Cstring, Cint, Cint), sys, path, 1, -1) == 0 || error(path)
    n_cages = Ref{Csize_t}(0)
    n_atoms = Ref{Csize_t}(0)
    ccall((:seams_cages_signature, LIB), Cint,
          (Ptr{Cvoid}, Cdouble, Cint, Cint, Cstring, Ptr{Csize_t}, Ptr{Csize_t}),
          sys, 3.5, 0, 0, spec, n_cages, n_atoms) == 0 || error(spec)
    ccall((:seams_system_destroy, LIB), Cvoid, (Ptr{Cvoid},), sys)
    Int(n_cages[])
end

function read_dump(path)
    pos = NTuple{3,Float64}[]
    open(path) do io
        reading = false
        for line in eachline(io)
            if startswith(line, "ITEM: ATOMS")
                reading = true
            elseif reading
                parts = split(line)
                length(parts) < 5 && break
                push!(pos, (parse(Float64, parts[3]), parse(Float64, parts[4]), parse(Float64, parts[5])))
            end
        end
    end
    pos
end

function layer_period(projection)
    best_score, best = -1.0, 3.6
    mean = sum(projection) / length(projection)
    for step in 0:15
        period = 3.30 + (4.05 - 3.30) * step / 15
        re = im = 0.0
        for value in projection
            angle = 2π / period * (value - mean)
            re += cos(angle)
            im += sin(angle)
        end
        score = hypot(re, im) / length(projection)
        if score > best_score
            best_score = score
            best = period
        end
    end
    best
end

function layer_string(pos, labels, axis)
    chosen = [i for i in eachindex(labels) if labels[i] == 0 || labels[i] == 1]
    projection = [pos[i][axis] for i in chosen]
    kinds = [labels[i] for i in chosen]
    period = layer_period(projection)
    lo = minimum(projection)
    bins = [floor(Int, (value - lo) / period) for value in projection]
    letters = Char[]
    for index in 0:maximum(bins)
        selected = [kinds[j] for j in eachindex(bins) if bins[j] == index]
        isempty(selected) && continue
        hexagonal = count(==(1), selected)
        cubic = count(==(0), selected)
        hexagonal + cubic == 0 && continue
        push!(letters, hexagonal >= cubic ? 'H' : 'C')
    end
    String(letters)
end

function layer_score(projection, period)
    mean = sum(projection) / length(projection)
    re = im = 0.0
    for value in projection
        angle = 2π / period * (value - mean)
        re += cos(angle)
        im += sin(angle)
    end
    hypot(re, im) / length(projection)
end

root = ARGS[1]
ccall((:seams_set_num_threads, LIB), Cint, (Cint,), 1)
iface = joinpath(root, "input", "traj", "interface_f0.lammpstrj")
labels, counts = chill_file(iface)
stack = layer_string(read_dump(iface), labels, 3)
counts[2] == 352 && counts[3] == 434 && counts[4] == 49 && counts[1] == 1 || error(counts)
stack == "HHHHHHHHC" || error(stack)
println("interface hexagonal 352 water 434 interfacial 49 cubic 1")
println("interface layers ", stack)

function strongest(pos, labels)
    best_axis, best_score = 1, -1.0
    for axis in 1:3
        projection = [pos[i][axis] for i in eachindex(labels) if labels[i] == 0 || labels[i] == 1]
        score = layer_score(projection, layer_period(projection))
        if score > best_score
            best_score = score
            best_axis = axis
        end
    end
    best_axis
end

nuc = joinpath(root, "input", "traj", "nucleation_f1.lammpstrj")
nlabels, ncounts = chill_file(nuc)
npos = read_dump(nuc)
best_axis = strongest(npos, nlabels)
nstack = layer_string(npos, nlabels, best_axis)
ncounts[1] == 47 && ncounts[2] == 62 && ncounts[3] == 3589 && ncounts[4] == 240 || error(ncounts)
best_axis == 1 && nstack == "CHHHCHHHHHHHHH" || error((best_axis, nstack))
println("nucleus cubic 47 hexagonal 62 water 3589 interfacial 240")
println("nucleus axis x layers ", nstack)

sii = joinpath(root, "input", "traj", "genice_sII.lammpstrj")
sh = joinpath(root, "input", "traj", "genice_sH.lammpstrj")
cages_of(sii, "512") == 16 || error("sII 512")
cages_of(sii, "51264") == 8 || error("sII 51264")
cages_of(sh, "512") == 12 || error("sH 512")
cages_of(sh, "sh") == 8 || error("sH sh")
cages_of(sh, "51268") == 4 || error("sH 51268")
println("cages sII 512 16")
println("cages sII 51264 8")
println("cages sH 512 12")
println("cages sH sh 8")
println("cages sH 51268 4")

function hexagonal_oxygens(path)
    count = 0
    for line in eachline(path)
        startswith(line, "O ") || continue
        split(line)[5] == "hexagonal" && (count += 1)
    end
    count
end

if length(ARGS) < 3
    println("brine tutorial skipped: cannot fetch d-SEAMS/brine-front-230K at 4828e33ecf419c70643b2d466bfa4f283b4138dc")
else
    hexagonal_oxygens(ARGS[2]) == 733 || error("brine t0")
    hexagonal_oxygens(ARGS[3]) == 741 || error("brine tend")
    println("brine t0 hexagonal 733")
    println("brine tend hexagonal 741")
end
