using CSV 
using DataFrames 
using DelimitedFiles 
using "misc.jl"

# SWIP_052_S1 has EMG Peak to peak 1 and 2. Ask Jen.

const ISIORD120 = [
    10, 0, 8, 0, -1, 10, 0, -1, 20, 0, -1, 15, 0, 10,
    0, 20, 0, 8, 0, 10, 0, 15, 0, -1, 8, 0, -1, 10,
    0, -1, 5, 0, -1, 5, 0, 10, 0, 8, 0, 4, 0, 8, 0,
    -1, 4, 0, -1, 4, 0, 5, 0, -1, 20, 0, -1, -1, -1,
    5, 0, 5, 0, 10, 0, -1, -1, -1, 4, 0, -1, 20, 0, 8,
    0, 8, 0, 5, 0, -1, 15, 0, -1, 15, 0, -1, 5, 0, 15,
    0, -1, 15, 0, 15, 0, 4, 0, 8, 0, 20, 0, 4, 0, 20,
    0, 4, 0, -1, -1, 10, 0, 20, 0, -1, 15, 0, 5, 0, 20, 0, 4, 0
]

const ISIORDER72 = [
    10, 8, -1, 10, -1, 20, -1, 15, 10, 20, 8, 10, 15, -1, 8, -1, 10, -1,
    5, -1, 5, 10, 8, 4, 8, -1, 4, -1, 4, 5, -1, 20, -1, -1, -1, 5, 5,
    10, -1, -1, -1, 4, -1, 20, 8, 8, 5, -1, 15, -1, 15, -1, 5, 15, -1,
    15, 15, 4, 8, 20, 4, 20, 4, -1, -1, 10, 20, -1, 15, 5, 20, 4
]


function read_swipf(f::String)::Matrix
    a = readdlm(f, '\t')
    if size(a, 2) != 14
        @warn "Inconsistent feature count on" f
    end
    a = a[:, setdiff(1:end, collect(2:size(a, 2) - 1))]
    a[:, 1] = map(uppercase, a[:, 1])
    from = findfirst(x -> occursin("MEP", x) || occursin("EVP", x), a[:, 1])
    to = findlast(x -> occursin("END", x) || occursin("LAST", x), a[:, 1])
    a = a[from:to, :]
    a = a[:, setdiff(1:end, collect(2:size(a, 2) - 1))]
    validate(a, f)
    add_categorical_features(a, f)
end

function validate(X::Matrix, f::String)
    from = findfirst(x -> occursin("ICF", x), X[:, 1]) 
    to = findlast(x -> occursin("ICF", x), X[:, 1])
    icf_count = to - from + 1
    if icf_count != 120 
        @warn "WARNING: Erroneous ICF count on " f
    end
    if any(x -> !isa(x, Real), X[:, 2])
        error( "Non-numeric EMG" )
    end
end

function add_categorical_features(X::Matrix, file::String)::Matrix
    session_vector = fill(parse(Int, file[11 + 5]), size(X, 1))
    subject_vector = fill(parse(Int, file[6+5:8+5]), size(X, 1))
    hcat(X[:, 1], subject_vector, session_vector, X[:, 2:end])
end

function read_swipf()::Matrix 
    fnames = readdir("Data/")
    X = Matrix{Real}(undef, 0, 4)
    for file in fnames 
        print(file, "\n")
        fdata = read_swipf("Data/" * file)
        X = vcat(X, fdata)
    end
    return X
end

read_swipf()
read_swipf("Data/SWIP_076_S2.txt")

