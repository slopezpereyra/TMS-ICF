using CSV 
using DataFrames 
using DelimitedFiles 

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

const PTINFO = CSV.read("TMS_ptinfo.csv", DataFrame)
const COLNAMES = ["Sample", "Subject", "Group", "Session", "EMGPeakToPeak"]


"""
    read_swipf : String ↦ Matrix 

    Read a .txt tab-delimited file with SWIP data as exported from BrainSight,
    remove unwanted columns, and add categorical features.
    
    # Arguments
    - `f::String`: The name of the file to read
    
    # Examples
    ```julia-repl
    julia> read_swipf("Data/SWIP_052_S1.txt")

    93×7 Matrix{Any}:
     "MEP_FIRST"   70  1  2  "BL"    nothing   232.593
     "SAMPLE 15"   70  1  2  "BL"    nothing    95.7737
     "SAMPLE 16"   70  1  2  "BL"    nothing    72.9705
     "SAMPLE 17"   70  1  2  "BL"    nothing   201.809
     "SAMPLE 18"   70  1  2  "BL"    nothing    77.5311
     "SAMPLE 19"   70  1  2  "BL"    nothing   367.133
     "SAMPLE 20"   70  1  2  "BL"    nothing   140.24
     "SAMPLE 21"   70  1  2  "BL"    nothing   285.041
     "SAMPLE 22"   70  1  2  "BL"    nothing   205.229
     ⋮                                       ⋮         
     "SAMPLE 99"   70  1  2  "BL"  -1          339.769
     "SAMPLE 100"  70  1  2  "BL"  10         1014.75
     "SAMPLE 101"  70  1  2  "BL"  20          281.62
     "SAMPLE 102"  70  1  2  "BL"  -1          287.321
     "SAMPLE 103"  70  1  2  "BL"  15          358.011
     "SAMPLE 104"  70  1  2  "BL"   5           92.3532
     "SAMPLE 105"  70  1  2  "BL"  20          107.175
     "ICF_LAST"    70  1  2  "BL"   4           92.3532
"""
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
    a = add_isi_feature(a, f)
    add_categorical_features(a, f)
end

"""
    gen_isi_vector(first_index::Int, isi_order::Vector)

    Generates a vector of the form [nothing, …, first_isi, …, last_isi] to serve 
    as column vector in a SWIP data Matrix.

    # Arguments 
    `first_index::Int` : Index of the first ICF pulse in a SWIP file. 
    `isi_order::Vector` : One of the ISIORD constants.
    # Examples 
    julia> gen_isi_vector(3, ISIORDER72)
    74-element Vector{Union{Nothing, Int64}}:
       nothing
       nothing
     10
      8
     -1
      ⋮
      5
     20
      4
"""
function gen_isi_vector(first_index::Int, isi_order::Vector)
    vcat(fill(nothing, first_index - 1), isi_order)
end     

"""
    add_isi_feature(X::Matrix, f::String)::Matrix
    
    Given a SWIP data matrix X and a file from which that data originates, 
    add to it a column vector x⃑ such that xᵢis the ISI under which the ith ICF pulse 
    was elicited.

    # Arguments 
    `X::Matrix` : A matrix with SWIP data.
    `f::String` : A file from which the data in X originates.
"""
function add_isi_feature(X::Matrix, f::String)
    from = findfirst(x -> occursin("ICF", x), X[:, 1]) 
    to = findlast(x -> occursin("ICF", x), X[:, 1])
    icf_count = to - from + 1
    if icf_count != 120 && icf_count != 72
        print(f, "\n")
        @warn "WARNING: Erroneous ICF count on " f
    end
    if any(x -> !isa(x, Real), X[:, 2])
        error( "Non-numeric EMG ", f )
    end
    isi_order = icf_count == 120 ? ISIORD120 : ISIORDER72
    isi_vector = gen_isi_vector(from, isi_order)
                                    
    hcat(X[:, 1], isi_vector, X[:, 2:end])
end

"""
    add_categorical_features(X::Matrix, file::String)::Matrix
    
    Given a SWIP data matrix X and a file from which that data originates, 
    infer all categorical variables (Subject, Session Number, Session Type, 
    etc.), and add to X the column vectors corresponding to them.

    # Arguments 
    `X::Matrix` : A matrix with SWIP data.
    `file::String` : A file from which the data in X originates.
"""
function add_categorical_features(X::Matrix, file::String)::Matrix
    session = parse(Int, file[11+5]); subject = parse(Int, file[6+5:8+5]);
    session_vector = fill(session, size(X, 1))
    subject_vector = fill(subject, size(X, 1))
    sgroup = filter([:patientid] => x -> x == subject, PTINFO).group
    session_type::Vector = []
    if session == 1 
        session_type = filter(:patientid => x -> x == subject, PTINFO).session1
    else 
        session_type = filter(:patientid => x -> x == subject, PTINFO).session2
    end
    group_vector = fill(sgroup[1], size(X, 1))
    type_vector = fill(session_type[1], size(X, 1))
    hcat(X[:, 1], subject_vector, group_vector, session_vector, type_vector, X[:, 2:end])
end

"""
    read_swipf(to_df::Bool = true)::Union{Matrix, DataFrame}

    Read all SWIP files in the Data with swipf(f::file)::Matrix and 
    concatenate the matrix resulting from each of them into a single 
    large matrix. 

    This is the main data reading function. The resulting matrix is our fully 
    formatted and complete data. 

    # Arguments 
    `to_df::Bool` : Convert the resulting Matrix to DataFrame? Defaults to true.
"""
function read_swipf(to_df::Bool = true)::Union{Matrix, DataFrame}
    fnames = readdir("Data/")
    X = Matrix{Real}(undef, 0, 7)
    for file in fnames 
        fdata = read_swipf("Data/" * file)
        X = vcat(X, fdata)
    end
    return to_df ? DataFrame(X, :auto) : X
end

read_swipf("Data/SWIP_070_S2.txt")
X = read_swipf()

