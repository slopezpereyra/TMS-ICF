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

const S09_S1_BAD = [35, 37, 38, 40, 41, 43, 44, 46, 48, 50, 52, 54, 56, 57, 58, 59, 60,
    62:68, 70:85, 87:94, 96:192, 104:111, 113, 114,
    116:119, 121, 122, 124, 126:130, 132:134, 136:140,
    142, 144, 145, 147, 148, 149, 151, 152]

const S18_S1_GOOD = [81, 82, 84, 87, 89, 90, 93, 112,
    124, 129, 131, 134, 136, 137, 140, 142, 147, 149, 151]
const S18_S2_BAD = [38, 40, 42, 49, 51, 56, 59,
    61:66, 68:72, 74:78, 82, 83,
    85, 87, 88, 90, 93, 96, 98, 99,
    100, 103, 105, 107] #107 is last icf

const S49_S2_BAD = [55, 57, 116, 117, 120,
    121, 124, 125, 134, 135,
    136, 137, 152, 153, 154,
    161, 162, 170, 171]

const S53_S2_BADS = [50, 53, 56, 58, 61, #MEPS thus far
    129, 155, 156] # ICFs

const S76_BADS = [22, 23, 24, 25, 26, 27, 28, 31, 34, 38, 39, 40, 41]

const S80_S1_BADS = [103, 104, 118, 119, 120, 134, 137, 139]
const S80_S2_BADS = [36, 37, 103, 104, 38, 112, 113, 116, 118, 119, 120, 128, 136, 137, 152, 153]

const S87_S2_BADS = [44, 45, 46, 47, 95, 97, 98, 102, 103,
    106, 107, 126:129, 143, 168, 169, 172, 173, 181:184]

const S81_S1_BADS = [59, 60, 73, 74, 82, 83]
const S81_S2_BADS = [60, 61, 68, 110, 111]

const S88_S1_BADS = [102, 103, 107, 108, 126, 127, 130, 131, 171, 172]
const S88_S2_BADS = [76, 83, 84, 118, 169, 170]

const S90_S1_BADS = [52, 83, 84, 137, 138, 139] # 139 is Last ICF
const S90_S2_BADS = [37, 60, 61, 74, 85, 87, 90]

const S94_S1_BADS = [42, 58, 63, 64, 65, 66, 67, 71, 91]
const S95_S1_BADS = [34, 35, 36, 37, 53, 59:64, 68, 69, 114, 115] # 59 is ICF Start

const BADS = vcat(S49_S2_BAD, S53_S2_BADS, S76_BADS,
            S80_S1_BADS, S80_S2_BADS, S81_S1_BADS, S81_S2_BADS, 
            S87_S2_BADS, S88_S1_BADS, S88_S2_BADS, S90_S1_BADS, S90_S2_BADS, 
            S94_S1_BADS, S95_S1_BADS)

const PTINFO = CSV.read("TMS_ptinfo.csv", DataFrame)
const COLNAMES = ["Sample", "Subject", "Group", "Session", "Type", "Label", "ISI", "EMGPeakToPeak"]

const TOTAL_BADS = length(BADS) 

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
    a = a[:, setdiff(1:end, collect(2:size(a, 2)-1))]
    a[:, 1] = map(uppercase, a[:, 1])
    from = findfirst(x -> occursin("MEP", x) || occursin("EVP", x), a[:, 1])
    to = findlast(x -> occursin("END", x) || occursin("LAST", x), a[:, 1])
    a = a[from:to, :]
    a = a[:, setdiff(1:end, collect(2:size(a, 2)-1))]
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
    vcat(fill(missing, first_index - 1), isi_order)
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
        error("Non-numeric EMG ", f)
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
    session = parse(Int, file[11+5])
    subject = parse(Int, file[6+5:8+5])
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
    label_vector = get_label.(group_vector, type_vector)
    hcat(X[:, 1], subject_vector, group_vector, session_vector, type_vector, label_vector, X[:, 2:end])
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
function read_swipf(to_df::Bool=true)::Union{Matrix,DataFrame}
    fnames = readdir("Data/")
    X = Matrix{Real}(undef, 0, 8)
    count = 0
    for file in fnames
        fdata = read_swipf("Data/" * file)
        X = vcat(X, fdata)
        count += 1
    end
    print("COUNT ", count, "\n")
    return to_df ? DataFrame(X, :auto) : X
end

function get_label(group, type)::String 
    if group == 1 
        return type == "BL" ? "HC BL" : "HC SWD"
    end 
    if group == 2 
        return type == "BL" ? "MDD BL" : "MDD SWD"
    end
end

get_label.(X.Group, X.Type)

read_swipf("Data/SWIP_070_S2.txt")
X = read_swipf()

rename!(X, COLNAMES)
CSV.write("JL_DF.csv", X)
