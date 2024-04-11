module OPTinv

using DrWatson, Unitful, Revise, Statistics, BLUEs, LinearAlgebra,
    SparseArrays, TMI, Interpolations, UnitfulLinearAlgebra, Measurements,
    PyPlot, CSV, DataFrames#, TMItransient
import Interpolations.deduplicate_knots! as deduplicate! 
import Interpolations.LinearInterpolation as LI
import DelimitedFiles.readdlm
import Measurements.measurement, Measurements.value, Measurements.uncertainty
export readbaconfile, interpolatebacon, meanbacon, covariancebacon,
    formatbacon, covariancedims, steadystateM, firstguess,
    fillcovariance, fillcovariance!, linearleastsquares,
    solvesystem, core_locations, formattransientM, subsampletransientM, reshape,
    fancyscatter

yr = u"yr"
permil = Unitful.FixedUnits(u"permille")
K = u"K"

using DimensionalData
using DimensionalData: @dim, YDim, XDim, TimeDim, ZDim

@dim Cores "Cores" #can't use "Core" because that's the import for base Julia
@dim Modes "Modes"
@dim StateVar "state variable"

export yr, permil, K,
    Ti, Cores, Modes, StateVar,
    DiagRule, OffDiagRule

include(srcdir("UnitfulPlots.jl"))

"""
struct CovMat

holds covariance matrix information, including labelled axes

# Arguments
    - `ax`: axes as N-dim Tuple, e.g. [(T₁, M₁), (T₁, M₂), ...]
    - `mat`: Symmetric matrix, can be sparse + unitful, holds all info. 
"""
struct CovMat
    ax::Vector{Tuple} #maps on to axes 
    mat 
end

"""
struct Est

similar to BLUEs.DimEstimate, but y is stored as an array, not a vector
Cnn is CovMat type
Cnn.ax should be the same as the axes of vec(y)

# Arguments
    - `y`: DimArray of data
    - `Cnn`: associated CovMat 
"""
struct Est
    y::DimArray 
    Cnn::CovMat
end


"""
struct DiagRule

instructions for filling the diagonal of a matrix

# Arguments
    - `val`: value to fill diagonal with (could also be a function?)
    - `rule`: Tuple indicating what values to assign `val` to

# Example
    If target matrix has dimensions (Time × Mode × Variable), and we only want
    to fill diagonals with Mode = 11, we would set DiagRule.rule = (:, 11, :)
"""
struct DiagRule
    val::Quantity
    rule::Tuple
end

"""
struct OffDiagRule

instructions for filling the off-diagonal of a matrix
Note that, because CovMats are symmetric, x- and y- dims are arbitrary

# Arguments
    - `val`: value to fill diag.
    - `rule1`: Tuple indicating x-dims of off-diagonals 
    - `rule2`: Tuple indicating y-dims of off-diagonals 

# Example
    1. for a matrix with dimensions (Time × Mode × Variable), if we want to fill
    the off-diagonal indicated by (T=T, M=M, Variable1 = Variable2), we
    would supply (:, :, Var. 1) and (:, :, Var. 2) as our rules.
    This will fill the off-diagonal elements where time and mode are equivalent,
    but the state variables are different 
"""
struct OffDiagRule
    val::Quantity
    rule1::Tuple
    rule2::Tuple
end

"""
function readbacon(age_file::String, d18O_file::String)

read a bacon age model and associated data file in as a matrix with units

# Arguments
    - `age_file`: String of age file location - should be a .csv
    - `d18O_file`: String of d18O file lcoation - should be a .csv 
"""
function readbacon(age_file::String, d18O_file::String)
    ages = readdlm(age_file, ',', Float64, header = true, skipblanks = true)[1][:, 2:end]yr
    d18O = readdlm(d18O_file, ',', Float64, header = true)[1][:, 2]permil
    return ages, d18O
end


"""
function formatbacon(age_file::Vector{String}, d18O_file::Vector{String},
                       cores::Vector{Symbol}; res = 5yr)
go from a list of age_file, d18O_file to a formatted Est structure
will compute means first, and then figure out the common time axis
covariance matrices will only be computed over common time axis

# Arguments
    - `age_file`: Vector of paths to age_files
    - `d18O_file`: Vector of paths to d18O_files, should have same order as age_file
    - `cores`: Vector of core symbols, should have same order as d18O_file
# Optional Arguments
    - `res`: time resolution to interpolate to, default is 5 yr 
"""
function formatbacon(age_file::Vector{String}, d18O_file::Vector{String},
                       cores; res = 5yr)
    #compute the means first
    #because we want to know exactly what indices to compute Cnn for 
    function mean_shortcut(af, df)
        a, d = readbacon(af, df)
        T, int = interpolatebacon(a, d, res)
        μ, T = meanbacon(int, T, a)
        return μ, T, int, a
    end

    #compute all outputs 
    output = [mean_shortcut(f...) for f in zip(age_file, d18O_file)]
    #find common time axis
    T = intersect([o[2] for o in output]...)
    #target output dimensions
    dims = (Ti(T), Cores(cores))
    #trim all means to common time axis 
    mus = [mu[findall(x->x∈T, Ts_)] for (mu, Ts_) in output]
    #compute all Cnns over common time axis, keep as list 
    Cnn_vector = [covariancebacon(output[i][3], T, output[i][1], output[i][4], cores[i]) for i in 1:length(cores)]

    #make data into dimarry
    data = DimArray(hcat(mus...), dims)
    #generate CovMat.ax vector 
    cdims = covariancedims(dims)
    #safely combine Cnns according to generated CovMat.ax vector 
    Cnn = combine(Cnn_vector, vec(cdims))
    #return Est object 
    return Est(data, Cnn) 
    
end

"""
function combine(Cnn_vector::Vector{CovMat}, target_dims::Vector{Tuple})

combine a list of CovMats according to some specified axis

# Arguments
    - `Cnn_vector`: vector of CovMats
    - `target_dims`: target CovMat.ax
"""
function combine(Cnn_vector::Vector{CovMat}, target_dims::Vector{Tuple})
    N = length(target_dims)
    units = repeat([u for u in unitrange(Cnn_vector[1].mat)], length(Cnn_vector))
    #this potentially could be faster by using blockdiag and sparse
    #however, this is pretty flexible to non blockdiagonal covariance matrices
    mat = zeros(N, N)
    helper_da = DimArray(1:N, X(target_dims))
    for Cnn in Cnn_vector
        ind = helper_da[At(Cnn.ax)]
        mat[ind, ind] .= parent(Cnn.mat)
    end

    return CovMat(target_dims, UnitfulMatrix(Symmetric(mat), units, units .^ -1))
end

"""
function covariancedims(d::Tuple)

make a CovMat.dims axis according to a Tuple of dims

# Arguments
    - `d`: Tuple of dimensions, can be found by doing x.dims for a DimArray x 
"""
function covariancedims(d::Tuple)
    arr = Array{Tuple}(undef, length.(d)...)
    #assign tuples according to column-major vectorization 
    for ind in CartesianIndices(arr)
        arr[ind] = Tuple([d[i][ind[i]] for i in 1:length(ind)])
    end
    return arr
end

"""
function interpolatebacon(ages::Matrix, d18O::Vector,
                          res::Quantity)
for a matrix age model realizations (depth × realizations), generate a vector
of interpolation objects onto a common time axis for all realizations

# Arguments
    - `ages`: Matrix of age model realizations (depth × realizations), can be Unitful
    - `d18O`: Vector of d18O according to depth
    - `res`: resolution of target common time axis 

# Output
    - `T`: common time axis
    - `d`: vector of Interpolation objects 
"""
function interpolatebacon(ages::Matrix, d18O::Vector,
                          res::Quantity)
    utype = typeof(first(ages))
    numr = size(ages)[2]
    Tmin = minimum(ages)
    Tremain = Tmin % res
    Tfloor = Tremain ≥ 0.5res ? Tmin + (res - Tremain) : Tmin - Tremain
    T = Tfloor:res:floor(utype, maximum(ages)) 
    d = [LI(deduplicate!(reverse(ages[:, j]), move_knots = true), reverse(d18O)) for j in 1:numr]

    return T, d 
end

"""
function meanbacon(d::Vector, T, ages::Matrix)

for a vector of interpolation objects and a common time axis, compute the
mean d18O value at every common time axis point

if there are less than 10 realizations at a certain time value, don't provide a value

# Arguments
    - `d`: Vector of Interpolation objects
    - `T`: target common time axis
    - `ages`: matrix of age model realizations 
"""
function meanbacon(d::Vector, T, ages::Matrix)
    numr = size(ages)[2]
    μ = fill(NaN * permil, length(T))
    for (i, t) in enumerate(T)
        vals = [d[j](t) for j in 1:numr if ages[end, j] < t < ages[begin, j]]
        if length(vals) > 10 
            μ[i] = mean(vals)
        end
    end
        #remove part of μ, T where there are less than 2 realizations 
    hasdata = (!).(isnan.(μ))
    μ = μ[hasdata]
    T = T[hasdata]
    return μ, T 
end

"""
function covariancebacon(d, T, μ, ages, c)
compute the covariance matrix according to interpolation objects and common time axis

# Arguments
    - `d`: Vector of interpolation objects
    - `T`: target common time axis
    - `μ`: mean values from `meanbacon`
    - `ages`: Matrix of age model realizations
    - `c`: core name

# Output
    - CovMat object 
"""
function covariancebacon(d, T, μ, ages, c)
    Cnn = fill(NaN * permil^2, length(T), length(T))
    for (μ1, (i, t1)) in zip(μ, enumerate(T)) 
        Cnn[i, i:end] .= [covary(d, ages, t1, t2, μ1, μ2) for (μ2, t2) in zip(μ[i:end], T[i:end])]
    end
    axis = [(t,c) for t in T]
    units = sqrt.(unit.(diag(Cnn)))
    return CovMat(axis, UnitfulMatrix(Symmetric(ustrip.(Cnn)), units, units .^ -1))
end

"""
function fillcovariance(units::Vector, rules::Vector, udims::Tuple)

fill a covariance matrix according to a DiagRule or OffDiagRule

# Arguments
    - `units`: Vector of Units
    - `rules`: Vector of rules to fill matrix with
    - `udims`: Tuple of dimensions 
"""
function fillcovariance(units::Vector, rules::Vector, udims::Tuple)
    cdims = covariancedims(udims) #compute vector of covariance dims 
    N = length(cdims)
    #array of indices of vector indexed by dimensions
    #this may not be super efficient, but is very safe
    #allows us to take advantage of At() indexing in DimensionalData 
    index_array = DimArray(reshape(1:N, length.(udims)), udims)
    #init UnitfulMatrix Cnn - this allows us to be safe about setting Unitful vals. 
    Cnn = UnitfulMatrix(zeros(N, N), units, units .^-1)
    #iterate through all rules
    for r in rules
        #for a diagonal rule, find all indices that satisfy rule 
        if r isa DiagRule
            indices = index_array[r.rule...]
            [safesetindex!(Cnn, r.val, ind, ind) for ind in indices]
        #for an offdiagonal rule, just set one of the off-diagonals, Symmetric will reflect it 
        elseif r isa OffDiagRule
            indices1 = index_array[r.rule1...]
            indices2 = index_array[r.rule2...]
            for (ind1, ind2) in zip(indices1, indices2)
                safesetindex!(Cnn, r.val, ind1, ind2)
                safesetindex!(Cnn, r.val, ind2, ind1) 
            end
            
        end
            
    end
    return CovMat(vec(cdims), Cnn)
end

function fillcovariance!(c::CovMat, rules, udims)
    units = [u for u in unitrange(c.mat)]
    addme = fillcovariance(units, rules, udims)
    c = CovMat(c.ax, c.mat + addme.mat)
end

function fillcovariance!(e::Est, rules, udims)
    newC = fillcovariance!(e.Cnn, rules, udims) 
    e = Est(e.y, newC)
end

"""
function safesetindex!(mat::UnitfulMatrix, val::Quantity, i::Int, j::Int)

    version of setindex! from ULA that checks to make sure we have the right units
"""
function safesetindex!(mat::UnitfulMatrix, val::Quantity, i::Int, j::Int)
    targetunit = unitrange(mat)[i] / unitdomain(mat)[j] 
    if targetunit == unit.(val)
        setindex!(mat, ustrip.(val), i, j) 
    else
        error("units do not match")
    end
end


"""
function covary

    computes covariance across all interpolation objects 

    #Arguments 
    `d`: Vector of Interpolation objects for each realization
    `ages`: Matrix of ages, columns refer to realization
    `t₁`: time 1
    `t₂`: time 2 
    `μ₁`: mean value across all realizations at time 1 
    `μ₂`: mean value across all realizations at time 2
"""
function covary(d::Vector, ages::Matrix, t₁, t₂, μ₁, μ₂)
    vals  = [(d[k](t₁) - μ₁) * (d[k](t₂) - μ₂) for k in 1:length(d)
                 if ages[end, k] ≤ t₁ ≤ ages[begin, k] &&
                     ages[end, k] ≤ t₂ ≤ ages[begin, k]]
    if !isempty(vals) return sum(vals)/length(vals) else return 0permil^2 end 
end

"""
function formatbacon(μ, T, Cnn, c)

make an Est object, assumes you already computed mean and Cnn

# Arguments
    - `μ`: Vector of means
    - `T`: time
    - `Cnn`: CovMat
"""
function formatbacon(μ::Vector{Quantity}, T, Cnn::CovMat, c::Vector{Symbol})
    if length(size(μ)) !== 2
        μ = reshape(μ, (length(μ), 1))
    end
    if c isa Symbol
        c = [c]
    end
    
    da = DimArray(μ, (Ti(T), Cores(c)))
    return Est(da, Cnn) 
end

"""
function formatbacon(af::String, df::String, core::Symbol; res = 5yr)

start at age_file, d18O_file, go to Est object

# Arguments
    - `af`: age_file String
    - `df`: d18O_File String
    - `core`: Symbol, core name
# Optional Arguments
    - `res`: target time resolution, default = 5yr 
"""
function formatbacon(af::String, df::String, core::Symbol; res = 5yr)
    ages, d18O = readbacon(af, df) 
    T, int = interpolatebacon(ages, d18O, res);
    μ, T = meanbacon(int, T, ages)
    Cnn = covariancebacon(int, T, μ, ages, core)
    de = formatbacon(μ, T, Cnn, core)
end

"""
function steadystateM(locs, cores, TMIversion)

for some list of lat, lon boxes, compute the steady-state M matrix for core locs

# Arguments
    - `locs`: NamedTuple of core locations (lon, lat, depth) 
    - `patches`: NamedTuple of surface patch grid boxes 
"""
function steadystateM(locs::NamedTuple, patches::NamedTuple; TMIversion = "modern_180x90x33_GH11_GH12")
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    patches = NamedTuple{keys(patches)}([surfacepatch(patches[k]..., γ) for k in keys(patches)])
    si = [steadyinversion(Alu, p, γ) for p in patches]
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]
    obs = [observe(s, wis, γ) for s in si]
    #θ = readfield(TMIfile, "θ", γ) #potential temp
    #d18O = readfield(TMIfile, "δ¹⁸Ow", γ) #.- 0.27
    modes = [m for m in keys(patches)]
    cores = [c for c in keys(locs)]
    return DimArray(hcat(obs...)', (Modes(modes), Cores(cores)))
end

"""
function firstguess(T, modes, σθ, σδ, c; fill_val = [0K, 0permil])

    compute u₀ and Cᵤᵤ
    The off-diagonals of Cᵤᵤ will be populated by either
    σθ^2/c or σδ^2c
    If c < σθ/σδ, then the off-diagonals must be populated by σδ^2c
    Otherwise, populated by σθ^2/c
"""
function firstguess(T, modes, σθ, σδ, ρ; fill_val = [0K, 0permil])
    z = zeros(length(T), length(modes))
    dims = (Ti(T), Modes(modes), StateVar([:θ, :δ]))
    da = DimArray(cat(z * K .+ fill_val[1], z * permil .+ fill_val[2], dims = 3), dims)
    units = unit.(vec(da))
    diagrules = [DiagRule(σθ^2, (:, :, At(:θ))), DiagRule(σδ^2, (:, :, At(:δ)))]
    cov = ρ*σθ*σδ
    @show cov
    @show c = ρ*σθ/σδ
    rules = vcat(diagrules, OffDiagRule(cov, (:, :, At(:θ)), (:, :, At(:δ))))
    Cuu = fillcovariance(vec(units), rules, dims)
    @show isposdef(Cuu.mat)

    return Est(da, Cuu) 
end

function solvesystem(e::Est, u₀::Est, E::UnitfulMatrix, predict::Function)
    yanom = e.y .- mean(e.y, dims = Ti)
    yvals = vec(yanom)
    y = UnitfulMatrix(ustrip(yvals)[:], unit.(yvals))
    up = UnderdeterminedProblem(y, E, e.Cnn.mat, u₀.Cnn.mat, u₀.y)
    ũ = solve(up)
    u₀de = DimEstimate(vec(u₀.y), u₀.Cnn.mat, u₀.y.dims)
    ỹ = predict(value.(ũ.x))
    Cññ = E * ũ.C * E'
    yde = DimEstimate(yvals, parent(e.Cnn.mat) .* permil ^2 , e.y.dims)
    ỹde = DimEstimate(E*ũ.v, parent(Cññ) .* permil^2, e.y.dims)
    return yde, ỹde, u₀de, ũ
end


"""
function linearleastsquares(x,y)

    use linear algebra representation to find the linear least squares solution of
    y = a*x+b

    returns lls: vector where a = lls[1], b = lls[2]
"""
function linearleastsquares(x,y)
    E = ones(length(x), 2)
    E[:, 1] .= x
    F = (E'*E)\I*E'
    lls = F*y
    return lls
end

function  formattransientM(arr::Array, τ, modes, cores)
    ℳ = DimArray(arr, (Ti(τ), Modes(modes), Cores(cores)))
    return ℳ
end

function subsampletransientM(ℳ::DimArray, newres::Quantity)
    dims = ℳ.dims
    τ = [t for t in dims[1]]
    sum_ind = [a:1yr:a+newres-1yr for a in 1yr:newres:length(τ)yr-newres+1yr]
    subℳ = cat([sum(ℳ[Ti = At(ind)], dims = Ti) for ind in sum_ind]..., dims = Ti)

    #unfortunately its really hard to make 
    Mbegin = Array(reshape(ℳ[At(0yr), :, :], (1, size(ℳ)[2:3]...)))
    arr = cat(Mbegin, Array(subℳ), dims = 1)
    newτ = 0yr:newres:1000yr
    ℳnew = DimArray(arr, (Ti(newτ), Modes([m for m in dims[2]]), Cores([c for c in dims[3]])))
    return ℳnew, newτ
end


function core_locations()
    locations_file = "../data/EN539_locations.csv"
    cores_df = CSV.read(locations_file,DataFrame)
    locs = [(cores_df.Lon[i], cores_df.Lat[i], cores_df.Depth[i]) for i in 1:length(cores_df.Core)]
    cores = replace.(cores_df.Core, "EN539-"=>"")
    return NamedTuple{Tuple(Symbol.(cores))}(locs)
    
end


end
