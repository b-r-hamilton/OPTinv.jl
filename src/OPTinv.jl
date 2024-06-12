module OPTinv

using DrWatson, Unitful, Revise, Statistics, BLUEs, LinearAlgebra,
    SparseArrays, TMI, Interpolations, UnitfulLinearAlgebra, Measurements,
    PyPlot, CSV, DataFrames, JLD2, NCDatasets, NaNMath, XLSX, Downloads,
    Dates
import Interpolations.deduplicate_knots! as deduplicate! 
import Interpolations.LinearInterpolation as LI
import DelimitedFiles.readdlm
import Measurements.measurement, Measurements.value, Measurements.uncertainty

export readbacon, interpolatebacon, meanbacon, covariancebacon,
    formatbacon, covariancedims, steadystateM, firstguess,
    fillcovariance, fillcovariance!, linearleastsquares,
    solvesystem, core_locations, formattransientM, subsampletransientM,
    γreshape, invert, inversion, solution

export loadM, loadcores, loadE, ptobserve, fancytext, boxmean, reconstructsurface,
    monthlyannual, γbox, geospatialsubset

export loadLMR, loadOcean2k, loadThornalley, loadEN4, loadCESMLME

export record, loadEN4MLD 

export orthographic_axes

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

struct inversion
    cores::Vector{Symbol}
    modetype::String
    name::String
    color::String
end

struct solution
    y::DimEstimate
    ỹ::DimEstimate
    u₀::DimEstimate
    ũ::DimEstimate
    θ::Vector
    δ::Vector
    γ::TMI.Grid
    spatialmodes::Matrix
    name::String
    color::String
end


function invert(x::inversion)
    corenums_full = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
    corenums_sorted = [i for i in corenums_full if i in x.cores]
    filename = x.modetype * ".jld2"
    ℳ, spatialmodes = loadM(filename, corenums_sorted)
    y = loadcores(corenums_sorted)
    ρ = 0.99
    T = Array(y.y.dims[1])
    res = unique(diff(T))[1]
    Tᵤ = 500.0yr:res:T[end]
    jld = jldopen(datadir("modemags.jld2"))
    mags = jld[x.modetype * "mags"]
    
    σθ = vec(std(mags, dims = 1)) * K

    if x.modetype == "svd"    
        σθ *= 5
    else
        σθ *= 0.1
    end
    
    σδ = σθ ./ 10 .* permil/K

    u₀ = firstguess(Tᵤ, ℳ.dims[2][:], σθ, σδ, ρ)
    E, predict = loadE(filename, ℳ, Tᵤ, T, σθ, σδ, ρ, corenums_sorted, y.Cnn.ax)
    yde, ỹde, u₀de, ũ = solvesystem(y, u₀, E, predict)
    θ, δ = reconstructsurface(spatialmodes, ũ)
    
    TMIversion = "modern_180x90x33_GH11_GH12"
    γ = Grid(download_ncfile(TMIversion))

    return solution(yde, ỹde, u₀de, ũ, θ, δ, γ,spatialmodes, x.name, x.color) 
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
        return means(μ, T, int, a)
    end

    #compute all outputs 
    ml = NamedTuple{Tuple(cores)}([mean_shortcut(f...) for f in zip(age_file, d18O_file)])
    #find common time axis
    T = intersect([m.T for m in ml]...)
    #target output dimensions
    dims = (Ti(T), Cores(cores))
    #trim all means to common time axis 
    mus = NamedTuple{Tuple(cores)}([ml[c].μ[findall(x->x∈T, ml[c].T)] for c in cores])
    #compute all Cnns over common time axis, keep as list 
    Cnn_vector = [covariancebacon(ml[c].int, T, mus[c], ml[c].ages, c) for c in cores]
    Cnns = NamedTuple{Tuple(cores)}(Cnn_vector) 
    #make data into dimarry
    data = DimArray(hcat([mus[c] for c in cores]...), dims)
    #generate CovMat.ax vector 
    cdims = covariancedims(dims)
    #safely combine Cnns according to generated CovMat.ax vector 
    Cnn = combine(Cnns, vec(cdims))
    #return Est object 
    return Est(data, Cnn) 
    
end

struct means
    μ
    T
    int
    ages
end


"""
function combine(Cnn_vector::Vector{CovMat}, target_dims::Vector{Tuple})

combine a list of CovMats according to some specified axis

# Arguments
    - `Cnn_vector`: vector of CovMats
    - `target_dims`: target CovMat.ax
"""
function combine(Cnns::NamedTuple, target_dims::Vector{Tuple})
    N = length(target_dims)
    units = repeat(Array(unitrange(Cnns[1].mat)), length(Cnns))
    #this potentially could be faster by using blockdiag and sparse
    #however, this is pretty flexible to non blockdiagonal covariance matrices
    mat = zeros(N, N)
    helper_da = DimArray(1:N, X(target_dims))
    for Cnn in Cnns
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
        if length(vals) > 100
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
    return CovMat(axis, UnitfulMatrix(Matrix(Symmetric(ustrip.(Cnn))), units, units .^ -1))
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
        
        error("unit = " * string(unit(val)) * ", should be " * string(targetunit))
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
function formatbacon(μ::Vector, T, Cnn::CovMat, c::Symbol)
    if length(size(μ)) !== 2
        μ = reshape(μ, (length(μ), 1))
    end
    if c isa Symbol
        c = [c]
    end
    
    da = DimArray(μ, (Ti(T), Cores(vec(c))))
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
function firstguess(T, modes, σθ::Quantity, σδ::Quantity, ρ; fill_val = [0K, 0permil])
    z = zeros(length(T), length(modes))
    dims = (Ti(T), Modes(modes), StateVar([:θ, :δ]))
    da = DimArray(cat(z * K .+ fill_val[1], z * permil .+ fill_val[2], dims = 3), dims)
    units = unit.(vec(da))
    diagrules = [DiagRule(σθ^2, (:, :, At(:θ))), DiagRule(σδ^2, (:, :, At(:δ)))]
    cov = ρ*σθ*σδ
    rules = vcat(diagrules,
                 OffDiagRule(cov, (:, :, At(:θ)), (:, :, At(:δ))))
    Cuu = fillcovariance(vec(units), rules, dims)
    @show isposdef(Cuu.mat)

    return Est(da, Cuu) 
end

function firstguess(T, modes, σθ::Vector, σδ::Vector, ρ; fill_val = [0K, 0permil]) 
    z = zeros(length(T), length(modes))
    dims = (Ti(T), modes, StateVar([:θ, :δ]))
    da = DimArray(cat(z * K .+ fill_val[1], z * permil .+ fill_val[2], dims = 3), dims)
    units = unit.(vec(da))
    N = length(modes)
    diagrules = vcat([DiagRule(σθ[i]^2, (:, At(i), At(:θ))) for i in modes], [DiagRule(σδ[i]^2, (:, At(i), At(:δ))) for i in modes])

    offdiagrules = [OffDiagRule(ρ*σθ[i]*σδ[i], (:, At(i), At(:θ)), (:, At(i), At(:δ))) for i in modes]
    Cuu = fillcovariance(vec(units), vcat(diagrules, offdiagrules), dims)
    @show isposdef(Cuu.mat) 
    return Est(da, Cuu) 
end

function solvesystem(e::Est, u₀::Est, E::UnitfulMatrix, predict::Function)
    if e.y.dims[1][begin] < u₀.y.dims[1][begin]
        error("Tᵤ must precede T")
    end
    
    yanom = e.y .- mean(e.y[DimensionalData.Between(1850.0yr, 1980.0yr), :], dims = Ti)
    #yanom = e.y .- mean(e.y, dims = Ti)
    yvals = vec(yanom)
    y = UnitfulMatrix(ustrip(yvals)[:], unit.(yvals))
    up = UnderdeterminedProblem(y, E, e.Cnn.mat, u₀.Cnn.mat, u₀.y)
    ũ = solve(up)
    u₀de = DimEstimate(vec(u₀.y), u₀.Cnn.mat, u₀.y.dims)
    ỹ = predict(value.(ũ.x))
    Cññ = E * ũ.C * E'
    yde = DimEstimate(yvals, parent(e.Cnn.mat) .* permil ^2 , e.y.dims)
    ỹde = DimEstimate(E*ũ.v, parent(Cññ) .* permil^2, e.y.dims)

    
    println("Cost fncts") 
    @show cost(ũ, up)
    @show datacost(ũ, up)
    @show controlcost(ũ, up)
    @show BLUEs.rmserror(ũ, up)
    @show BLUEs.rmscontrol(ũ, up, :θ)
    @show BLUEs.rmscontrol(ũ, up, :δ)

    return yde, ỹde, u₀de, ũ
end


"""
function linearleastsquares(x,y)

    use linear algebra representation to find the linear least squares solution of
    y = a*x+b

    returns lls: vector where a = lls[1], b = lls[2]
"""
function linearleastsquares(x,y;C=nothing)
    E = ones(length(x), 2)
    E[:, 1] .= x
    F = (E'*E)\I*E'
    lls = F*y
    if !isnothing(C)
        Cnew = F*C*F'
        return lls, Cnew
    else
        return lls
    end
end

"""
doesnt work rn 
"""
function linearleastsquares(x::UnitfulMatrix, y::UnitfulMatrix)
    E = ones(length(x), 2)
    E[:, 1] .= parent(x)
    E = UnitfulMatrix(E, unitrange(x), [unit(y[1])/unit(x[1]), unit(y[1])])
    F = (E'*E)\UnitfulMatrix(I)*E'
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
    sum_ind = [a:1yr:a+newres-1yr for a in 1yr:newres:length(τ)yr-newres]
    subℳ = cat([sum(ℳ[Ti = At(ind)], dims = Ti) for ind in sum_ind]..., dims = Ti)

    #unfortunately its really hard to make 
    Mbegin = Array(reshape(ℳ[At(0yr), :, :], (1, size(ℳ)[2:3]...)))
    arr = cat(Mbegin, Array(subℳ), dims = 1)
    newτ = 0yr:newres:(size(arr)[1]-1)*newres
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


function γreshape(v::Vector{T}, γ) where T
    template = Array{T}(undef, size.(γ.axes)[1][1], size.(γ.axes)[2][1])
    template .= NaN
    template[γ.wet[:, :, 1]] .= v
    return template
end

function loadM(filename::String, core_list::Vector{Symbol};  res = 10yr) 
    filepath = joinpath("../data/M", filename)

    #load in variables from ex3.svdmodes.jl 
    jld = jldopen(filepath)
    M = jld["arr"][begin:end-1, :, :]
    τ = jld["τ"][begin:end-1]yr
    spatialmodes = jld["modes"]
    modes = 1:size(M)[2]
    cores = keys(core_locations())
    close(jld)

    #format into DimArray 
    ℳ = formattransientM(M, τ, [m for m in modes], [c for c in cores])

    
    #normalize modes and ℳ so that the surface mode sums to 1
    norm = sum(spatialmodes, dims = 2)
    [spatialmodes[i, :] ./= norm[i] for i in 1:11]
    [ℳ[:, At(i), :] ./= norm[i] for i in 1:11]

    
    #M has 1 yr resolution
    ℳ, τ = subsampletransientM(ℳ, res)
    return ℳ[:, :, At(core_list)], spatialmodes
    #return ℳ, spatialmodes
end

function loadcores(core_list::Vector{Symbol}, res = 10yr)
    files = readdir(datadir())
    ae_files = [f for f in files if occursin("ae", f)]
    d18O_files = [f for f in files if occursin("d18O", f)]
    directory_cores = Tuple([Symbol(split(f, ".")[1]) for f in ae_files])

    #re-order directory cores so they line up with core_list
    #@show sort_indices = [findall(x->x== c, core_list)[1] for c in directory_cores if c in core_list]
    sort_indices = []
    for c in core_list
        ind = findall(x->x == c, directory_cores)
        if ! isempty(ind) 
            push!(sort_indices, (ind[1]))
        end
    end
    #=
    @show ae_files[sort_indices]
    @show d18O_files[sort_indices]
    @show core_list
    =#
    println("Pulling in core data, calculating covariance matrices") 
    @time e = formatbacon(datadir.(ae_files)[sort_indices], datadir.(d18O_files)[sort_indices], core_list, res = res)
    #add in measurement uncertainty 
    e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))
    @show isposdef(e.Cnn.mat)
    return e 
end

function loadE(filename::String, ℳ::DimArray, Tᵤtarget, Ttarget, u₀::Est,
               core_list::Vector{Symbol}, yax::Vector{Tuple})
    sv = (:θ, :δ)
    coeffs = NamedTuple{sv}([-0.224permil/K, 1])
    τ = [t for t in ℳ.dims[1]]

    #T = collect(0yr:10yr:2020yr)
    #Tᵤ = -500yr:10yr:2020yr
    
    #this is the same as BLUEs.convolve, I just did it explicitly here 
    predict(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ[t .- τ .> minimum(Tᵤ)]]...) for s in sv]...) for t in T]..., dims = 2)', (Ti(T), Cores(core_list)))

    filename = split(filename, ".")[1] .* "_E.jld2"
    filepath = joinpath("../data/M", filename) 
    if !isfile(filepath)
        println("Computing E, will take a little while!") 
        @time E = impulseresponse(predict, u₀.y) #328 seconds
        jldsave(filepath; E)
    else
    
        println("Loading in pre-saved file") 
        jld = jldopen(filepath)
        E = jld["E"]
        close(jld)
    end
    
    return E, predict
    #following subsetting doesn't work but might be moot because we need to calc. over a new time period anyways... 

    #=
    allcores = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9, 13, 14]) .* "A")
    
    Erange = vec(covariancedims((Ti(T), Cores(allcores))))
    Edomain = u₀.Cnn.ax 
    xindices = [x[1] ∈ Ttarget && x[2]∈core_list for x in Erange]
    yindices = [y[1] ∈ Tᵤtarget for y in Edomain] 
    return UnitfulMatrix(parent(E[indices, :]), unitrange(E)[xindices], unitdomain(E)[yindices]), predict
   =# 
end

#alternative version, dont provide firstguess
#will calculate over longer time period 
function loadE(filename::String, ℳ::DimArray, Tᵤtarget, Ttarget,
               σθ::Vector, σδ::Vector, ρ,
               core_list::Vector{Symbol}, yax::Vector{Tuple})
    sv = (:θ, :δ)
    coeffs = NamedTuple{sv}([-0.224permil/K, 1])
    τ = [t for t in ℳ.dims[1]]

    #maximum time period we're interested in for this problem 
    T = collect(0yr:10yr:2020yr)
    Tᵤ = -500yr:10yr:2020yr
    u₀ = firstguess(Tᵤ, ℳ.dims[2][:], σθ, σδ, ρ) 
    
    #this is the same as BLUEs.convolve, I just did it explicitly here 
    predict(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ[t .- τ .> minimum(Tᵤ)]]...) for s in sv]...) for t in T]..., dims = 2)', (Ti(T), Cores(core_list)))

    filename = split(filename, ".")[1] .* "_E.jld2"
    filepath = joinpath("../data/M", filename) 
    if !isfile(filepath)
        println("Computing E, will take a little while!") 
        @time E = impulseresponse(predict, u₀.y) 
        jldsave(filepath; E)
    else
    
        println("Loading in pre-saved file") 
        jld = jldopen(filepath)
        E = jld["E"]
        close(jld)
    end
    
    allcores = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9, 13, 14]) .* "A")
    Erange = vec(covariancedims((Ti(T), Cores(allcores))))
    modes = ℳ.dims[2]
    Edomain = vec(covariancedims((Ti(Tᵤ), Modes(1:11), StateVar([s for s in sv]))))
    xindices = [x[1] ∈ Ttarget && x[2]∈core_list for x in Erange]
    
    yindices = [y[1] ∈ Tᵤtarget && y[2] ∈ modes for y in Edomain]
    
    predict_subset(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ[t .- τ .> minimum(Tᵤtarget)]]...) for s in sv]...) for t in Ttarget]..., dims = 2)', (Ti(Ttarget), Cores(core_list)))
    return UnitfulMatrix(parent(E[xindices, yindices]), fill(permil, sum(xindices)), unitdomain(E)[yindices]), predict_subset
 
end

function ptobserve(pt::Tuple, γ, spatialmodes::Matrix, ũ)
    surface = convert.(Int64, zeros(180,90))
    surface[γ.wet[:, :, 1]] .= 1:11113

    γind = [findmin(abs.(γl .- p))[2] for (γl, p) in zip([γ.lon, γ.lat], pt)]
    surfind = surface[γind...]#ind in 1:11113
    #Spatial Mode mags. assoc. with this pt 
    sm_pt = DimArray(spatialmodes[:, surfind], (ũ.dims[2])) 

    #we need to construct a matrix O that can multiply ũ
    #to generate an obs. ts at a pt
    #get covariance dims associated with this matrix (Tᵤ x Sv) 
    x = vec(covariancedims((ũ.dims[1], ũ.dims[3])))
    y = vec(covariancedims(ũ.dims))
    O = zeros(length(x), length(y))
    
    #now, at the correct indices (T = T, Sv = Sv), put the mode value 
    for (i, xi) in enumerate(x)
        for (j, yj) in enumerate(y)
            if xi[1] == yj[1] 
                if xi[2] == yj[3]
                    O[i,j] = sm_pt[yj[2]]
                end
            end
        end
    end

    #compute Oũ (Estimate), then make a DimEstimate
    Tl = length(ũ.dims[1])
    Oũ = UnitfulMatrix(O, vcat(fill(K, Tl), fill(permil, Tl)), unitrange(ũ.v)) * ũ
    return DimEstimate(Oũ.v, Oũ.C, (ũ.dims[1], ũ.dims[3]))
end

#assumes fld is (lon x lat x time) 
function boxmean(fld::Array{T, N}, lat, lon, latsouth, latnorth, lonwest, loneast) where {T, N}
    latnorth = findmin(abs.(latnorth .- lat))[2]
    latsouth = findmin(abs.(latsouth .- lat))[2]
    lonwest = findmin(abs.(lonwest .- lon))[2]
    loneast = findmin(abs.(loneast .- lon))[2]
    fldμbox = Vector{T}(undef, size(fld)[end])
    box = Vector{Int64}
    if lonwest > loneast 
        west = fld[lonwest:end, latsouth:latnorth, :]
        east = fld[begin:loneast, latsouth:latnorth, :]
        box = cat(west, east, dims = 1)
        #println("concatenating") 
    else
        box = fld[lonwest:loneast, latsouth:latnorth, :]
    end
    if length(size(fld)) > 2 
        [fldμbox[i] = mean(box[:, :, i][(!).(isnan.(box[:, :, i]))]) for i in 1:size(fld)[3]]
        return fldμbox
    else
        return NaNMath.mean(convert(Vector{Float64}, box[(!).(ismissing.(box))]))
    end
end

"""
function γbox(γ::TMI.Grid, latsouth, latnorth, lonwest, loneast)

    find indices of surfacevector corresponding to a box 
"""
function γbox(γ::TMI.Grid, latsouth, latnorth, lonwest, loneast)
    surfind = findall(x->x == 1, γ.wet[:, :, 1])
    lon_index = [i[1] for i in surfind]
    lat_index = [i[2] for i in surfind]
    lat = γ.lat
    lon = γ.lon 
    latnorth = findmin(abs.(latnorth .- lat))[2]
    latsouth = findmin(abs.(latsouth .- lat))[2]
    lonwest = findmin(abs.(lonwest .- lon))[2]
    loneast = findmin(abs.(loneast .- lon))[2]

    if loneast > lonwest
        return intersect(findall(x->lonwest<x<loneast, lon_index), findall(x->latsouth<x<latnorth, lat_index))
    elseif lonwest > loneast #this means that we wrap around the centerline
        return intersect(vcat(findall(x->lonwest<x, lon_index),
                         findall(x->x<loneast, lon_index)),
                         findall(x->latsouth<x<latnorth, lat_index))
    end
end

#calculates full covariance matrix spatially, very expensive 
#=
function reconstructsurface(spatialmodes, ũ)
    
    θ = Vector{Estimate}(undef, length(ũ.dims[1]))
    δ = Vector{Estimate}(undef, length(ũ.dims[1]))
    covdims = vec(covariancedims(ũ.dims))
    convertθ = UnitfulMatrix(Matrix{Float32}(undef, 11113, length(covdims)), fill(K, 11113), unit.(vec(ũ.x)));
    convertδ = UnitfulMatrix(Matrix{Float32}(undef, 11113, length(covdims)), fill(permil, 11113), unit.(vec(ũ.x)));
    for (i, t) in enumerate(ũ.dims[1])
        @show t 
        indθ = findall(x->x[1] == t && x[3] == :θ, covdims)
        indδ = findall(x->x[1] == t && x[3] == :δ, covdims)
        println("found indices") 
        [convertθ[:, i] = spatialmodes[covdims[i][2], :] for i in indθ]
        [convertδ[:, i] = spatialmodes[covdims[i][2], :] for i in indδ]
        println("filled mat") 
        θ[i] = convertθ * ũ
        println("first conversion")
        δ[i] = convertδ * ũ
    end
    return θ, δ
end

=#

function reconstructsurface(spatialmodes, ũ)
    θvals = [spatialmodes' * vec(ũ.x[At(t), :, At(:θ)]) for t in ũ.dims[1]]
    δvals = [spatialmodes' * vec(ũ.x[At(t), :, At(:δ)]) for t in ũ.dims[1]]
    return θvals, δvals
end

"""
function load_EN4(analysis_path, pt)

    assumes you have downloaded all EN4 analyses locally
    this is a bit of a pain, used a wget script
    https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html

    # Arguments
    - analysis_path: where the objective analyses are located
    - lon: approximate lon pt to observe
    - lat
"""
function loadEN4(analysis_path, lonpt, latpt)
    files = readdir(analysis_path)
    nc = NCDataset(joinpath(analysis_path,files[1]))
    lon = nc["lon"][:]
    lat = nc["lat"][:]
    
    @show lonindex = [findmin(abs.(lon .- lp))[2] for lp in lonpt] 
    @show latindex = [findmin(abs.(lat .- lp))[2] for lp in latpt] 
    close(nc)
    time = Vector{DateTime}(undef, length(files))
    temp = Vector{Float64}(undef, length(files)); salinity = copy(temp)
    for (i,f) in enumerate(files)
        nc = NCDataset(joinpath(analysis_path, files[i]))
        time[i] = nc["time"][1]
        if length(lonindex) == 1 
            temp[i] =  nc["temperature"][lonindex[1], latindex[1], 1, 1]
            salinity[i] = nc["salinity"][lonindex[1], latindex[1], 1, 1]
        elseif length(lonindex) == 2 
            temp[i] = boxmean(nc["temperature"][:, :, 1, 1], lat, lon, minimum(latindex), maximum(latindex), lonpt[1], lonpt[2])
            salinity[i] = boxmean(nc["salinity"][:, :, 1, 1], lat, lon, minimum(latindex), maximum(latindex), lonpt[1], lonpt[2])
        end
        close(nc)
    end
    return time, temp, salinity 
end


"""
    function winterannual(x::Vector{T}, time::Vector{DateTime})

"""
function monthlyannual(x::Vector{T}, time::Vector{DateTime}, months = 1:12) where T
    months = collect(1:12)
    data = Vector{T}(undef, length(unique(year.(time))))
    for (i, y) in enumerate(unique(year.(time)))
        ind = intersect(findall(x->x == y, year.(time)), findall(x->x ∈ months, month.(time)))
        data[i] = NaNMath.mean(x[ind]) 
    end
    return data 
end

function orthographic_axes(latsouth,latnorth, lonwest, loneast, inc)
    southx = lonwest:inc:loneast-inc
    southy = latsouth * ones(length(southx))
    easty = latsouth:inc:latnorth-inc
    eastx = loneast * ones(length(easty))
    northx = loneast:-inc:lonwest+inc
    northy = latnorth * ones(length(northx))
    westy = latnorth:-inc:latsouth
    westx = lonwest * ones(length(westy))
    return vcat(southx, eastx, northx, westx), vcat(southy, easty, northy, westy)
end

function geospatialsubset(mat, lat, lon, lat_target, lon_target)
    lat_ind = [findmin(abs.(lat .- lt))[2] for lt in lat_target]
    lon_ind = [findmin(abs.(lon .- lt))[2] for lt in lon_target]
    if lon_ind[1] < lon_ind[2]
        mat[begin:lon_ind[1], :] .= 0
        mat[lon_ind[1]:lon_ind[2], begin:lat_ind[1]] .= 0
        mat[lon_ind[1]:lon_ind[2], lat_ind[2]:end] .= 0 
        mat[lon_ind[2]:end, :] .= 0
        return mat 
    else
        mat[lon_ind[1]:end, begin:lat_ind[1]] .= 0
        mat[lon_ind[1]:end, lat_ind[2]:end] .= 0
        mat[begin:lon_ind[2],  begin:lat_ind[1]] .= 0
        mat[begin:lon_ind[2],  lat_ind[2]:end] .= 0
        mat[lon_ind[2]:lon_ind[1], :] .= 0
        return mat 
        
    end
end

function loadCESMLME(d, start, stop)
    f = readdir(d)
    startyear = [parse(Int64, f_[end-15:end-12]) for f_ in f]
    endyear = [parse(Int64, f_[end-8:end-5]) for f_ in f]
    #find first file that the start is greater than the start index 
    ind = findall(x->x>start, startyear)[1]-1

    nc = NCDataset(joinpath(d, f[ind]))
    time = year.(nc["time"][:])
    tind = findall(x->start < x < stop, time)
    
    sst = mean(nc["SST"][:, :, 1, tind], dims = 3)[:, :, 1]

    lon = nc["TLONG"][:, :]
    lat = nc["TLAT"][:, :]
    return sst, lon, lat 
    
end

"""
struct record

    structure for dealing with EN4 CTD profiles 
"""
struct record
    lat::Number
    lon::Number
    depth::Vector
    t::DateTime
    T::Vector
    S::Vector
    p::Vector 
    CT::Vector
    SA::Vector
    ρ::Vector
    MLDT::Number #temperature value at MLD-temperature 
    MLDS::Number #salinity value at MLD-salinity
end

function loadEN4MLD()
    jld = jldopen(datadir("EN4_march_mld.jld2"))
    mldtemp = jld["mldtemp"]
    mldsalinity = jld["mldsalinity"]
    mlddepthtemp = jld["mlddepthtemp"]
    mlddepthsalinity = jld["mlddepthsalinity"]
    lon_rs = jld["lon_rs"]
    lat_rs = jld["lat_rs"]
    z = jld["z"]
    close(jld)
    jld = jldopen(datadir("en4_recs.jld2")); recs = jld["recs"]; recstime = jld["recstime"];close(jld)
    return mldtemp, mldsalinity, lon_rs, lat_rs, recs, recstime
end





end

