module OPTinv

using DrWatson, Unitful, Revise, Statistics, BLUEs, LinearAlgebra,
    SparseArrays, TMI, Interpolations, UnitfulLinearAlgebra, Measurements,
    CSV, DataFrames, JLD2, NCDatasets, NaNMath, Downloads,
    Dates, DimensionalData, TMItransient, NMF, PythonCall, Distances 

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

export estimate

export loadLMR, loadOcean2k, loadThornalley, loadEN4, loadCESMLME

export record, loadEN4MLD 

export orthographic_axes, datadir

export computearea

#import units 
yr = u"yr"
kyr = u"kyr"
permil = Unitful.FixedUnits(u"permille")
K = u"K"
using DimensionalData: @dim, YDim, XDim, TimeDim, ZDim

@dim Cores "Cores" 
@dim Modes "Modes"
@dim StateVar "state variable"
export ccrs, cm, cfeature, cu, mpath, patches, cmap, mal, mticker, inset_axes

#import Python packages
const ccrs = pyimport("cartopy.crs")
const cm = pyimport("cmocean.cm")
const cfeature = pyimport("cartopy.feature")
const cu = pyimport("cartopy.util")
const mpath = pyimport("matplotlib.path")
const patches = pyimport("matplotlib.patches")
const cmap = pyimport("matplotlib.cm")
const mal = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable
const mticker = pyimport("matplotlib.ticker")
const inset_axes = pyimport("mpl_toolkits.axes_grid1.inset_locator").inset_axes

function __init__()

    PythonCall.pycopy!(ccrs,pyimport("cartopy.crs"))
    PythonCall.pycopy!(cm,pyimport("cmocean.cm"))
    PythonCall.pycopy!(cfeature,pyimport("cartopy.feature"))
    PythonCall.pycopy!(cu,pyimport("cartopy.util"))
    PythonCall.pycopy!(mpath,pyimport("matplotlib.path"))
    PythonCall.pycopy!(patches,pyimport("matplotlib.patches"))
    PythonCall.pycopy!(cmap,pyimport("matplotlib.cm"))
    PythonCall.pycopy!(mal,pyimport("mpl_toolkits.axes_grid1").make_axes_locatable)
    PythonCall.pycopy!(mticker, pyimport("matplotlib.ticker"))
    PythonCall.pycopy!(inset_axes, pyimport("mpl_toolkits.axes_grid1.inset_locator").inset_axes)

    println("OPTinv.jl: Python libraries installed")
end

export yr, permil, K,
    Ti, Cores, Modes, StateVar,
    DiagRule, OffDiagRule

datadir(x) = DrWatson.datadir(x)

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
struct inversion

details for inversion 

# Arguments
    - `cores`: Vector of symbols of cores to include
    - `modetype`: mostly deprecated, "svd" or "nnmf" (non-negative matrix fact.)
    - `name`: name of inversion
    - `color`: color for plotting 

# Example
    ```
    oldcores = Symbol.("MC" .* string.([26, 22, 13]) .* "A")
    inversion(oldcores, "svd", "old", "blue")
    ```
"""
struct inversion
    cores::Vector{Symbol}
    modetype::String
    name::String
    color::String
end

"""
struct solution

very big structure holding all estimate parameters and reconstructed fields 

# Arguments
    - `y`: DimEstimate of original sediment core data
    - `ỹ`: DimEstimate of reconstructed sediment core data according to solution
    - `ỹ₀`: DimEstimate of reconstructed sediment core data according to prior 
    - `u₀`: DimEstimate of prior 
    - `ũ`: DimEstimate of estimated mode magnitudes
    - `θ`: Vector of vectorized surface temperature ± 1σ unc. for each timestep 
    - `δ`: Vector of vectorized surface δ¹⁸O ± 1σ unc. for each timestep
    - `γ`: TMI.Grid
    - `spatialmodes`: V matrix 
    - `name`: String name of inversion 
    - `color`: String color for plotting
    - `E`: UnitfulMatrix 
    - `predict`: Function version of E
"""
struct solution
    y::DimEstimate
    ỹ::DimEstimate
    ỹ₀::DimEstimate
    u₀::DimEstimate
    ũ::DimEstimate
    θ::Vector
    δ::Vector
    γ::TMI.Grid
    spatialmodes::Matrix
    name::String
    color::String
    E::UnitfulMatrix
    predict::Function
end

"""
function invert

    completes inversion, checks to see if M, spatialmodes, modemags, exist, calculates
    if they don't, and runs inversion 

# Arguments
    - `x`: `inversion` type
"""
function invert(x::inversion; res = 10yr)::solution
    corenums_full = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
    corenums_sorted = [i for i in corenums_full if i in x.cores]
    filename = x.modetype * ".jld2"
    ℳ, spatialmodes = loadM(filename, corenums_sorted, res = res)
    display(ℳ)
    y = loadcores(corenums_sorted,res)
    ρ = 0.99
    T = Array(y.y.dims[1])
    Tᵤ = 500.0yr:res:T[end]

    if !isfile(DrWatson.datadir("modemags.jld2"))
        include(DrWatson.scriptsdir("CEdatanalysis/sigma.jl"))
    end
    jld = jldopen(DrWatson.datadir("modemags.jld2"))
    mags = jld[x.modetype * "mags"]
    
    σθ = vec(std(mags, dims = 1)) * K
    println("multiplying σθ by 4")  #fitting at 4 σ level
    σθ .*= 4
    #widening this up to 28 overfits the data! 
    
    #σδ = σθ ./ 14.84 .* permil/K
    σδ = σθ ./ 14.93 .* permil/K
    u₀ = firstguess(Tᵤ, ℳ.dims[2][:], σθ, σδ, ρ) #u₀ with correct T_u
    E, predict = loadE(filename, ℳ, Tᵤ, T, σθ, σδ, ρ, corenums_sorted, y.Cnn.ax; res = res)
    
    yde, ỹde, u₀de, ũ = solvesystem(y, u₀, E, predict)
    θ, δ = reconstructsurface(spatialmodes, ũ)
    
    TMIversion = "modern_180x90x33_GH11_GH12"
    γ = Grid(download_ncfile(TMIversion))
    
    ỹ₀ = DimEstimate(E*u₀de.v, parent(E*u₀de.C*E')*permil^2, y.y.dims)

    return solution(yde, ỹde, ỹ₀, u₀de, ũ, θ, δ, γ,spatialmodes, x.name, x.color, E, predict) 
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

"""
struct means

    structure for holding age model uncertainty information
"""
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

"""
function firstguess(T, modes, σθ::DimArray, σδ::Vector, ρ; fill_val = [0K, 0permil])

version where σθ is a DimArray to allow for Ocean2k constraint 
"""
function firstguess(T, modes, σθ::DimArray, σδ::Vector, ρ; fill_val = [0K, 0permil])
    z = zeros(length(T), length(modes))
    dims = (Ti(T), modes, StateVar([:θ, :δ]))
    da = DimArray(cat(z * K .+ fill_val[1], z * permil .+ fill_val[2], dims = 3), dims)
    units = unit.(vec(da))

    δdiagrules = [DiagRule(σδ[i]^2, (:, At(i), At(:δ))) for i in modes]
    θdiagrules = Vector{DiagRule}(undef, length(modes) * length(T))
    offdiagrules = Vector{OffDiagRule}(undef, length(modes) * length(T))
    for (i, t) in enumerate(T)
        da[At(t), :, At(:θ)] .= value.(vec(σθ[At(t), :]))
        inds = (i-1)*11+1:i*11
        σθ_ = vec(Measurements.uncertainty.(σθ[At(t), :, At(:θ)]))
        θdiagrules[inds] = [DiagRule(σθ_[j]^2, (At(t), At(j), At(:θ))) for j in modes]
        offdiagrules[inds] = [OffDiagRule(ρ*σθ_[j]*σδ[j], (At(t), At(j), At(:θ)), (At(t), At(j), At(:δ))) for j in modes]       
    end
    
    Cuu = fillcovariance(vec(units), vcat(θdiagrules, δdiagrules, offdiagrules), dims)
    @show isposdef(Cuu.mat) 
    return Est(da, Cuu) 
end

"""
function solvesystem

    optimizes system, prints cost function, returns optimized values 

# Arguments
    - `e`: Est object of observations
    - `u₀`: Est object of prior
    - `E`: UnitfulMatrix that maps between prior and observations
    - `predict`: function version of E 
"""
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
function linearleastsquares(x::UnitfulMatrix,y::UnitfulMatrix;C=nothing)
    urange = unitdomain(inv(C))
    udomain = [unitrange(y)[1]/unitrange(x)[1], unitrange(y)[1]]    
    E = UnitfulMatrix(hcat(parent(x), ones(length(x))), urange, udomain)
    F = inv(transpose(E)*inv(C)*E)*transpose(E)*inv(C)
    lls = F*y
    if !isnothing(C)
        Cnew = F*C*F'
        return lls, Cnew
    else
        return lls
    end
end

"""
function linearleastsquares(x,y)

    use linear algebra representation to find the linear least squares solution of
    y = a*x+b

    returns lls: vector where a = lls[1], b = lls[2]
"""
function linearleastsquares(x::Vector{T},y::Vector{T};C=nothing) where T
    E = hcat(parent(x), ones(length(x)))
    C = isnothing(C) ? I(length(x)) : C 
    F = inv(transpose(E)*inv(C)*E)*transpose(E)*inv(C)
    lls = F*y
    if !isnothing(C)
        Cnew = F*C*F'
        return lls, Cnew
    else
        return lls
    end
end

"""
function formattransientM

    returns ℳ as a DimArray
# Arguments
    - `arr`: Array
    - `τ`: lags vector
    - `modes`: vector of mode names
    - `cores`: vector of core names 
"""
function formattransientM(arr::Array, τ, modes, cores)
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

"""
function core_locations

    returns a NamedTuple with core names and locations 
"""
function core_locations()
    locations_file = datadir("EN539_locations.csv")
    cores_df = CSV.read(locations_file,DataFrame)
    locs = [(cores_df.Lon[i], cores_df.Lat[i], cores_df.Depth[i]) for i in 1:length(cores_df.Core)]
    cores = replace.(cores_df.Core, "EN539-"=>"")
    return NamedTuple{Tuple(Symbol.(cores))}(locs)
    
end

"""
function γreshape

    takes a vectorized wet surface points and reshapes them into a lat x lon grid 
"""
function γreshape(v::Vector{T}, γ) where T
    template = Array{T}(undef, size.(γ.axes)[1][1], size.(γ.axes)[2][1])
    template .= NaN
    template[γ.wet[:, :, 1]] .= v
    return template
end

"""
function loadM

    if a file is available, open, format, and return V and ℳ files
    otherwise, compute SVD and then propagate through TMItransient
"""
function loadM(filename::String, core_list::Vector{Symbol}; res = 10yr) 
    filepath = joinpath("../data/M", filename)
    !isdir("../data/M") && mkdir("../data/M")
    corelocs = core_locations()
    #load in variables from ex3.svdmodes.jl
    if isfile(filepath)
        println("opening pre-computed Vt and ℳ file") 
        jld = jldopen(filepath)
        ℳ = jld["ℳ"][begin:end-1, :, :]
        τ = jld["τ"][begin:end-1]yr
        spatialmodes = jld["SVD"].Vt
        close(jld)
    else
        surforigin, SVD = generatemodes(corelocs)
        res = 1
        τ = 0:res:1000
        println("computing ℳ, will take a while!") 
        ℳ = transientM(corelocs, SVD.Vt, τ)
        jldsave(filepath; ℳ, τ, res, surforigin, SVD)
        spatialmodes = SVD.Vt
        res = 10yr
    end
    cores = keys(core_locations())    
    modes = 1:size(spatialmodes)[1]
    #format into DimArray 
    ℳ_ = formattransientM(ℳ, τ, [m for m in modes], [c for c in cores])
    
    #M has 1 yr resolution
    ℳ_, τ = subsampletransientM(ℳ_, res)
    return ℳ_[:, :, At(core_list)], spatialmodes

end

"""
function loadcores

    generate Est objects for a list of cores 
"""
function loadcores(core_list::Vector{Symbol}, res = 10yr; dir = nothing, rules = nothing)
    dir = isnothing(dir) ? DrWatson.datadir() : dir 
    files = readdir(dir)
    ae_files = [f for f in files if occursin("ae", f)]
    d18O_files = [f for f in files if occursin("d18O", f)]
    if !isnothing(rules)
        rule_files = union([[f for f in files if occursin(r, f)] for r in rules]...)
        ae_files = intersect(ae_files, rule_files)
        d18O_files = intersect(d18O_files, rule_files)
    end
    
    directory_cores = Tuple([Symbol(split(f, ".")[1]) for f in ae_files])

    #re-order directory cores so they line up with core_list
    sort_indices = []
    for c in core_list
        ind = findall(x->x == c, directory_cores)
        if ! isempty(ind) 
            push!(sort_indices, ind[1])
        end
    end
    
    println("Pulling in core data, calculating covariance matrices") 
    @time e = formatbacon(joinpath.(dir, ae_files)[sort_indices], joinpath.(dir, d18O_files)[sort_indices], core_list, res = res)
    #add in measurement uncertainty 
    e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))
    @show isposdef(e.Cnn.mat)
    return e 
end

"""
function loadE(filename::String, ℳ::DimArray, Tᵤtarget, Ttarget,
               σθ, σδ::Vector, ρ,
               core_list::Vector{Symbol}, yax::Vector{Tuple})

    Calculates E matrix from ℳ array using an impulse response technique
    E is computed over an excessive time period (-500:2020) and then saved
    Each time we load in E, we will subset it to an appropriate time period 

# Arguments
    - `filename`
    - `ℳ`
    - `Tᵤtarget`
    - `Ttarget`
    - `σθ` : DimArray (changes over time) or vector (constant for each mode) 
    - `σδ` : vector (right now it can't change throughout time)
    - `ρ`: correlation coefficient
    - `core_list`: vector of core names as symbols
    - `yax`: Vector{Tuple}: covariance dims of Cnn 

# Output
    - `E`: UnitfulMatrix that maps vec(y) = E vec(x)
    - `predict`: functional form of `E` 
"""
function loadE(filename::String, ℳ::DimArray, Tᵤtarget, Ttarget,
               σθ, σδ::Vector, ρ,
               core_list::Vector{Symbol}, yax::Vector{Tuple}; res = 10yr)
    sv = (:θ, :δ)
    coeffs = NamedTuple{sv}([-0.224permil/K, 1])
    τ = [t for t in ℳ.dims[1]]

    #maximum time period we're interested in for this problem 
    T = collect(0yr:res:2020yr)
    Tᵤ = -500yr:res:2020yr
    u₀ = firstguess(Tᵤ, ℳ.dims[2][:], σθ, σδ, ρ)  #just for computing here 
    
    #this is the same as BLUEs.convolve, I just did it explicitly here 
    predict(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ[t .- τ .> minimum(Tᵤ)]]...) for s in sv]...) for t in T]..., dims = 2)', (Ti(T), Cores(core_list)))

    if res == 10yr
        filename = split(filename, ".")[1] .* "_E.jld2"
    else
        filename = split(filename, ".")[1] .* "_E" * string(res) * ".jld2"
    end
    
    filepath = joinpath("../data/M", filename)

    if !isfile(filepath)
        println("Computing E, will take a little while!") 
        @time E = impulseresponse(predict, u₀.y) #2 hrs
        jldsave(filepath; E)
    else
        println("Loading in pre-saved E file") 
        jld = jldopen(filepath)
        E = jld["E"]
        close(jld)
    end
    
    #this E matrix will only work for old cores 
    allcores = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9, 13, 14]) .* "A")
    Erange = vec(covariancedims((Ti(T), Cores(allcores))))
    modes = Array(ℳ.dims[2])
    Edomain = vec(covariancedims((Ti(Tᵤ), Modes(modes), StateVar([s for s in sv]))))
    xindices = [x[1] ∈ Ttarget && x[2]∈ core_list for x in Erange]
    yindices = [y[1] ∈ Tᵤtarget && y[2] ∈ modes for y in Edomain]
    predict_subset(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ[t .- τ .> minimum(Tᵤtarget)]]...) for s in sv]...) for t in Ttarget]..., dims = 2)', (Ti(Ttarget), Cores(core_list)))
    return UnitfulMatrix(parent(E[xindices, yindices]), fill(permil, sum(xindices)), unitdomain(E)[yindices]), predict_subset
end

"""
function ptobserve(pt::Tuple, γ, spatialmodes::Matrix, ũ)

    For a solution ũ, compute the value of surface variables at some location

# Arguments
    - `pt`: (lon, lat)
    - `γ`: associated TMI Grid
    - `spatialmodes`: Matrix
    - `ũ`: DimEstimate
# Output 
    - `de`: DimEstimate with dimensions (Time × Sv) 
"""
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

"""
function boxmean(fld::Array{T, N}, lat, lon, latsouth, latnorth, lonwest, loneast) where {T, N}

    Computes a spatial mean for a 3D array fld

# Arguments
    - `fld`: 3D Array, must be (lon x lat x time)
    - `lat`: vector of latitude values
    - `lon`: vector of longitude values
    - `latsouth`: southern most latitude value for chosen bb
    - `latnorth`: ^
    - `lonwest`: ^
    - `loneast`: ^
"""
function boxmean(fld::Array{T, N}, lat, lon, latsouth, latnorth, lonwest, loneast) where {T, N}
    #println("Using boxmean, but this should be replaced by `estimate`") 
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

    # Arguments
    - `γ`
    - `latsouth`: southern most latitude value for chosen bb
    - `latnorth`: ^
    - `lonwest`: ^
    - `loneast`: ^
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

"""
function orthographic_axes(latsouth,latnorth, lonwest, loneast, inc)

    Helper function for creating orthographic N. Atl. plots

    # Arguments
    - `latsouth`: southern most latitude value for chosen bb
    - `latnorth`: ^
    - `lonwest`: ^
    - `loneast`: ^

    # Output 
"""
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

"""
function geospatialsubset(mat, lat, lon, lat_target, lon_target)

    Subset a 3D matrix, accounts for roll-around 
"""
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
#=
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
    jld = jldopen(DrWatson.datadir("EN4_march_mld.jld2"))
    mldtemp = jld["mldtemp"]
    mldsalinity = jld["mldsalinity"]
    mlddepthtemp = jld["mlddepthtemp"]
    mlddepthsalinity = jld["mlddepthsalinity"]
    lon_rs = jld["lon_rs"]
    lat_rs = jld["lat_rs"]
    z = jld["z"]
    close(jld)
    jld = jldopen(DrWatson.datadir("en4_recs.jld2")); recs = jld["recs"]; recstime = jld["recstime"];close(jld)
    return mldtemp, mldsalinity, lon_rs, lat_rs, recs, recstime
end
=#
"""
function estimate(ũ::DimEstimate, spatialmodes::Matrix, variable::Symbol; spatialinds = :)

    Generate a mean surface estimate within some region
    Accounts for full covariance matrix

    # Arguments
    - `ũ`
    - `spatialmodes`
    - `variable`: Symbol
    - `spatialinds`: Vector of chosen indices within spatialmodes spatial indexing

    # Output
    - `da`: DimArray of estimated quantity 
    
"""
function estimate(ũ::DimEstimate, spatialmodes::Matrix, variable::Symbol; spatialinds = :, rolling = nothing, weights = nothing )    
    cdims = vec(covariancedims(ũ.dims))
    Tᵤ = Array(ũ.dims[1])
    if isnothing(weights) 
        Edag = spatialmodes[:, spatialinds] * fill(1/length(spatialinds), length(spatialinds))
    else
        Edag = spatialmodes[:, spatialinds] * weights
    end
    
    mat = fill(0.0, length(Tᵤ), length(cdims))
    for (i, t) in enumerate(Tᵤ)
        xinds = findall(x->x[1] == t && x[3] == variable, cdims)
        if !isnothing(rolling)
            tinds = i-rolling:1:i+rolling
            tinds = tinds[0 .< tinds .<= length(Tᵤ)]
            [mat[it, xinds] .= Edag ./ length(tinds) for it in tinds]
        else
            
            mat[i, xinds] = Edag
        end
        
    end
    un = variable == :θ ? K : permil
    
    mat = UnitfulMatrix(mat, fill(un, length(Tᵤ)), unitrange(ũ.v))
    y = vec(mat * ũ.v)
    C = mat * ũ.C * mat'
    #unit_unc = sqrt.(unitrange(C) ./ unitdomain(C))
    #unc = sqrt.(diag(parent(C)))
    return DimEstimate(y, C, (Ti(Tᵤ),))
end

"""
function estimate(ũ::DimEstimate, spatialmodes::Matrix, variable::Symbol, t1::Quantity, t2::Quantity; spatialinds = :)    

    Estimate the difference between a spatially-averaged quantity at 2 timesteps

    # Arguments
    - `ũ`
    - `spatialmodes`
    - `variable`: Symbol
    - `t1`: timestep 1
    - `t2`: timestep 2 
    - `spatialinds`: Vector of chosen indices within spatialmodes spatial indexing

    # Output
    - `da`: DimArray of estimated quantity 
    
"""
function estimate(ũ::DimEstimate, spatialmodes::Matrix, variable::Symbol, t1::Quantity, t2::Quantity, Tm1::Number, Tm2::Number; spatialinds = :)
    
    #helper 
    function convert_timestep(t)
        t = t == Tm1*yr ? t + 20yr : t
        t = t == Tm1*yr+10yr ? t + 10yr : t
        t = t == Tm2*yr ? t - 20yr : t
        t = t == Tm2*yr-10yr ? t - 10yr : t
        return t
    end
    t1 = convert_timestep(t1)
    t2 = convert_timestep(t2)

    cdims = vec(covariancedims(ũ.dims))
    Tᵤ = Array(ũ.dims[1])
    mat = fill(0.0, (1,length(cdims)))
    L = length(spatialinds)
    offsets = collect(-20yr:10yr:20yr)
    Edag = spatialmodes[:, spatialinds]*fill(1/L, L) ./ length(offsets)

    for offset in offsets
        mat[findall(x->x[1] == t1+offset && x[3] == variable, cdims)] .= -Edag
        mat[findall(x->x[1] == t2+offset && x[3] == variable, cdims)] .= Edag
    end
    
    un = variable == :θ ? K : permil
    mat = UnitfulMatrix(mat,fill(un, size(mat)[1]), unitrange(ũ.v))
    y = mat * ũ.v
    Css = mat * ũ.C * transpose(mat)
    return getindexqty(y,1) ± sqrt(getindexqty(Css,1,1))
end

"""
function transientM

    propagates some set of surface boundary conditions through the transient TMI matrix
    and observes at an interior location

    # Arguments
    - `corelocs`: NamedTuple of core locations
    - `patches`: NamedTuple of boundary conditions (should be TMI.BoundaryCondition)
    - `τ`: vector of lagged times to compute at 
"""
function transientM(corelocs, patches::NamedTuple, τ; TMIversion = "modern_180x90x33_GH11_GH12")
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    patches = NamedTuple{keys(patches)}([surfacepatch(patches[k]..., γ) for k in keys(patches)])
    N = length(corelocs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(corelocs[i],γ) for i in 1:N]
    modes = [k for k in keys(patches)]
    cores = [c for c in keys(corelocs)]
    arr = Array{Float64}(undef, length(τ), length(modes), length(cores)) # dimensions = Ti x cores x modes
    for (i, key) in enumerate(keys(patches))
        bc = patches[key]
        @time D̄ = stepresponse(TMIversion, bc, γ, L, B, τ, eval_func = observe, args = (wis, γ))
        Ḡ, _ = globalmean_impulseresponse(D̄, τ)
        Ḡ = hcat(Ḡ...)'    
        arr[begin:size(Ḡ)[1], i, :] = Ḡ
    end 
    return arr 
end

"""
function generatemodes

    based on core location water source maps, generate generate SVD or NNMF
# Arguments
    - `corelocs`: NamedTuple of core locations
# Optional Arguments
    - `TMIversion`: String of TMIversion to use, defaults to 2° modern
    - `func`: decomposition function to use, defaults to svd 
"""
function generatemodes(corelocs; TMIversion = "modern_180x90x33_GH11_GH12", func = svd)
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    tracers = hcat([surfaceorigin(corelocs[c], Alu, γ).tracer[γ.wet[:, :, 1]] for c in keys(corelocs)]...)
    tracers = 10 .^ tracers
    𝒩 = gaussiandistancematrix(γ, 1, 1000.0)
    mat = tracers' * 𝒩
    if func == nnmf 
        return mat, nnmf(mat, length(corelocs)).H
    elseif func == svd
        return mat,svd(mat)
    else
        return mat,func(mat)
    end    
end

"""
function transientM

    propagate matrix of boundary conditions through transient TMI
    and observe at core sites

# Arguments
    - `corelocs`: NamedTuple of core locations
    - `mat`: V matrix
    - `τ`: lags to compute over
# Optional Arguments
    - `TMIversion`
"""
function transientM(corelocs, mat::Matrix, τ; TMIversion = "modern_180x90x33_GH11_GH12")
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    
    bcs = [BoundaryCondition(γreshape(mat[i, :],γ), γ.axes[1:2], 0.0, 3, 1, γ.wet[:, :, 1], :i, string(i), string(i)) for i in 1:size(mat)[1]] 
    N = length(corelocs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(corelocs[i],γ) for i in 1:N]
    arr = Array{Float64}(undef, length(τ), length(bcs), length(corelocs)) # dimensions = Ti x cores x modes
    for (i, bc) in enumerate(bcs)
        @time D̄ = stepresponse(TMIversion, bc, γ, L, B, τ, eval_func = observe, args = (wis, γ))
        Ḡ, _ = globalmean_impulseresponse(D̄, τ)
        Ḡ = hcat(Ḡ...)'    
        arr[begin:size(Ḡ)[1], i, :] = Ḡ
    end
    return arr 
end

"""
function ginterp

    interpolate z(x,z) to (x2,y2) grid
"""
function ginterp(x,y,z::Array{Float32, 2}, x2, y2)
    x2 = matchvec(x2, x)
    y2 = matchvec(y2, y) 
    itp = LinearInterpolation((x, y), z)
    z2 = [itp(x,y) for x in x2, y in y2]
    return z2, x2, y2
end
"""
matchvec: trim a vector x1 so that it does not exceed extrema of x2
"""
matchvec(x1, x2) =  x1[x1 .< maximum(x2) .&& x1 .> minimum(x2)]



# Code to compute cell area, borrowed from TMI.cellarea, rephrased here
# to not need a Grid 
function zonalgriddist(lat, lon)
    dx = similar(lat)
    for j in eachindex(lat)
        dx[j] = haversine((lon[1],lat[j])
                         ,(lon[2],lat[j]))
    end
    return dx
end

function computearea(lon, lat, I)
    dx = zonalgriddist(lat,lon)
    dy = haversine((lon[1],lat[1]),(lon[1],lat[2]))
    area = Matrix{Float64}(undef,length(lon),length(lat))
    fill!(area,0.0)
    [area[I[ii][1],I[ii][2]] = dx[I[ii][2]] * dy for ii ∈ eachindex(I)]
    return area
end
end

