module OPTinv

using DrWatson, Unitful, Revise, Statistics, BLUEs, LinearAlgebra,
    SparseArrays
import Interpolations.deduplicate_knots! as deduplicate! 
import Interpolations.LinearInterpolation as LI
import DelimitedFiles.readdlm 
export readbaconfile, interpolatebacon, meanbacon, covariancebacon,
    formatbacon, covariancedims

yr = u"yr"
permil = u"permille"


using DimensionalData
using DimensionalData: @dim, YDim, XDim, TimeDim, ZDim

@dim Cores XDim "Cores" #can't use "Core" because that's the import for base Julia

struct CovMat
    ax::Vector{Tuple} #maps on to axes 
    mat::Symmetric #to be efficient 
end

struct Est
    y::DimArray 
    Cnn::CovMat
end


function readbaconfile(age_file::String, d18O_file::String; matrix=false, res = 5yr, core = :replaceme)
    ages = readdlm(age_file, ',', Float64, header = true, skipblanks = true)[1][:, 2:end]yr
    d18O = readdlm(d18O_file, ',', Float64, header = true)[1][:, 2]permil

    if matrix 
        return ages, d18O
    else
        T, int = interpolatebacon(ages, d18O, res);
        μ, T = meanbacon(int, T, ages)
        Cnn = covariancebacon(int, T, μ, ages, core)
        de = formatbacon(μ, T, Cnn, core)
        return de 
    end
end

function readbaconfile(age_file::Vector{String}, d18O_file::Vector{String},
                       cores::Vector{Symbol}; res = 5yr)
    #compute the means first
    #because we want to know exactly what indices to compute Cnn for 
    function mean_shortcut(af, df)
        a, d = readbaconfile(af, df, matrix = true)
        T, int = interpolatebacon(a, d, res)
        μ, T = meanbacon(int, T, a)
        return μ, T, int, a
    end
    
    output = [mean_shortcut(f...) for f in zip(age_file, d18O_file)]
    T = intersect([o[2] for o in output]...)
    mus = [mu[findall(x->x∈T, Ts_)] for (mu, Ts_) in output]

    dim = (Ti(T), Cores(cores))
    Cnn_vector = [covariancebacon(output[i][3], T, output[i][1], output[i][4], cores[i]) for i in 1:length(cores)]

    dims = (Ti(T), Cores(cores))
    
    data = DimArray(hcat(mus...), dims)
    cdims = covariancedims(dims)

    Cnn = combine(Cnn_vector, vec(cdims))
    
    return Est(data, Cnn) 
    
end

function combine(Cnn_vector::Vector{CovMat}, target_dims::Vector{Tuple})
    N = length(target_dims)

    #this potentially could be faster by using blockdiag and sparse
    
    mat = zeros(N, N)
    helper_da = DimArray(1:N, X(target_dims))
    for Cnn in Cnn_vector
        ind = helper_da[At(Cnn.ax)]
        mat[ind, ind] .= Cnn.mat
    end
    return CovMat(target_dims, Symmetric(sparse(mat)))
end

function covariancedims(d::Tuple)
    arr = Array{Tuple}(undef, length.(d)...)
    for ind in CartesianIndices(arr)
        arr[ind] = Tuple([d[i][ind[i]] for i in 1:length(ind)])
    end
    return arr
end


function interpolatebacon(ages::Matrix, d18O::Vector,
                          res::Quantity)
    utype = typeof(first(ages))
    numr = size(ages)[2]
    Tmin = minimum(ages)
    Tremain = Tmin % res
    Tfloor = Tremain ≥ 0.5res ? Tmin + (res - Tremain) : Tmin - Tremain
    T = Tfloor:res:floor(utype, maximum(ages)) 
    d = [LI(deduplicate!(reverse(ages[:, j])), reverse(d18O)) for j in 1:numr]

    return T, d 
end

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

function covariancebacon(d, T, μ, ages, c)
    Cnn = fill(NaN * permil, length(T), length(T))
    for (μ1, (i, t1)) in zip(μ, enumerate(T)) 
        Cnn[i, i:end] .= [covary(d, ages, t1, t2, μ1, μ2) for (μ2, t2) in zip(μ[i:end], T[i:end])]
    end
    axis = [(t,c) for t in T] 
    return CovMat(axis, Symmetric(Cnn))
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
    
    #vals = [(d[k](t₁) - μ₁) * (d[k](t₂) - μ₂) for k in collect(1:length(d))[yes]]
    if !isempty(vals) return sum(vals)/length(vals) else return 0permil^2  end 
end

function formatbacon(μ, T, Cnn, c)
    if length(size(μ)) !== 2
        μ = reshape(μ, (length(μ), 1))
    end
    if c isa Symbol
        c = [c]
    end
    
    da = DimArray(μ, (Ti(T), Cores(c)))
    return Est(da, Cnn) 
end

end
