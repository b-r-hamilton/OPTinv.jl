import Pkg
Pkg.activate("../")
using PyPlot, OPTinv, UnitfulLinearAlgebra, Unitful, DimensionalData, LinearAlgebra, DrWatson

res = 10yr
corr = [-0.5, 0.25, 0.99]#[0.2K/permil,5K/permil, 15K/permil]
σθ = 11K * 2
σδ = 0.8permil / 1 *2   #std of global surface d18O in WOCE
num_iter = 1
Tᵤ = 1000yr:res:2000yr
modes = (:Mode1, :Mode2, :Mode3)

utest = firstguess(Tᵤ, [m for m in modes], σθ, σδ,0.5)


figure(figsize = (15,4))
for (c, color) in zip(corr, ["red", "blue", "green"])

    Tᵤ = 1000yr:res:2000yr
    modes = (:Mode1, :Mode2, :Mode3)
    u₀ = firstguess(Tᵤ,[m for m in modes], σθ, σδ, c)

    U = cholesky(u₀.Cnn.mat).U
    #sol = Dict()
    N = first(size(U))
    for n in 1:num_iter
        x = UnitfulMatrix(reshape(randn(N), (1, N)), fill(Unitful.NoUnits, 1), unitrange(U))
        xU = x*U #row vector
        arr = reshape(vec(parent(xU)) .* unit.(vec(xU)), size(u₀.y))
        da = DimArray(arr, u₀.y.dims)
        for (i,s) in enumerate([:θ, :δ])
            for (j, m) in enumerate(modes)
                subplot(2, length(modes), (i-1)*length(modes)+j)
                if i == 1
                    title(m)
                end
                plot(da[:, At(m), At(s)], color = color)
                #legend()
            end
        end
    end
end
tight_layout()
savefig(plotsdir("ex3_ts.png"))

figure(figsize = (15,4))
for (c, color) in zip(corr, ["red", "blue", "green"])

    Tᵤ = 1000yr:res:2000yr
    modes = (:Mode1, :Mode2, :Mode3)
    u₀ = firstguess(Tᵤ,[m for m in modes], σθ, σδ, c)

    U = cholesky(u₀.Cnn.mat).U
    #sol = Dict()
    N = first(size(U))
    @show c 
    for n in 1:num_iter
        x = UnitfulMatrix(reshape(randn(N), (1, N)), fill(Unitful.NoUnits, 1), unitrange(U))
        xU = x*U #row vector
        arr = reshape(vec(parent(xU)) .* unit.(vec(xU)), size(u₀.y))
        da = DimArray(arr, u₀.y.dims)

        for (j, m) in enumerate(modes)
            @show m 
            subplot(1, length(modes), j)
            title(m)
            slope = linearleastsquares(ustrip.(vec(da[:, At(m), At(:δ)])), ustrip.(vec(da[:, At(m), At(:θ)])))[1]            
            scatter(vec(da[:, At(m), At(:δ)]), vec(da[:, At(m), At(:θ)]), color = color, label = "ρ = " * string(ustrip.(c)) * ", slope = " * string(round.(slope, sigdigits = 2)) * "K/‰")

            println(slope)
            legend(loc = "lower left")
        end
    end
    println("")
end
tight_layout()
savefig(plotsdir("ex3_scatter.png"))
