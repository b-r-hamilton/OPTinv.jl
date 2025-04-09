#=
Script to generate OPT-3-GLOBAL and OPT-3-MODE2 estimates and custom plots 
=#
import Pkg;Pkg.activate("../")
using OPTinv, PythonPlot, JLD2, DrWatson, Statistics, Unitful, DimensionalData, PythonPlotExt, BLUEs, TMI


oldcores = Symbol.("MC" .* string.([26, 22, 13]) .* "A")
xall = inversion(oldcores, "global", "all", "#5D3A9B")
xold = inversion(oldcores, "mode2", "old", "#E66100")
corenums_full = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
solutions = Vector{solution}(undef, 2)
for (i,x) in enumerate([xall, xold])
    corenums_sorted = [i for i in corenums_full if i in x.cores]
    filename = x.modetype * ".jld2"
    #template for loadM
    ℳtemp, spatialmodestemp = loadM("svd.jld2", corenums_sorted)
    N = 11113
    #global

    if x.modetype == "global"
    #sum(spatialmodestemp.^2, dims = 2) should be 1, dims should be 1 x 11113
        spatialmodes = reshape(fill(sqrt(1/N), N), (1, N))
    elseif x.modetype == "mode2"
    #just mode 2 solution 
        spatialmodes = reshape(spatialmodestemp[2, :], (1,N))
    end
    
    τ = 0:1:1000
    filepath = joinpath("../data/M", filename)
    if !isfile(filepath)
        ℳ = OPTinv.transientM(core_locations(), spatialmodes, τ)
        jldsave(filepath; ℳ, τ, spatialmodes)
    else
        jld = jldopen(filepath)
        ℳ = jld["ℳ"]
    end
    cores = keys(core_locations())    
    modes = 1:size(spatialmodes)[1]
    #format into DimArray 
    ℳ = formattransientM(ℳ, τ * yr, [m for m in modes], [c for c in cores])
    #M has 1 yr resolution
    ℳ, τ = subsampletransientM(ℳ, 10yr)
    ℳ = ℳ[:, :, At(x.cores)]
    y = loadcores(corenums_sorted)
    ρ = 0.99
    T = Array(y.y.dims[1])
    res = unique(diff(T))[1]
    Tᵤ = 500.0yr:res:T[end]
    
    if x.modetype == "global"
        jld = jldopen(DrWatson.datadir("modemags_global.jld2"))
        mags = jld["svdmags"]
        σθ = vec(std(mags, dims = 1)) * K
        σθ .*= 8
        σδ = σθ ./ 14.84 .* permil/K
        u₀ = firstguess(Tᵤ, ℳ.dims[2][:], σθ, σδ, ρ) #u₀ with correct T_u
    elseif x.modetype == "mode2"
        jld = jldopen(DrWatson.datadir("modemags.jld2"))
        mags = jld["svdmags"]

        σθ = vec(std(mags, dims = 1)) * K
        σθ .*= 8
        σδ = σθ ./ 14.84 .* permil/K
        u₀ = firstguess(Tᵤ, ℳ.dims[2][:], [σθ[2]], [σδ[2]], ρ) #u₀ with correct T_u
    end
    
    E, predict = loadE(filename, ℳ, Tᵤ, T, σθ, σδ, ρ, corenums_sorted, y.Cnn.ax)

    yde, ỹde, u₀de, ũ = solvesystem(y, u₀, E, predict)
    θ, δ = reconstructsurface(spatialmodes, ũ)

    TMIversion = "modern_180x90x33_GH11_GH12"
    γ = Grid(download_ncfile(TMIversion))

    #yes this is bad to force the permil, but I'll come back to it...
    ỹ₀ = DimEstimate(E*u₀de.v, parent(E*u₀de.C*E')*permil^2, y.y.dims)#predict(u₀de.x) 
    solutions[i] = solution(yde, ỹde, ỹ₀, u₀de, ũ, θ, δ, γ,spatialmodes, x.name, x.color, E, predict)
end
suffix = "globmode2"

min_age = [minimum(s.y.dims[1]) for s in solutions]
max_age = [maximum(s.y.dims[1]) for s in solutions]
xl = [minimum(min_age), maximum(max_age)]

locs = core_locations()
#recon @ cores
figure(figsize = (10,3))
for (i, c) in enumerate(oldcores)
    ax = subplot(1,3,i)
    plot(solutions[1].ỹ.x[:, At(c)], color = solutions[1].color)
    if c ∉ oldcores
        plot(solutions[1].y.x[:, At(c)], color = "black",lwcentral = 2, alpha = 0.2)
    else
        plot(solutions[2].ỹ.x[:, At(c)], color = solutions[2].color, lwcentral = 2)
        plot(solutions[2].y.x[:, At(c)], color = "black", lwcentral = 2)
    end
    
    xticks(1000:250:1750, fontsize = 12)
    xlabel("Time [years CE]", fontsize = 15)
    
    if i ∈ [1]
        yticks(-0.2:0.1:0.1, fontsize = 12) 
        ylabel(L"\mathrm{\delta}^{18}\mathrm{O}_\mathrm{calcite}" * " [‰]", fontsize = 15)
    else
        ylabel("")
        yticks(-0.2:0.1:0.1, ["", "", "", ""])
    end
        
    title(string(c) * ", " * string(locs[c][3]) * "m", fontsize = 15)
    ax.invert_yaxis()
end
tight_layout()
savefig(DrWatson.plotsdir("reconsol" * "global" * ".png"), dpi = 600)


figure();
#for sol in solutions
plot(estimate(solutions[1].ũ, solutions[1].spatialmodes, :θ, spatialinds = 1:1:11113, rolling = 2).x, color = solutions[1].color) #this line only works for `global` because `Mode2` doesn't have a unique T value at each timestep 
#end
xlim(800,1970)
grid()
ylabel("Global Surf. Temp. Anom. [K]", fontsize = 15)
xlabel("Time [years CE]", fontsize = 15) 


figure(figsize = (8,4))
for (i, file) in enumerate(["global.jld2", "mode2.jld2"])
    subplot(1,2,i) 
    filepath = joinpath("../data/M", file)
    jld = jldopen(filepath)
    ℳ = jld["ℳ"]
    cores = keys(core_locations())    
    modes = 1:1
    #format into DimArray 
    ℳ = formattransientM(ℳ, collect(1:1:1001) * yr, [m for m in modes], [c for c in cores])

    select_cores = [:MC26A, :MC22A, :MC13A]
    colors = ["orange", "blue", "purple"]
    for (j,(c, color)) in enumerate(zip(select_cores, colors))
        plot(ℳ[:, At(1), At(c)], color = color, label = c, zorder = 11-j, linewidth = (j+1)*2, alpha = 1)
    end
    xlim(0,200)
    ylim(-0.0004, 0.00025)
    xlabel("Lagged Time [years]",fontsize = 15)
    if i == 1 
        ylabel("Mode Mag. []", fontsize = 15)
    end
    text(x = 175, y =  -0.00037, s = ["A", "B"][i], fontsize = 30, fontweight = "bold")
end
tight_layout()
savefig(plotsdir("Mode2.png"))

