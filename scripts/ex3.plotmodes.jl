
import Pkg
Pkg.activate("../")


using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions
using PyCall
cmocean = pyimport("cmocean.cm")
mal = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable
close("all")
import Measurements.value as value
import OPTinv.Est


ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")
mticker = pyimport("matplotlib.ticker")

core_list = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13, 14]) .* "A")
#load in ℳ and surface spatial modes
#=
filename = "nmf.jld2"
ℳnmf, spatialmodesnmf = loadM(filename, core_list)
=#
filename = "svd.jld2"
ℳsvd, spatialmodessvd = loadM(filename, core_list)

filepath = joinpath("../data/M", filename)
#load in variables from ex3.svdmodes.jl 
jld = jldopen(filepath)
S = jld["SVD"].S
pervar = S.^2/sum(S.^2)


if ! @isdefined A 
    TMIversion = "modern_180x90x33_GH11_GH12"
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
end

lon = ustrip.(γ.lon); lat = ustrip.(γ.lat)
#lev = -0.05:0.01:0.05
#for (sm, t) in zip([spatialmodessvd, spatialmodesnmf], ["SVD", "NMF"])
for (sm, t) in zip([spatialmodessvd], ["SVD"])
    fig = figure(figsize = (9,12))
    for i in 1:11 
        ax = fig.add_subplot(3,4,i,projection = ccrs.PlateCarree())
        if t == "SVD"
            title(string(i) * ": " * string(round.(pervar .* 100, sigdigits = 3)[i]) * "%")
        else
            title(i)
            
        end
        
        plotme = γreshape(sm[i, :], γ)'
        plotme, lon = cu.add_cyclic_point(plotme, coord = ustrip.(γ.lon))
        pm = plotme[(!).(isnan.(plotme))]
        vm = maximum(abs.([minimum(pm), maximum(pm)]))

        cf = ax.contourf(lon, lat, plotme, cmap = cm.balance, vmin = -vm, vmax = vm, levels = 8)


        
        ax.set_extent([-80, 30, -90, 90])
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))


        gl = ax.gridlines(draw_labels = true)
        gl.top_labels = false
        gl.right_labels = false
        gl.bottom_labels = false
        gl.left_labels = false
        if i ∈ [1,5,9]
            gl.left_labels = true
        end
        if i ∈ [9, 10, 11, 8]
            gl.bottom_labels = true
        end
        #ticks = floor.(range(minimum(pm), maximum(pm), length = 4), sigdigits = 1)
        cb = colorbar(cf,ax = ax, location = "bottom", fraction = 0.046, pad = 0.1)
        cb.ax.tick_params(rotation = 45, labelsize = 10)
    end

    tight_layout()
    savefig(plotsdir(t * "_modes.png"))
end
tight_layout()

select_cores = [:MC26A, :MC22A, :MC13A]
colors = ["orange", "blue", "purple"]
cm = get_cmap("rainbow")
figure(figsize = (10,10))
for i in 1:11
    ax = subplot(4,3,i)
    title(string(i), fontsize = 15)
    #=
    for (j,c) in enumerate(ℳsvd.dims[3])
        plot(ℳsvd[:, At(i), At(c)], color = cm((-j+11)/12), label = c, zorder = 11-j)
    end
    =#
    for (j,(c, color)) in enumerate(zip(select_cores, colors))
        plot(ℳsvd[:, At(i), At(c)], color = color, label = c, zorder = 11-j, linewidth = (j+1)*2, alpha = 1)
    end

    #legend()
    xlim(0,200)
    #ylim(-0.002, 0.004)
    yt = ax.get_yticks()
    ax.set_yticks(yt, yt, fontsize = 12) 
    if i ∉ [10, 11, 9]
        xticks([])
        xlabel("")
    else
        xlabel("Lagged Time [years]", fontsize = 15)
        xt = ax.get_xticks()
        xt = convert(Vector{Int}, xt) 
        ax.set_xticks(xt, xt, fontsize = 12)
    end
    if i ∉ [1,4,7,10]
        #yticks([])
    else
        ylabel("Mode Mag. []", fontsize = 15)
    end
    #grid(axis = "x")
end
tight_layout() 
savefig(plotsdir("modemags.png"))


figure()
U = jld["U"]
Vt = jld["modes"]
s = U * diagm(S) * Vt
loading = transpose(U)*s

figure()
cm = get_cmap("rainbow") 
for i in 1:11
    plot(U[:, i], Array(ℳsvd.dims[3]), label = i, color = cm((11-i)/11))
    gca().invert_yaxis()
end
legend()


figure(figsize = (8,10))
for i in 1:11
    subplot(4,3,i,projection = ccrs.PlateCarree())
    plotme = γreshape(loading[i, :], γ)'
    pm = maximum(abs.([NaNMath.minimum(plotme),NaNMath.maximum(plotme)]))
    cf = contourf(γ.lon, γ.lat, plotme, vmin = -pm, vmax = pm, cmap = "coolwarm")
    colorbar(cf)
end

diagm(S)
mc28a = γreshape(vec(U[1, :]' * diagm(S)*Vt), γ)
figure()
contourf(γ.lon, γ.lat, mc28a')

vec(U[1, :]' * diagm(S)*Vt)
U*diagm(S)*Vt

mat = randn(11,11113)
u,s,v = svd(mat)
u*diagm(s)*v'
mat
