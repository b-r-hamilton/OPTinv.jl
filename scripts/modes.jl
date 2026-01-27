import Pkg
Pkg.activate("../")
using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PythonPlot, Revise, DrWatson, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions, PythonPlotExt, PythonCall
import Measurements.value as value
import OPTinv.Est

# ========================== SURFACE SPATIAL MODE MAPS ====================== #
core_list = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13, 14]) .* "A")
filename = "svd.jld2"
ℳsvd, spatialmodessvd = loadM(filename, core_list)

filepath = joinpath("../data/M", filename)
#load in variables from ex3.svdmodes.jl 
jld = jldopen(filepath)
S = jld["SVD"].S
pervar = S.^2/sum(S.^2)

pvtrans = jldopen(joinpath("../data/", "pcvartrans.jld2"))["pcvar"]

if ! @isdefined A 
    TMIversion = "modern_180x90x33_GH11_GH12"
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
end

lon = ustrip.(γ.lon); lat = ustrip.(γ.lat)
for (sm, t) in zip([spatialmodessvd], ["SVD"])
    fig = figure(figsize = (11,12))
    for i in 1:11 
        ax = fig.add_subplot(3,4,i,projection = ccrs.PlateCarree())
        if t == "SVD"
            title(string(i) * ": ss = " * string(round.(pervar .* 100, sigdigits = 2)[i]) * "%" *"; t = " * string(round.(pvtrans, sigdigits = 2)[i]) * "%", fontsize = 12)
        else
            title(i) 
        end
        plotme = γreshape(sm[i, :], γ)'
        plotme, lon = cu.add_cyclic_point(plotme, coord = ustrip.(γ.lon))
        pm = plotme
        vm = maximum(abs.([minimum(sm[i,:]), maximum(sm[i, :])]))
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
        cb = colorbar(cf,ax = ax, location = "bottom", fraction = 0.03, pad = 0.1)
        cb.ax.tick_params(rotation = 45, labelsize = 10)
    end
    tight_layout()
    savefig(plotsdir(t * "_modes.png"), dpi = 600)
end

# =========== PLOT THE IMPULSE RESPONSE OF EACH MODE AT THREE SELECT CORES ===== #
select_cores = [:MC26A, :MC22A, :MC13A]
colors = ["orange", "blue", "purple"]
figure(figsize = (10,10))
for i in 1:11
    ax = subplot(4,3,i)
    title(string(i), fontsize = 15)
    for (j,(c, color)) in enumerate(zip(select_cores, colors))
        plot(ℳsvd[:, At(i), At(c)], color = color, label = c, zorder = 11-j, linewidth = (j+1)*2, alpha = 1)
    end
    xlim(0,200)
    yt = ax.get_yticks()
    xt = ax.get_xticks()
    ax.set_yticks(yt, yt, fontsize = 12) 
    if i ∉ [10, 11, 9]
        ax.set_xticks(xt, fill("", length(xt)))
        xlabel("")
    else
        xlabel("Lagged Time [years]", fontsize = 15)
        xt = pyconvert(Vector{Int}, xt) 
        ax.set_xticks(xt, xt, fontsize = 12)
    end
    if i ∉ [1,4,7,10]
    else
        ylabel("Mode Mag. []", fontsize = 15)
    end
end
tight_layout() 
savefig(plotsdir("modemags.png"), dpi = 600)

# ================================ PLOT SINGULAR VALUES ================= #
s = jldopen("../data/M/svd.jld2")["SVD"]
U = s.U #cores × modes
locs = core_locations()
depths = [locs[c][3] for c in keys(locs)]
figure(figsize = (10, 4))
for i in 1:11
    ax = subplot(1, 11, i) 
    scatter(U[:, i], depths, c = U[:, i], cmap = cm.balance, vmin = -1, vmax = 1,edgecolors = "gray", 100)
    ax.invert_yaxis()
    xlim(-1, 1)
    xl = ax.get_xlim()
    yt = collect(1000:200:2200)
    hlines(y = depths, xmin = xl[0], xmax = xl[1], color = "gray", zorder = 0)
    tx = ax.twinx()
    tx.set_xlim(xl)
    
    if i != 1
        ax.set_yticks(yt, fill("", length(yt)))
        yl = ax.get_ylim()

    else
        ax.set_ylabel("Depth [m]", fontsize = 15)
        yt = convert(Vector{Int}, yt)
        ax.set_yticks(yt, yt, fontsize = 12)
        yl = ax.get_ylim()
    end
    if i == 11
        tx.set_yticks(depths, string.(keys(locs)), fontsize = 12)
        tx.set_ylabel("Sediment Core", fontsize = 15)
    else
        tx.set_yticks(depths, fill("", length(depths)))
end
    tx.set_ylim(yl)
    title(string(i), fontsize = 15)
    ax.set_xticks([-1, 0, 1], [-1, 0, 1], fontsize = 12)
end
tight_layout()

savefig(plotsdir("U.png"), dpi = 600)
