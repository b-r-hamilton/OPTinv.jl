#=
Are the sediment cores interpretable in the context of modern-day observations?
E.g., can we quantify the depth-structure of d18Oc?

This creates Figure 1, which shows
- the range of CE variability of each sediment core
- d18Oc computed from the steady-state TMI solution for T and d18Ow, using the Marchitto 2014 eqtns 
=#
using Pkg
Pkg.activate("../")
using OPTinv, Unitful, DimensionalData, DrWatson, LinearAlgebra, PythonPlot, TMI, DataFrames, GibbsSeaWater, Measurements, UnitfulLinearAlgebra

import OPTinv.yr, OPTinv.permil

#import sediment core data, interpolated to 10 yr resoltuion 
corenums_full = [28, 26, 25, 22, 21, 20, 19,10, 9,13,14]
core_list_full = Symbol.("MC".* string.(corenums_full) .* "A")

y = loadcores(core_list_full)
ycm = y.Cnn.ax
yuncd = Dict()
for c in core_list_full
    Cnnind = findall(x->x[1] == 1980yr && x[2] == c, ycm)[1]
    yuncd[c] = sqrt(getindexqty(y.Cnn.mat, Cnnind, Cnnind))
end

locs = core_locations()
depths = [locs[c][3] for c in keys(locs)]
lats = [locs[c][2] for c in keys(locs)]
lons = [locs[c][1] for c in keys(locs)]

#"observe" temperature and d18O at the sediment core sites in TMIss 
TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
θ = readfield(TMIfile, "θ", γ) #potential temp, degC
σθ = readfield(TMIfile, "σθ", γ)
d18O = readfield(TMIfile, "δ¹⁸Ow", γ) #permil VSMOW
σd18O = readfield(TMIfile, "σδ¹⁸Ow", γ) #permil VSMOW
Sp = readfield(TMIfile, "Sp", γ)

loctuple = [convert(Tuple{Float64, Float64, Float64}, locs[c]) for c in keys(locs)]
θobs = observe(θ, loctuple, γ)
σθobs = observe(σθ, loctuple, γ) 
d18Oobs = observe(d18O, loctuple, γ) #.- 0.27
σd18Oobs = observe(σd18O, loctuple, γ)

#=
Spobs = observe(Sp, loctuple, γ)
pobs = gsw_p_from_z.(-1 .* depths, lats)
SAobs = gsw_sa_from_sp.(Spobs, pobs, lons, lats)
CTobs = gsw_ct_from_pt.(SAobs, θobs)
T = gsw_t_from_ct.(SAobs, CTobs,pobs)  

effd18Oc_cibs = @. -0.224 * T + 3.53 + d18Oobs - 0.27
effd18Oc_uvi = @. -0.231 * T + 4.03 + d18Oobs - 0.27 - 0.47
effd18Oc_cibs_fs = @. -0.225 * T + 3.50 + d18Oobs - 0.27
effd18Oc_uvi_fs = @. -0.207 * T + 3.75 + d18Oobs - 0.27 - 0.47
=#

#equation to follow the Marchitto 2014 recommendation of applying the Cibs eqtns to both Cibs and Uvi records 
effd18Oc_cibs = (-0.224 ± 0.002) * (θobs .± σθobs) .+ (3.53 ± 0.02) .+ (d18Oobs .± σd18Oobs) .- 0.27
#effd18Oc_uvi = (-0.231 ± 0.004) * (θobs .± σθobs) .+ (4.03 ± 0.03) .+ (d18Oobs .± σd18Oobs) .- 0.27 .- (0.47 ± 0.04)

#effd18Oc_cibs_fs =  (-0.225 ±0.006) * (θobs .± σθobs) .+ (3.50 ± 0.07) .+ (d18Oobs .± σd18Oobs) .- 0.27
#effd18Oc_uvi_fs = (-0.207 ± 0.007) * (θobs .± σθobs) .+ (3.75 ± 0.08) .+ (d18Oobs .± σd18Oobs) .- 0.27 .- (0.47 ± 0.04)

#plot d18Oc v. depth 
fig = figure(figsize = (5,7)) 
ax = gca()
T = ustrip.(Array(dims(y.y)[1]))
for c in keys(locs)
    ytmp = ustrip(y.y[:, At(c)])
    ytmp = c ∈ [:MC28A, :MC26A, :MC25A] ? ytmp .- 0.47 : ytmp
    d = fill(locs[c][3], length(ytmp)) 
    s = scatter(ytmp, d, c = T, cmap = "rainbow", s = (2050 .- T) ./ 5, zorder = 10)
    #s = scatter([ytmp], [d], c = "red", s = 50)
    if c == keys(locs)[end]
        cb = colorbar(s, location = "top")
        cb.set_label("Time [years CE]", fontsize = 12)
        
    end
end
T
ms = 15
scatter(effd18Oc_cibs, depths, color = "black", label = "Effective " * L"\mathrm{\delta}^{18}\mathrm{O}_\mathrm{calcite}", markersize = ms, capsize = 3, fmt = "X", zorder = 0, linewidth = 3)
#scatter(effd18Oc_uvi_fs[begin:3], depths[begin:3], color = "lightgray", label = "Effective " * L"\mathrm{\delta}^{18}\mathrm{O}_\mathrm{calcite}", markersize = ms, marker = "X", markeredgecolor = "black", capsize = 5)
#scatter(effd18Oc_cibs_fs[4:end], depths[4:end], color = "lightgray", label = "Effective " * L"\mathrm{\delta}^{18}\mathrm{O}_\mathrm{calcite}", markersize = ms, marker = "X", markeredgecolor = "black", capsize = 5)
#plot(effd18Oc_uvi, depths)
ax.invert_yaxis()
depths
xl = ax.get_xlim()
yt = ax.get_yticks()
yl = ax.get_ylim()
#hlines(y = depths, xmin = xl[1], xmax = xl[2], color = "gray", zorder = 0)
tx = ax.twinx()
tx.set_xlim(xl)
tx.set_yticks(depths, string.(keys(locs)), fontsize = 12)
tx.set_ylim(yl)

#=

ctd_dir = "/home/brynn/Documents/JP/rdata/EN539CTD/"
d18Ofp = joinpath(ctd_dir, "EN539_O18water.xlsx")
function get_d18O(filepath, index, ctd_col::Char, depth_col::Char, data_col, names) 
    xf = XLSX.readxlsx(filepath)
    sn = XLSX.sheetnames(xf)
    @show ctd_num = xf[sn[1] * "!" * ctd_col*':'*ctd_col][:, 1]
    depth = xf[sn[1] * "!" * depth_col * ':' * depth_col][:, 1]
    
    recs = [xf[sn[1] * "!" * dc*':' *dc][:, 1] for dc in data_col]
    nm = intersect(findall(x->(!).(ismissing(x)), ctd_num), findall(x->x isa Number, ctd_num))
    indices = [findall(x->x == i, ctd_num[nm])[1] for i in index if i ∈ ctd_num[nm]]

    return DataFrame(hcat(ctd_num[nm][indices], depth[nm][indices],[r[nm][indices] for r in recs]...), ["CTD","depth", names...])
end

function get_ctd(filepath) 
    xf = XLSX.readxlsx(filepath)
    sn = XLSX.sheetnames(xf)
    depth = xf[sn[1] * "!" * "U:U"][2:end, 1]
    temp = xf[sn[1] * "!" * "Z:Z"][2:end, 1]
    sal = xf[sn[1] * "!" * "V:V"][2:end, 1]
    return DataFrame(hcat(depth, temp, sal), ["depth", "temp", "sal"])
end

ctd12 = get_ctd(joinpath(ctd_dir, "EN539_ctd12.xlsx"))
ctd17 = get_ctd(joinpath(ctd_dir, "EN539ctd17.xlsx"))
ctd27 = get_ctd(joinpath(ctd_dir, "EN539ctd27.xlsx"))


core_nums = [28, 26, 25, 22, 21, 20, 19, 10, 9, 13, 14] 
d18O = get_d18O(d18Ofp, core_nums, 'H', 'I', ['O'], ["d18O"])
θmulticore = get_d18O(joinpath(ctd_dir, "en539_salinity2.xlsx"), core_nums, 'B', 'C', ['H', 'I'], ["θ", "S"])
display(d18O)
display(θmulticore)


inds = vcat([findall(x->x==c, d18O[!, "CTD"]) for c in θmulticore[!, "CTD"]]...)

multicored18Oc = -0.224 * θmulticore[inds,"θ"] .+ d18O[inds, "d18O"] .+ 3.53 .- 0.27
#scatter(multicored18Oc, θmulticore[inds, "depth"], color = "blue", s = 25, marker = "x")
=#
ax.set_xticks(2.2:0.2:2.9, 2.2:0.2:2.9, fontsize = 12)
ax.set_xlim(xl)
ax.set_xlabel(L"\delta^{18}\mathrm{O}_\mathrm{calcite}" *" [‰]", fontsize = 15)
ax.set_yticks(2200:-200:1000, 2200:-200:1000, fontsize = 12)
ax.set_ylabel("Depth [m]", fontsize = 15)
tight_layout()
savefig(plotsdir("offsets.png"), dpi = 600)

using Statistics
@show (mean(y.y, dims = Ti))[1, :] .- effd18Oc_cibs .* permil
@show (maximum(y.y, dims = Ti) .- minimum(y.y, dims = Ti))[1,:]

#=
figure()
subplot(1,2,1)
using GibbsSeaWater
stringcores = "MC" .* string.(θmulticore[!, "CTD"]) .* "A"
lats = [locs[Symbol(sc)][2] for sc in stringcores]
lons = [locs[Symbol(sc)][1] for sc in stringcores]
p = gsw_p_from_z.(-1 .* θmulticore[!, "depth"], lats)
sa = gsw_sa_from_sp.(θmulticore[!, "S"], p, lons, lats)
pt = gsw_pt_from_t.(sa, θmulticore[!, "θ"], p, 0)
scatter(pt, θmulticore[!, "depth"])
[plot(c[!, "temp"], c[!, "depth"], label = lab ) for (c, lab) in zip([ctd12, ctd17, ctd27], "CTD" .* string.([12,17,27]))]
scatter(θmulticore[!, "θ"], θmulticore[!, "depth"], label = "Multicore", color = "black")
gca().invert_yaxis()
ylabel("Depth [m]")
xlabel("Potential Temperature [degC]")
xl = gca().get_xlim()

legend()
subplot(1,2,2)
[plot(c[!, "sal"], c[!, "depth"], label = lab ) for (c, lab) in zip([ctd12, ctd17, ctd27], "CTD" .* string.([12,17,27]))]
scatter(θmulticore[!, "S"], θmulticore[!, "depth"], label = "Multicore", color = "black")
gca().invert_yaxis()
ylabel("Depth [m]")
xlabel("Salinity [psu]")
legend()
tight_layout()


=#
