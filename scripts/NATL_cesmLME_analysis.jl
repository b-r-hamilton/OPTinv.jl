import Pkg;Pkg.activate("../")
using JLD2, PyPlot, Dates, Statistics, DrWatson, NCDatasets, OPTinv, RollingFunctions, DateFormats, Measurements, LinearAlgebra

jld = jldopen(datadir("NATL_cesmLME.jld2"))
mat = jld["mat"]
exp = jld["exp"]
inds = jld["inds"]
indnames = ["N. Atl.", "SPNA", "NS"]
ts = collect(DateTimeNoLeap(850, 2,1):Month(1):DateTimeNoLeap(2006, 1,1))

ts_10yr_left = collect(DateTimeNoLeap(850, 2,1):Year(10):DateTimeNoLeap(1996, 1,1))
ts_10yr_right = collect(DateTimeNoLeap(860, 2,1):Year(10):DateTimeNoLeap(2006, 1,1))
mat_binned = Array{Float64}(undef, size(mat)[1:2]..., length(ts_10yr_left))
for i in 1:length(ts_10yr_left) 
    inds = findall(x->ts_10yr_left[i] < x < ts_10yr_right[i], ts)
    summer = findall(x->x ∈ [7,8,9], month.(ts))
    inds = intersect(inds, summer) 
    mat_binned[:, :, i] .= mean(mat[:, :, inds], dims = 3)
end

ts_plot = collect(DateTime(850,2,1):Year(10):DateTime(1996,1,1))

#calculate rolling mean of 10-yr binned records 
ts_plot_r = rolling(mean, yeardecimal.(ts_plot), 5)
mat_br = Array{Float64}(undef, size(mat)[1:2]..., length(ts_plot_r))
for i in 1:length(exp)
    for j in 1:length(inds)
        mat_br[i, j, :] .= rolling(mean, mat_binned[i, j, :], 5)
    end
end

exp_inds = NamedTuple{(:ff, :orb, :volc)}([1:13, 14:16, 17:21])

figure(figsize = (4,8))
for (i, name) in enumerate(indnames)
    subplot(3,1,i)
    title(name)
    @show name
    for (eind, color) in zip(keys(exp_inds), ["black", "blue", "red"])
        @show eind
        plot(ts_plot_r, mat_br[exp_inds[eind], i, :]', color = color, alpha = 0.2)
        ensmean = mean(mat_br[exp_inds[eind], i, :], dims = 1)[:]
        ensvar = var(mat_br[exp_inds[eind], i, :], dims = 1)[:]
        C = Diagonal(ensvar)
        lls, C = linearleastsquares(ts_plot_r, ensmean, C = C)
        slope = (lls[1] ± sqrt(C[1,1])) * 1000
        println(string(slope) * "°C/kyr")
        plot(ts_plot_r, ensmean, color = color, linewidth = 3, label = string(eind) * ": " * string(slope) *  " °C/kyr")
        legend()
        ylabel("Mean Regional SST [°C]")
        yl = gca().get_ylim()
        #text(x = 1950, y = (yl[2] - yl[1]) * 0.8 + yl[1], s = string(slope), color = color)
    end
end
xlabel("Time [years CE]") 
tight_layout()

savefig(plotsdir("NATL_cesmLME_summermean.png"))
#=
[plot(ts_plot, mat_binned[i,1,:], color = "gray") for i in 1:9]
subplot(3,12)
[plot(ts_plot, mat_binned[i,2,:], color = "red") for i in 1:9]
subplot(3,1,3)
[plot(ts_plot, mat_binned[i,3,:], color = "blue", alpha = 0.2) for i in 1:9]
plot(ts_plot, mean(mat_binned[1:9, 3, :], dims = 1)[:], color = "blue")
=#


