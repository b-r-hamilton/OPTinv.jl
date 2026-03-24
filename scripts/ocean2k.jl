#=
McGregor et al 2015 report an
"average anomaly trend of
-0.31°C/ky (median slope)  and -0.32°C/ky (area-weighted average anomaly)
from 1-2000"

The average anomaly method consists of
1) averaging each of the 57 SST reconstructions at 200-year resolutioin
2) subtracting the mean of each reconstruction to produce 57 anomaly time series, with units of °C anom.
3) Averaging the 57 anomaly time series to produce a global mean anomaly trend estimate

Here, I use output from the provided MATLAB script to reconstruct these results

We don't get it perfectly, but will assume some differences due to algorithm selection

Then compute for the records in Figure 7 
=#
import Pkg;Pkg.activate("../")
using OPTinv, PaleoData, PythonPlot, SkipNan, Statistics
using XLSX
using LsqFit

using Measurements; import Measurements.value; import Measurements.uncertainty
plotclose("all")
_, _, _, _, mat = loadOcean2kBinned()
t = vec(mat["rangec"]) #bin center

#subset to 801-1800 CE bins 
t = t[begin+1:6]
sst_binned = mat["binnm"][begin+1:6, :]

μsst = [mean(skipnan(sst_binned[:, i])) for i in 1:size(sst_binned)[2]]
sst_b_mr = sst_binned' .- μsst
sst_b_mr_glob = [mean(skipnan(sst_b_mr[:,i])) for i in 1:length(t)]
sst_b_mr_glob_w = [std(skipnan(sst_b_mr[:,i])) for i in 1:length(t)]
figure()
[plot(t, sst_b_mr[i,:], color = "gray") for i in 1:size(sst_binned)[2]]
plot(t, sst_b_mr_glob, color = "black", linewidth = 5)
@. model(x, p) = p[1]*x + p[2]
#weight by inverse variance
fit = curve_fit(model, t, sst_b_mr_glob, 1 ./ (sst_b_mr_glob_w.^2), [-0.5/1000,0.1])
#fit = curve_fit(model, t, sst_b_mr_glob, [-0.5/1000,0.1])
println("Global Ocean2k Mean Trend in °C/kyr [1000-2000]") 
@show coef(fit)[1] .* 1000
@show stderror(fit)[1] .* 1000
@show confint(fit, level = 0.95)[1] .* 1000


_, oc2k, = loadOcean2k()
# load in the header of the Excel sheet - this gives us core names 
n = names(oc2k)
# get the names of all records, in order
name_ind = unique([x[begin:end-4] for x in n])
recs_of_interest = ["Arctic1147Bonnet2010_JM-06-WP-04-MCB",
                    "Arctic1148Spielhagen2011_MSM5/5-712",
                    "Atlantic0235Sicre2011__RAPiD-21-3K",
                    "Atlantic1565Calvo2002_ MD95-2011",
                    "Atlantic0234Sicre2011__MD99-2275",
                    "Atlantic0058Richter2009__ENAM9606,M200309",
                    "Atlantic0219Came2007__ODP-162-984"] 
nnatl_sub =[findall(x->x==r, name_ind)[1] for r  in recs_of_interest]
sst_b_mr_natl = [mean(skipnan(sst_b_mr[nnatl_sub,i])) for i in 1:length(t)]
sst_b_mr_natl_w = [std(skipnan(sst_b_mr[nnatl_sub,i])) for i in 1:length(t)]
[plot(t, sst_b_mr[i,:], color = "blue") for i in nnatl_sub]
plot(t, sst_b_mr_natl, color = "blue", linewidth = 5)

fit = curve_fit(model, t, sst_b_mr_natl,1 ./ (sst_b_mr_natl_w.^2), [-0.5/1000,0.1])
println("Northern North Atlantic Ocean2k Mean Trend in °C/kyr [1000-2000]") 
#fit = curve_fit(model, t, sst_b_mr_natl, [-0.5/1000,0.1])
@show coef(fit)[1] .* 1000
@show stderror(fit)[1] .* 1000
@show confint(fit, level = 0.95)[1] .* 1000

xlabel("Time [yr CE]", fontsize = 15)
ylabel("Temperature Anomaly [°C]", fontsize = 15)
tight_layout()
savefig("../plots/oceank2_vals.png", dpi = 600)
