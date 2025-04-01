#=
This script consists of two parts
(1) a script that can be copy and pasted into a Julia REPL on Casper to
download a timeseries of T, S, δ¹⁸O from a CESM and iCESM realization
(2) an analysis script that will make Supplemental Figure B1 
=#

# ============== CODE TO RUN ON CASPER ==================== #
#=
using NCDatasets, JLD2, Statistics, Dates

const DTM = Union{Date, DateTime}
# assuming you want 2000-01-01 to become 2000.0
yfrac(dtm::DTM) = (dayofyear(dtm) - 1) / daysinyear(dtm)
decimaldate(dtm::DTM) = year(dtm) + yfrac(dtm)

yearstringsstart = vcat(["085001", "090001"], string.(1000:100:1800) .* "01" )
yearstringsend = vcat(["089912", "099912"], string.(1099:100:1799) .* "12", ["184912"])

compset = ["b.e11.BLMTRC5CN.f19_g16.006.pop.h.", "b.ie12.B1850C5CN.f19_g16.LME.002.pop.h."]
varnames = [["SALT", "TEMP"], ["SALT", "TEMP", "R18O"]]

for (i, c) in enumerate(compset)
    for var in varnames[i]
        path = "/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/ocn/proc/tseries/monthly/"
        fnames = path * var * "/" .* [c * var * "." * yearstringsstart[j] * "-" * yearstringsend[j] * ".nc" for j in 1:length(yearstringsend)]
        d = Dict()
        for f in fnames
            @show f
            nc = NCDataset(f)
            d[f] = nc[var][21, 355, 45, :]
            #maybe? 
            d[f*"_tlong"] = nc["TLONG"][355, 21]
            d[f*"_tlat"] = nc["TLAT"][355,21]
        end
        savename = c * "_" * var * "_EN539.jld2"
        #just save in SALT dir, otherwise change var name 
        jldsave("/glade/u/home/glennliu/scratch/CESM_LME_Crop/SALT/" * savename; d) 
    end
end
=#

# ================== CODE TO GENERATE FIGURE B1 and SUPPLEMENTAL FIGURE 4 =========== #
import Pkg;Pkg.activate("../../")
using NCDatasets, JLD2, Statistics, Dates

const DTM = Union{Date, DateTime}
# assuming you want 2000-01-01 to become 2000.0
yfrac(dtm::DTM) = (dayofyear(dtm) - 1) / daysinyear(dtm)
decimaldate(dtm::DTM) = year(dtm) + yfrac(dtm)

yearstringsstart = vcat(["085001", "090001"], string.(1000:100:1800) .* "01" )
yearstringsend = vcat(["089912", "099912"], string.(1099:100:1799) .* "12", ["184912"])

compset = ["b.e11.BLMTRC5CN.f19_g16.006.pop.h.", "b.ie12.B1850C5CN.f19_g16.LME.002.pop.h."]
varnames = [["SALT", "TEMP"], ["SALT", "TEMP", "R18O"]]

using PythonPlot, DrWatson, Dates, LinearAlgebra, PythonCall

function linearleastsquares(x,y)
    E = ones(length(x), 2)
    E[:, 1] .= x
    F = (E'*E)\I*E'
    lls = F*y
    return lls #slope, intercept 
end

datapath = ("../../data/iCESM")
files = readdir(datapath)
t = collect(Date(850,2,1):Month(1):Date(1850, 1, 1))

tind = findall(x->x>=850, year.(t))

compsetnames = ["CESM1 LME Realization 006", "iCESM1.2"]
titles = [["CESM1 LME SALT", "CESM1 LME TEMP"], ["iCESM1.2 SALT", "iCESM1.2 TEMP", "iCESM1.2 δ¹⁸O"]]
ylabels = [["S [g/kg]","T [°C]"], ["S [g/kg]", "T [°C]", "δ¹⁸O [‰]"]]
vars = [["S", "T"], ["S", "T", "δ¹⁸O"]]
sd = [[4,2],[4,2,2]]
letters = [["A", "B"], ["C", "D", "E"]]
num = 5*12
figure(figsize = (13,7))
for (i, c) in enumerate(compset)
    for (j, var) in enumerate(varnames[i])
        savename = c * "_" * var * "_EN539.jld2"
        d = jldopen(joinpath(datapath, savename))["d"]
        path = "/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/ocn/proc/tseries/monthly/"
        fnames = path * var * "/" .* [c * var * "." * yearstringsstart[j] * "-" * yearstringsend[j] * ".nc" for j in 1:length(yearstringsend)]
        data = vcat([d[fnames[k]] for k in 1:length(fnames)]...)
        subplot(2,3, (i - 1) * 3 + j)
        if i == 2 && j == 3 #R18O
            #display(data)
            x = convert(Vector{Float64}, decimaldate.(t)) 
            y = convert(Vector{Float64}, data .- 1) .* 1000
            plot(x, y, color = "gray", linewidth = 1)
            x = [mean(x[i:i+num]) for i in collect(1:num:length(x)-num)][20:end]
            y = [mean(y[i:i+num]) for i in collect(1:num:length(y)-num)][20:end]
            plot(x, y, color = "black", linewidth = 3)
            yt = round.(pyconvert(Vector{Float64}, gca().get_yticks()), sigdigits = sd[i][j])
            yticks(yt, yt, fontsize = 12)

        else
            x = convert(Vector{Float64}, decimaldate.(t))
            y = convert(Vector{Float64}, data)
            plot(x, y, color = "gray", linewidth = 1) 
            x = [mean(x[i:i+num]) for i in collect(1:num:length(x)-num)][20:end]
            y = [mean(y[i:i+num]) for i in collect(1:num:length(y)-num)][20:end]
            plot(x, y, color = "black", linewidth =3)
            yt = round.(pyconvert(Vector{Float64}, gca().get_yticks()), sigdigits = sd[i][j])
            yticks(yt, yt, fontsize = 12)
        end
        title(titles[i][j], fontsize = 20, weight = "bold")
        xlabel("Time [yr CE]", fontsize = 15)
        ylabel(ylabels[i][j], fontsize = 15)
        xticks(1000:200:1800, fontsize = 12)
        xl = gca().get_xlim()
        yl = gca().get_ylim()
        text(x = xl[0] + (xl[1] - xl[0])*0.88, y = yl[0] + (yl[1] - yl[0])*0.85, s = letters[i][j], fontweight = "bold", fontsize = 30)
    end
end
tight_layout()
savefig(plotsdir("cesm_timeseries.png"), dpi = 600)
#x = datadict[1]

function plotme(x,y,ax, xlab,ylab,t,varx, vary, letter, cb = true)
    x = convert(Vector{Float32}, x[tind])
    y = convert(Vector{Float32}, y[tind])
    num = 5*12
    x = [mean(x[i:i+num]) for i in collect(1:num:length(x)-num)][20:end]
    y = [mean(y[i:i+num]) for i in collect(1:num:length(y)-num)][20:end]
    
    s = scatter(x, y, c = "black")    
    ylabel(ylab, fontsize = 12)
    xlabel(xlab, fontsize = 12)

    lls = linearleastsquares(x,y)
    ỹ = lls[1] .* x .+ lls[2]


    xt = round.(pyconvert(Vector{Float64}, gca().get_xticks()), sigdigits = 5)
    yt = round.(pyconvert(Vector{Float64}, gca().get_yticks()), sigdigits = 5)
    xl = gca().get_xlim()
    yl = gca().get_ylim() 
    plot(x, ỹ, color = "red", label = "slope = " * string(round(lls[1], sigdigits = 3)) * "°C/g/kg", linewidth = 5)

    xticks(xt, xt, fontsize = 12)
    yticks(yt, yt, fontsize = 12)
    gca().set_xlim(xl)
    gca().set_ylim(yl)
    
    r2 = round(1 - sum((y .- ỹ).^2) /sum((y .- mean(y)).^2), sigdigits = 3)
    α = round(lls[1], sigdigits = 2)
    β = round(lls[2], sigdigits = 2) 
    text(x = xl[0] + (xl[1] - xl[0])*0.03, y = yl[0] + (yl[1] - yl[0])*0.85, s = vary * " = " * string(α) * varx * " + " * string(β) * "\nR² = " * string(r2), fontweight = "bold", fontsize = 12)
    text(x = xl[0] + (xl[1] - xl[0])*0.85, y = yl[0] + (yl[1] - yl[0])*0.05, s = letter, fontweight = "bold", fontsize = 30)
end

figure(figsize = (9, 4))
#for (i, c) in enumerate(compset)
i = 2
c = compset[2] 
    datadict = Dict() 
    for (j, var) in enumerate(varnames[i])
        savename = c * "_" * var * "_EN539.jld2"
        d = jldopen(joinpath(datapath, savename))["d"]
        path = "/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/ocn/proc/tseries/monthly/"
        fnames = path * var * "/" .* [c * var * "." * yearstringsstart[j] * "-" * yearstringsend[j] * ".nc" for j in 1:length(yearstringsend)]
        
        data = vcat([d[fnames[k]] for k in 1:length(fnames)]...)
        
        if i == 2 && j == 3 #R18O
            datadict[j] = (data .- 1) .* 1000 
        else
            datadict[j] = data
        end
        
    end

ax = subplot(1,2,1)
plotme(datadict[1], datadict[2], ax, ylabels[i][1], ylabels[i][2],compsetnames[i], vars[i][1], vars[i][2], "A", false) 
ax = subplot(1,2,2)
plotme(datadict[3], datadict[1], ax,ylabels[i][3], ylabels[i][1],compsetnames[i], vars[i][3], vars[i][1], "B",  false)
xtmp = range(start = 0.265 - 0.0025, stop = 0.283, length = 10)
ytmp022 = xtmp .* 1/0.22 .- xtmp[1] * 1/0.22 .+ 35.11
ytmp06 = xtmp .* 1/0.5 .- xtmp[1] * 1/0.5 .+ 35.11
plot(xtmp, ytmp06, color = "red", linestyle = "dashed")
tight_layout()
savefig(plotsdir("cesm_scatter.png"), dpi = 600) 
