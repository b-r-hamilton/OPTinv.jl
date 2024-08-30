
import Pkg;Pkg.activate("../")
using OPTinv, PyPlot, PaleoData, StringDistances, Measurements, Statistics, Unitful, LinearAlgebra, UnitfulLinearAlgebra

oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
oc2k, oc2kdata = loadOcean2k()

oc2k_binned #binned, std values
oc2k_binnedv #number of values in each bin

oc2k #locs
oc2kdata #actual data

#what I want is an estimate of cooling rate, in K/kyr of N. Atl. recs
atl_indices = findall(x->occursin("Atlantic", x), oc2k[!, "name"])
natl_indices = intersect(atl_indices, findall(x->x>0, oc2k[!, "lat"]))


binleft = collect(700:200:1700)
bincenter = collect(800:200:1800)
binright = collect(900:200:1900)

#=
binleft = collect(-100:200:1900)
bincenter = collect(0:200:2000)
binright = collect(100:200:2100)
=#


labels = ["N. Atl.", "Global"]
for (lab, ind) in zip(labels, [natl_indices, collect(1:57)])
    binnedtemps = Matrix{Measurement}(undef, length(ind), length(bincenter))
for (bti, oc2ki) in enumerate(ind)
    ageindex = findmax([compare(oc2k[oc2ki, "name"] * "_Age", n, Levenshtein()) for n in names(oc2kdata)])[2]
    dataindex = findmax([compare(oc2k[oc2ki, "name"] * "_SST", n, Levenshtein()) for n in names(oc2kdata)])[2]
    #@show ageindex, dataindex
    t = oc2kdata[!, ageindex]
    y = oc2kdata[!, dataindex]
    missingind = findall(x->x isa Number, t) 
    t = t[missingind]
    y = y[missingind]
    y .-= mean(y)
    for (j, (bl, br)) in enumerate(zip(binleft, binright))
        ind = findall(x->bl<x<br, t)
        if length(ind) == 0
            binnedtemps[bti,j] = NaN 
        elseif length(ind) == 1
            binnedtemps[bti,j] = Measurements.measurement(mean(y[ind]), 0.0)
        else
            binnedtemps[bti, j] = Measurements.measurement(mean(y[ind]), std(y[ind]))
        end 
    end    
end

#binnedtemps .-= binnedtemps[:, end]
atlmean = Vector{Measurement{Float64}}(undef, length(bincenter))
i = 1
for i in 1:length(bincenter)
    y = binnedtemps[:, i]
    atlmean[i] = mean(y[(!).(isnan.(y))])
end
figure()
    scatter(collect(bincenter), atlmean)

    C = UnitfulMatrix(Diagonal(Measurements.uncertainty.(atlmean).^2), fill(K, length(atlmean)), fill(K^-1, length(atlmean)))
    lls, Cnew = linearleastsquares(UnitfulMatrix(collect(bincenter)yr), UnitfulMatrix(Measurements.value.(atlmean)K), C = C)
    println(lab)
    println(string(Measurements.measurement(lls[1] * 1000,sqrt(Cnew[1,1]) * 1000)) * "°C/kyr")

    slopes = Vector{Measurement{Float64}}(undef, size(binnedtemps)[1])
    for i in 1:size(binnedtemps)[1]
        y = binnedtemps[i, :]
        mind = (!).(isnan.(y))
        if sum(mind) > 1
            C_ = UnitfulMatrix(Diagonal(Measurements.uncertainty.(y[mind]).^2), fill(K, sum(mind)), fill(K^-1, sum(mind)))
        lls_, Cnew_ = linearleastsquares(UnitfulMatrix(collect(bincenter[mind])yr), UnitfulMatrix(Measurements.value.(y[mind])K), C=C_)
            slopes[i] = Measurements.measurement(lls_[1]*1000, sqrt(Cnew_[1,1])*1000)
        else
            slopes[i] = NaN
        end
        
    end
    mind = (!).(isnan.(slopes))
    @show median(slopes[mind])

end




