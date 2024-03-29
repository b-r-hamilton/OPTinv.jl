#=
for geoChronR output files
- interpolate to given resolution
- compute mean across all resolutions at all interpolated time values
- compute Cnn matrix
- format nicely into a BLUEs.DimEstimate (contains mean, Cnn, and time in one obj) 
=#
using DrWatson
@quickactivate "OPTinv"
using OPTinv, Unitful, DimensionalData
import OPTinv.yr, OPTinv.permil

#for one record, compute 
age_file = datadir("MC26A.U.p.ae.csv")
d18O_file = datadir("MC26A.U.p.d18O.csv")

#each individual step to process a sed. core record and Bacon age model 
a, d = readbaconfile(age_file, d18O_file) 
T, int = interpolatebacon(a, d, 5yr);
μ, T = meanbacon(int, T, a)
println("Time to compute full covariance matrix for MC26A, 5yr res.") 
@time Cnn = covariancebacon(int, T,μ, a, :MC26A)
e = formatbacon(μ, T, Cnn, :MC26A)

#alternatively, format bacon can budle these all up 
println("Time to do all steps for MC26A, 10yr res.") 
@time de = formatbacon(age_file, d18O_file, :MC26A; res = 10yr)

#read multiple files
age_files = datadir.(["MC26A.U.p.ae.csv", "MC9A.C.w.ae.csv"])
d18O_files = datadir.(["MC26A.U.p.d18O.csv", "MC9A.C.w.d18O.csv"])
cores = [:MC26A, :MC9A]
println("Time to do all steps for MC26A and MC9A, 10yr res.") 
@time e = formatbacon(age_files, d18O_files, cores, res = 10yr) 

