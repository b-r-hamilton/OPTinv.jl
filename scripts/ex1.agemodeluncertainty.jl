#=
for geoChronR output files
- interpolate to given resolution
- compute mean across all resolutions at all interpolated time values
- compute Cnn matrix
- format nicely into a BLUEs.DimEstimate (contains mean, Cnn, and time in one obj) 
=#
using Pkg
Pkg.activate("../")
using OPTinv, Unitful, DimensionalData, DrWatson, LinearAlgebra, PyPlot
import OPTinv.yr, OPTinv.permil

#for one record, compute 
age_file = datadir("MC26A.U.p.ae.csv")
d18O_file = datadir("MC26A.U.p.d18O.csv")

#each individual step to process a sed. core record and Bacon age model 
a, d = readbacon(age_file, d18O_file) 
T, int = interpolatebacon(a, d, 5yr);
μ, T = meanbacon(int, T, a)
println("Time to compute full covariance matrix for MC26A, 5yr res.") 
@time Cnn = covariancebacon(int, T,μ, a, :MC26A)
#OH YOU ARE THE SILLIEST GOOSE IN THE WHOLE WORLD 
Cnn.mat[100:200, 100:200]
Cnn2 = covariancebacon(int, T[100:200], μ[100:200], a, :MC26A)
Cnn2.mat






e = formatbacon(μ, T, Cnn, :MC26A)
isposdef(e.Cnn.mat)
numneg = sum(eigvals(e.Cnn.mat) .< 0)
println("Number of negative eigenvalues: " * string(numneg))

#alternatively, format bacon can budle these all up 
println("Time to do all steps for MC26A, 10yr res.") 
@time e = formatbacon(age_file, d18O_file, :MC26A; res = 10yr)
numneg = sum(eigvals(e.Cnn.mat) .< 0)
println("Number of negative eigenvalues: " * string(numneg))

#read multiple files
age_files = datadir.(["MC26A.U.p.ae.csv", "MC9A.C.w.ae.csv"])
d18O_files = datadir.(["MC26A.U.p.d18O.csv", "MC9A.C.w.d18O.csv"])
cores = [:MC26A, :MC9A]
println("Time to do all steps for MC26A and MC9A, 10yr res.") 
@time e = formatbacon(age_files, d18O_files, cores, res = 10yr)
isposdef(e.Cnn.mat)
numneg = sum(eigvals(e.Cnn.mat) .< 0)
println("Number of negative eigenvalues: " * string(numneg))
e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))
isposdef(e.Cnn.mat)
numneg = sum(eigvals(e.Cnn.mat) .< 0)
println("Number of negative eigenvalues: " * string(numneg))


corenums_full = [28, 26, 25, 22, 21, 20, 19,10, 9,13,14]
core_list_full = Symbol.("MC".* string.(corenums_full) .* "A")

y = loadcores(core_list_full)
eigvals(y.Cnn.mat)
figure();plot(eigvecs(parent(y.Cnn.mat))[:, 1])
mc26a = loadcores([:MC26A])
y.y[At(1900.0yr), At(:MC26A)] == mc26a.y[At(1900.0yr), At(:MC26A)]
ind1 = (1800.0yr, :MC26A)
ind2 = (1900.0yr,:MC26A)

#clearly there's something fucked up going on here 
y.Cnn.mat[findall(x->x == ind1, y.Cnn.ax), findall(x->x == ind2, y.Cnn.ax)]
mc26a.Cnn.mat[findall(x->x == ind1, mc26a.Cnn.ax), findall(x->x == ind2, mc26a.Cnn.ax)]

corenums_subset = [13, 28, 26, 22, 20]
isposdef(y.Cnn.mat)

core_list_subset = Symbol.("MC".* string.(corenums_subset) .* "A")
y_subset = loadcores(core_list_subset)
eigvals(y_subset.Cnn.mat)
isposdef(y_subset.Cnn.mat)
figure();plot(eigvecs(parent(y_subset.Cnn.mat))[:, 4])


startyear = []
stopyear = []
for c in core_list_full 
    y = loadcores([c])
    println()
    push!(startyear, y.y.dims[1][1])
    push!(stopyear, y.y.dims[1][end])
    numneg = sum(eigvals(y.Cnn.mat) .< 0)
    println("Number of negative eigenvalues: " * string(numneg))
end

e1 = eigvals(loadcores([:MC10A]).Cnn.mat)
e1_ = eigvals(loadcores([:MC25A]).Cnn.mat)
e2 = eigvals(loadcores([:MC10A, :MC25A]).Cnn.mat)
figure();scatter(1:length(e1), e1)
scatter(1:length(e2), e2)

@show [e ∈ e2 for e in e1]
@show [e ∈ e2 for e in e1_]

figure()
for (i, c) in enumerate([:MC10A, :MC25A])
    evecs = eigvecs(parent(loadcores([c]).Cnn.mat))
    evals = eigvals(loadcores([c]).Cnn.mat)
    subplot(1,2,i) 
    for i in length(evals)-3:length(evals)
        plot(parent(evecs)[:, i], label = evals[i])
    end
    legend()
end

figure()
c = [:MC10A, :MC25A]
evecs = eigvecs(parent(loadcores(c).Cnn.mat))
evals = eigvals(loadcores(c).Cnn.mat)
for i in length(evals)-20:length(evals)
    plot(parent(evecs)[:, i], label = evals[i])
end
legend()

core_list_full[3], startyear[3], stopyear[3]
core_list_full[8], startyear[8], stopyear[8]
#1080 - 1980
#so determinants of MC10A should stay the same? 

mat = parent(loadcores([:MC25A]).Cnn.mat)
timeax = [i[1] for i in loadcores([:MC25A]).Cnn.ax]
timeax_match =[i[1] for i in loadcores(c).Cnn.ax]
ind_match = [t ∈ timeax_match for t in timeax]

evals = eigvals(mat[ind_match, ind_match])
evecs = eigvecs(parent(mat[ind_match, ind_match]))
figure()
for i in length(evals)-3:length(evals)
    plot(evecs[:, i], label = evals[i])
end
legend()

loadcores(c).Cnn.mat[92:end, 92:end]

loadcores([:MC25A]).Cnn.mat[ind_match, ind_match]
