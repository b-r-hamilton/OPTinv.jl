import Pkg; Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, DrWatson
import Measurements.value
close("all")

files = readdir(datadir())
ae_files = [f for f in files if occursin("ae", f)]
d18O_files = [f for f in files if occursin("d18O", f)]
cores = Tuple([Symbol(split(f, ".")[1]) for f in ae_files])

@time e = formatbacon(datadir.(ae_files), datadir.(d18O_files), [c for c in cores], res = 5yr)
#add in measurement uncertainty 
e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))

#generate steady-state M

surfaceboxes = [[[300, 320], [50,70]], [[320, 340], [50, 70]], [[340, 360],  [50, 70]], [[300, 360], [70, 90]], [[300,360], [30, 50]]]

modes = (:East,:Central, :West, :North, :South)
patches = NamedTuple{modes}(surfaceboxes) 
corelocs = core_locations()
M = steadystateM(corelocs, patches)

#get first guess and first guess covariance 
T = collect(e.y.dims[1])
σθ = 11K /2
σδ = 0.8permil / 1 #std of global surface d18O in WOCE
f = σθ^2/ (20K/permil) 
u₀ = firstguess(T,[m for m in modes], σθ, σδ, f)

#use impulseresponse technique to make E matrix 
predict(u) = u[:, :, At(:θ)] * M * 0.27permil/K + u[:, :, At(:δ)] * M
E = impulseresponse(predict, u₀.y)

@time yde, ỹde, u₀de, ũ = solvesystem(e, u₀, E, predict);

figure(figsize = (8, 8))
for (i, core) in enumerate(cores) 
    subplot(4,3, i)
    plot(yde.x[:, At(core)], label = "y", color = "black")
    plot(ỹde.x[:, At(core)], label = "ỹ", color = "red")
    ylim(-0.25, 0.25)
    title(core)
end
tight_layout()

figure(figsize = (8,4))
for (i, m) in enumerate(modes)
    for (j, s) in enumerate([:θ, :δ])
        subplot(2, length(modes), (j - 1) * length(modes) + i)
        title(string(m) * " " * string(s)) 
        plot(u₀de.x[:, At(m), At(s)], label = "u₀", color = "gray")
        plot(ũ.x[:, At(m), At(s)], label = "ũ", color = "red")    
    end
end
tight_layout()

figure()
scatter(vec(ũ.x[:, :, At(:θ)]), vec(ũ.x[:, :, At(:δ)]), capsize = 5)
lls = linearleastsquares(ustrip.(value.(vec(ũ.x[:, :, At(:δ)]))), ustrip.(value.(vec(ũ.x[:, :, At(:θ)]))))
println("slope [K/permil] = " * string( lls[1]))
