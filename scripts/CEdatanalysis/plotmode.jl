import Pkg;Pkg.activate("../../")
using JLD2, PythonPlot, PythonPlotExt, TMI, DrWatson, OPTinv, PythonCall, NaNMath

corenums_full = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
mode1 = jldopen(DrWatson.datadir("M/svd.jld2"))["SVD"].Vt[1, :]
TMIversion = "modern_180x90x33_GH11_GH12"
γ = Grid(download_ncfile(TMIversion))
m = γreshape(mode1, γ)

lat =γ.lat
lon = γ.lon

m, lon = cu.add_cyclic_point(m', coord = lon)
lon = pyconvert(Array{Float64}, lon)
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]


proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = 1, ncols = 1,
                    subplot_kw = Dict("projection"=>proj), figsize = (10, 10),
                    constrained_layout = true)
ax = axs
ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                 color="black", linewidth=0.5,
                 transform=noproj)
#PyPlot version
#tx_path = ax_hdl[1]._get_transformed_path()
#PythonPlot version
tx_path = first(ax_hdl)._get_transformed_path()

path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
polygon1s = mpath.Path(path_in_data_coords.vertices)
ax.set_boundary(polygon1s)
ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
ax.set_extent([lon_box..., lat_box[1]- 10, lat_box[2]], noproj)
lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
lon_box_ = [360-80 - 20, 30 + 20]

m = pyconvert(Array{Float64}, m)
#m = geospatialsubset(m', lat, lon, lat_box_, lon_box_)
nmin = NaNMath.minimum(m)
nmax = NaNMath.maximum(m)
nrange = nmax - nmin
lev = 8
cf = ax.contourf(lon, lat, m, cmap = "Reds_r", transform = noproj, levels = lev)
c = ax.contour(lon, lat, m, colors = "black", levels = lev, transform = noproj)
ax.clabel(c, fontsize = 15)
savefig(DrWatson.plotsdir("agu_mode1.png"), dpi = 600, transparent = true)
#colorbar(cf)
