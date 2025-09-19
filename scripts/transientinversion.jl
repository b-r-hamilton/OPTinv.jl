#=
Script to generate OPT-3 and OPT-11 estimates 
=# 
import Pkg
Pkg.activate("../")

using OPTinv 
using Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, Revise, DrWatson, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions, PaleoData, DateFormats, Measurements, PythonPlotExt, PythonPlot
import Measurements.value as value
import OPTinv.Est

allcores = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
oldcores = Symbol.("MC" .* string.([26, 22, 13]) .* "A")

# @time allc = invert(inversion(allcores, "svd", "all", "red"), res = 50yr)
# @time oldc = invert(inversion(oldcores, "svd", "old", "blue"), res = 50yr)
@time allc = invert(inversion(oldcores, "global", "all", "#5D3A9B"))
#@time oldc = invert(inversion(oldcores, "svd", "old", "blue"))

#solutions = [oldc, allc]
solutions = [allc]
suffix = "globmode2"#50yr"
