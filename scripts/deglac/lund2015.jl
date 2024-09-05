import Pkg
Pkg.activate("../../")
using OPTinv, PaleoData, LinearAlgebra

locs = loadLund2015()

mat, smat = OPTinv.generatemodes(locs, func = svd)
