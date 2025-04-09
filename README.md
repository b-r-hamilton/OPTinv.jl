# OPTinv
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://b-r-hamilton.github.io/OPTinv.jl/dev/)
[![Build Status](https://github.com/b-r-hamilton/OPTinv.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/b-r-hamilton/OPTinv.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/b-r-hamilton/OPTinv.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/b-r-hamilton/OPTinv.jl)
# What is this?
The source code for `OPTinv` holds methods to compute spatiotemporal surface maps of variables that explain the time evolution of some observed variable from benthic sediment core compilations. The scripts developed here focusing on applying this methodology to infer surface values of temperature and $\mathrm{\delta}^{18}\mathrm{O}\_{\mathrm{seawater}}$ from the [Lu et al, 2023](https://www.science.org/doi/10.1126/science.adf1646) dataset of 11 sediment core records of $\mathrm{\delta}^{18}\mathrm{O}\_{\mathrm{calcite}}$ from the Reykjanes Ridge in Iceland. 

# How to run 
(assuming in terminal) 
0. Get [Julia](https://julialang.org/downloads/) 
1. Navigate to this directory (`OPTinv.jl/scripts`)
2. Start a Julia session in terminal by typing `julia` 
3. Activate the Julia environment in `OPTinv.jl` by typing `]` to enter Package Manager mode, then type `activate ../`. The first time, you will have to instantiate (download) all supporting code, so type `instantiate`. This will take a bit to run. 
4. Exit package manager by pressing back space 
5. Run a script by typing `include("transientinversion.jl")` 

Note: there is a small issue with downloading the 2 degree TMI file used in this project. To circumvent this, download `TMI_modern_180x90x33_GH11_GH12` from [Google Drive](https://drive.google.com/drive/folders/1nFDSINbst84pK68aWwRGBVfYZkfN1mUR), then paste it into the data folder in the TMI package in your package manager system (e.g., on my Unix system, this is at `/home/brynn/.julia/packages/TMI/[version name]/data`.

# Scripts 
Transient Inversion of Lu et al, 2023 records 
- `transientinversion.jl`: main script that generates the 2 estimates (OPT-3 and OPT-11) detailed in the manuscript
- `transientinversionglobal.jl`: same as `transientinversion.jl` but generates the OPT-3-GLOBAL and OPT-3-MODE2 solutions in the manuscript (Figure 9, Figure 10) 

The first time that either `transientinversion.jl` or `transientinversionglobal.jl` are run on a machine, it will need to calculate the SVD modes, propagate them through the transient TMI matrix, and also compute the uncertainties of these modes in CESM. This will take 2+ hours, but the output will be saved for faster (~100 second) runtimes from there on out. 

The following scripts either need to be run after `transientinversion.jl` or `transientinversionglobal.jl`. If you attempt to run with no `solutions` already calculated, it will run `transientinversion.jl`

- `gamma.jl`: stand-alone script that explores sensitivity of theoretical North Atlantic temperature change to the correlation constant (Supplemental Figure S9) 
- `offsets.jl`: stand-alone script that generates depth-profile of effective modern-day d18Oc and compares to the full range of CE variability in recorded sediment core d18Oc (Figure 1) 
- `modes.jl`: stand-alone script that plots surface spatial mode plots, impulse response of each mode at 3 select cores, and singular values w.r.t. sediment core/depth (Figure 2, Supplemental Figure S1, Supplemental Figure S2) 
- `utilde.jl`: follow-on script, plots the mode solutions to `ex3` (Supplemental Figure S6, Figure 3)
- `ytilde.jl`: follow-on script, plots the reconstructions at the core sites for each solution produced in `ex3` (Figure 4) 
- `regionmean.jl`: follow-on script, plots the regional means (for N. Atl. region and subregions) for solutions in `ex3`. Compares to Ocean2k and LMR and HadISST (Figure 5)
- `maps.jl`: follow-on script, plots the surface reconstruction maps for temperature (Figure 6, 7) 
- `effd18Oc.jl`: follow-on script, computes the effective d18Oc at the surface and compares it to the planktic stack from the EN539 sediment cores (Figure 8) 

Analysis of Common Era Data Products
- `CEdatanalysis/NATL_iCESM.jl`: two-part script, first part is code to copy and paste into Casper to download a timeseries of T, S, and d18O from a iCESM and CESM realization. Second part makes Figure B1 and Supplemental Figure S4 
- `CEdatanalysis/T_v_d18O.jl`: generates Supplemental Figure S5
- `CEdatanalysis/sigma.jl`: computes the variance of calculated spatial modes in a CESM historical realization and saves them to a JLD2 file. Will automatically run the first time, shouldn't need to run after that. Also makes Supplemental Figure S3. 
- 
# Philosophy, and important dependencies
The code developed here attempts to preserve information about the dimensionality and units of all variables. This is accomplished by using the following Julia libraries 
- [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) to manage quantities, 
- [`DimensionalData`](https://github.com/rafaqz/DimensionalData.jl) for managing dimensioned data, 
- [`BLUEs.jl`](https://github.com/ggebbie/BLUEs.jl) for solving the linear underdetermined system presented here, 
- [`UnitfulLinearAlgebra.jl`](https://github.com/ggebbie/UnitfulLinearAlgebra.jl) for properly propagating physical quantities through linear algebra 

Information about the modern-day circulation is obtained from the Total Matrix Intercomparison (TMI) results from [Gebbie, 2010](https://journals.ametsoc.org/view/journals/phoc/40/8/2010jpo4272.1.xml) and [Gebbie, 2012](https://journals.ametsoc.org/view/journals/phoc/42/2/jpo-d-11-043.1.xml). These results are accessed through the following Julia libraries 
- [`TMI.jl`](https://github.com/ggebbie/TMI.jl)
- [`TMItransient.jl`](https://github.com/ggebbie/TMItransient.jl) 

# Future Structural Work (i.e., how could this code be better?) 
One of the cruxes of the problem presented here is the dimensionality of the output: we are inferring the temporal evolution of the magnitude of 11 distinct spatial modes for two different state variables, temperature and $\mathrm{\delta}^{18}\mathrm{O}\_{\mathrm{seawater}}$. This is a three-dimensional output that must be rearranged into a 1D vector. We also have a complete representation of the first and second moments of this variable, contained in a covariance matrix. Here, we use the `Est` and `CovMat` structures to manage these high-dimensional data structures. 

The `Est` structure contains a dimensional array `y`, and associated `CovMat`, which is a structure that contains a matrix, and its associated `ax`. The `ax` is a list of Tuples, where each Tuple is generated according to the column-major vectorization of the dimensions associated with the `y`. An `ax` can easily be generated using the `covariancedims` method. 

`CovMat`s can be populated using `fillcovariance` method, which relies upon the generation of `DiagRule` and `OffDiagRule`. This allows for a `CovMat` to be populated with certain values, that can be a constant, at particular dimensioned values, or functions of dimensioned values. 

Jake has attempted to wrap similar, experimental, functionality into `BLUEs.jl` on the [this branch](https://github.com/ggebbie/BLUEs.jl/tree/multipliable-dimarrays) (has not been incorporated into the main code, needs testing). Future work would include testing that methodology, and possibly incorporating elements of what I have developed here, into that branch. 


# Started package by doing
```
t = Template(; user = "b-r-hamilton", dir = "~/Code", authors = "B. Hamilton", julia = v"1.10", plugins = [License(; name = "MIT"), Git(; manifest = true, ssh=true), GitHubActions(; x86=false), Codecov(), Develop(), ],)
t("OPTinv.jl") 
```

