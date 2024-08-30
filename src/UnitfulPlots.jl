import PyPlot.plot
import PyPlot.scatter

function plot(y::DimArray; lwcentral=3, lwedges=0, kwargs...)
    unc = ustrip.(uncertainty.(vec(y)))    
    x = ustrip.(collect(first(dims(y))))
    xunit = unit(first(first(dims(y))))
    yunit = unit(first(y))
    y = ustrip.(value.(vec(y)))
    if sum(unc) != 0 
        f = fill_between(x = x, y1 = y .- unc, y2 = y .+ unc, alpha = 0.3, linewidth = lwedges; kwargs...)
    end
    
    plot(x, y, alpha = 0.7, linewidth = lwcentral; kwargs...)
    xlabel(string(xunit))
    ylabel(string(yunit))
    return x, y 
end

function scatter(x::Vector{Quantity{Measurement{T}}}, y::Vector{Quantity{Measurement{T}}}; kwargs...) where T <: AbstractFloat
    xerr = uncertainty.(ustrip.(x))
    yerr = uncertainty.(ustrip.(y))
    xunit = unit(first(x))
    yunit = unit(first(y))
    x = value.(ustrip.(x))
    y = value.(ustrip.(y))
    errorbar(x, y, xerr = xerr, yerr = yerr; kwargs...)
    xlabel(string(xunit))
    ylabel(string(yunit))
end

function scatter(x::Vector, y::Vector{Measurement{T}}; kwargs...) where T <: AbstractFloat
    #xerr = uncertainty.(ustrip.(x))
    yerr = uncertainty.(ustrip.(y))
    #xunit = unit(first(x))
    yunit = unit(first(y))
    #x = value.(ustrip.(x))
    y = value.(ustrip.(y))
    errorbar(x, y, yerr = yerr; kwargs...)
    #xlabel(string(xunit))
    #ylabel(string(yunit))
end

function scatter(x::Vector{Measurement{T}}, y::Vector; kwargs...) where T <: AbstractFloat
    
    xerr = uncertainty.(ustrip.(x))
    #yerr = uncertainty.(ustrip.(y))
    #xunit = unit(first(x))
    #yunit = unit(first(y))
    x = value.(ustrip.(x))
    y = value.(ustrip.(y))
    errorbar(x, y, xerr = xerr; kwargs...)
    #xlabel(string(xunit))
    #ylabel(string(yunit))
end


#this one def. works 
function scatter(x::Array{Quantity{Measurement{T}, D1, A1}}, y::Array{Quantity{Measurement{T}, D2, A2}}; kwargs...) where {T <: AbstractFloat, D1, A1, D2, A2}
    xerr = uncertainty.(ustrip.(x))
    yerr = uncertainty.(ustrip.(y))
    xunit = unit(first(x))
    yunit = unit(first(y))
    x = value.(ustrip.(x))
    y = value.(ustrip.(y))
    errorbar(x, y, xerr = xerr, yerr = yerr; kwargs...)
    xlabel(string(xunit))
    ylabel(string(yunit))
end

function scatter(x::Vector{Quantity{T, D1, A1}}, y::Vector{Quantity{T, D2, A2}}; kwargs...) where {T <: AbstractFloat, D1, A1, D2, A2}
    xunit = unit(first(x))
    yunit = unit(first(y))
    x = value.(ustrip.(x))
    y = value.(ustrip.(y))
    scatter(x, y; kwargs...)
    xlabel(string(xunit))
    ylabel(string(yunit))
end

function plot(x::Vector{Union{Missing, T1}}, y::Vector{Union{Missing, T2}}; kwargs...) where {T1, T2}
    x[ismissing.(x)] .= NaN
    y[ismissing.(y)] .= NaN
    x = convert(Vector{T}, x)
    y = convert(Vector{T}, y)
    plot(x, y; kwargs...) 
end

function plot(x::Vector{T1}, y::Vector{Union{Missing, T2}}; kwargs...) where {T1, T2}
    x[ismissing.(x)] .= NaN
    y[ismissing.(y)] .= NaN
    x = convert(Vector{T1}, x)
    y = convert(Vector{T2}, y)
    plot(x, y; kwargs...) 
end

function plot(x::Vector{Union{Missing, T1}}, y::Vector{T2}; kwargs...) where {T1, T2}
    x[ismissing.(x)] .= NaN
    y[ismissing.(y)] .= NaN
    x = convert(Vector{T1}, x)
    y = convert(Vector{T2}, y)
    plot(x, y; kwargs...) 
end

function plot(x::Vector{Quantity{T}}, y::Vector{Quantity{T}}; kwargs...) where T
    xu = unique(unit.(x))
    yu = unique(unit.(y))
    plot(ustrip.(x), ustrip.(y); kwargs...)
    xlabel(xu)
    ylabel(yu)
end

#this doesn't work, not sure why? 
#=
function contourf(x::Vector{Number}, y::Vector{Number}, z::Matrix{Union{Missing, T}}, levels; kwargs...) where T
    z[ismissing.(z)] .= NaN
    z = convert(Matrix{T}, z)
    contourf(x,y,z,levels;kwargs...) 
end

function contour(x, y, z::Matrix{Union{Missing, T}}, levels; kwargs...) where T
    z[ismissing.(z)] .= NaN
    z = convert(Matrix{T}, z)
    contour(x,y,z,levels;kwargs...) 
end
=#





