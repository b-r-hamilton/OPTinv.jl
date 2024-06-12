import PyPlot.plot
import PyPlot.scatter

function plot(y::DimArray; kwargs...)
    unc = ustrip.(uncertainty.(vec(y)))    
    x = ustrip.(collect(first(dims(y))))
    xunit = unit(first(first(dims(y))))
    yunit = unit(first(y))
    y = ustrip.(value.(vec(y)))
    if sum(unc) != 0 
        f = fill_between(x = x, y1 = y .- unc, y2 = y .+ unc, alpha = 0.5; kwargs...)
    end
    
    plot(x, y; kwargs...)
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


