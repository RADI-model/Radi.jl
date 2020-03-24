module tst

using Plots
Plots.default(show=true)
Plots.closeall()

include("RADI.jl")
import .RADI

depths, oxy, poc = RADI.model(5.0, 1/128000, 128000)

p1 = plot(depths.*100, oxy.*1e3, legend=false)
p2 = plot(depths.*100, poc, legend=false)
plot(p1, p2, layout=(2, 1))

RADI.say_RADI()

# # Trying out composite type but slower and too clever for its own good
# struct PorewaterVariable
#     depths::Array{Float64,1}
#     start::Array{Float64,1}
#     previous::Array{Float64,1}
#     now::Array{Float64,1}
#     npts::Int64
#     overlying::Float64
# end # struct PorewaterVariable
#
# function PorewaterVariable(depths::Array{Float64,1}, start::Array{Float64,1},
#         overlying::Float64)
#     npts = length(depths)
#     previous = fill(NaN, npts)
#     PorewaterVariable(depths, start, previous, start, npts, overlying)
# end # function PorewaterVariable
#
# function PorewaterVariable(depths::Array{Float64,1}, start::Float64,
#         overlying::Float64)
#     npts = length(depths)
#     start = fill(start, npts)
#     previous = fill(NaN, npts)
#     PorewaterVariable(depths, start, previous, start, npts, overlying)
# end # function PorewaterVariable
#
# function now2previous!(myh::PorewaterVariable)
#     myh.previous[:] = myh.now
# end # function now2previous!

end # module tst
