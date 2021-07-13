include("classical.jl")
using .Classical
using Plots

r = Classical.ring(100)
fig = Plots.plot()

for i in 1:100
    Plots.plot!(r.x_,r.σ_)
    for j in 1:10000
        Classical.diffuse!(r)
    end
    display(fig)
end
