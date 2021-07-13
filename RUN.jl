include("classical.jl")
using .Classical
using Plots

r = Classical.ring(100;xf=10.0)
fig = Plots.plot()

for i in 1:10
    Plots.plot!(r.x_,r.σ_;legend=false)
    for j in 1:10000
        Classical.diffuse!(r)
        Classical.advect!(r)
    end
    display(fig)
end
