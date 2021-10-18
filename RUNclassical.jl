include("classical.jl")
using .Classical
using Plots

r = Classical.ring(64;xf=4.0)
fig = Plots.plot()

for b in [0.0,1.0,2.0]#,-0.5,1.0,2.0,3.0]
    r.β = b
    for j in 1:1000000
        Classical.iterate!(r)
    end
    Plots.plot!(r.x_[1:24],r.σ_[1:24];legend=false)
    display(fig)
end
