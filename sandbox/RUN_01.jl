include("sand_01.jl")
using .Sand_01
using Plots
fig = Plots.plot()


NN = 200 # no of cells
tmax = 10.0 # max time
dt = 1.0 # how often to plot
tp = 0.0#5000 # when to start plotting


r=Sand_01.ring(NN; U=0.5)
for t in 0.0:dt:tmax
    if t > tp
        t = r.t
        x0_ = 0.5 * (r.X_[1:end-1] + r.X_[2:end])
        Plots.plot!(x0_[1:end],r.Σ_[1:end],label="t: $t",legend=:bottomright)#,title="Σ")#;legend=false)
        display(fig)
    end
    Sand_01.tstep!(r,dt)
end
