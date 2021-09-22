include("sand_02.jl")
using .Sand_02
using Plots
fig = Plots.plot()


NN = 100 # no of cells
tmax = 40.0 # max time
dt = 4.0 # how often to plot
tp = 0.0#5000 # when to start plotting


r=Sand_02.ring(NN; U=0.5)
for t in 0.0:dt:tmax
    if t > tp
        t = r.t
        x0_ = 0.5 * (r.X_[1:end-1] + r.X_[2:end])
        Plots.plot!(x0_[1:end],r.Σ_[1:end],label="t: $t",legend=false)#:bottomright)#,title="Σ")#;legend=false)
        display(fig)
    end
    Sand_02.tstep!(r,dt)
end
