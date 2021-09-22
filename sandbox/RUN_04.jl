include("sand_04.jl")
using .Sand_04
using Plots
fig = Plots.plot()


NN = 20 # no of cells
tmax = 0.4 # max time
dt = 0.4 # how often to plot
tp = 0.0#5000 # when to start plotting


r= Sand_04.ring(NN; U=5.5)
for t in 0.0:dt:tmax
#for i in 1:10
    t = r.t
    if t ≥ tp
        x0_ = 0.5 * (r.X_[1:end-1] + r.X_[2:end])
#        Plots.plot!(x0_[1:end],r.Σ_[1:end],label="t: $t",legend=false)#:bottomright)#,title="Σ")#;legend=false)
#        display(fig)
        Plots.plot!(x0_[1:end],r.J_[1:end],label=false)#,label="t: $t",legend=:topleft)#,title="J")#;legend=false)
        display(fig)
    end
    Sand_04.tstep!(r,dt)
    #Sand_03.step!(r)
end
