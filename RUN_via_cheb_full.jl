include("viaCheb.jl")
using .bycheb
using LinearAlgebra, Plots

#method="Explicit"
#ethod="Implicit"
method="Mixed"; μ = 1.0


α= 1.0
ν = 1.0/3
s = 0.01
Lrel = 10

#N=40
N = 23

M = 1

TINY = 0.001

#################

function eye(n::Integer)
    diagm(ones(n))
end

r = bycheb.ring(N;params=(α,ν,s),Lrel=Lrel)

v=[ones(r.nn);zeros(r.nn)]; v[r.nn+1]=-1

dt = min(0.001,r.parmd["srel"])
#dt = 0.5 * minimum( abs.(diff(r.y)./(abs.(r.u[1:end-1]).+ws)))
#Nm = trunc(Int,1.0/(dt*TINY))
fig=Plots.plot()
for i in 1:100#Nm÷10 #M*Nm
    global v
    local w 
    if method=="Explicit"
        #w = [abs.(v[1:N]); -abs.(v[N+1:end])]; w[1]=1.0; w[N]=1.0;
        #w[N+1]=-1.0; w[end]=1.0; w./=sum(w[1:N]);
        #w[1:N] = 0.5*(w[1:N]+w[N:-1:1])
        #w[N+1:2N] = 0.5 * (w[N+1:2N]-w[2N:-1:N+1])
        #v = w + dt*(-r.AAf+r.Bf)*w
    elseif method=="Implicit"
        #w = [abs.(v[1:N]); -abs.(v[N+1:end])]; w[1]=0.0; w[N]=0.0;
        #w[N+1]=0.0; w[end]=0.0; w./=sum(w[1:nn]);
        #v = (eye(2*N) - dt*(-r.AA+r.B)) \ w
    elseif method=="Mixed"
        w = ( eye(2r.nn) + (1-μ)*dt*(r.AA+r.B) )*v
        w[1:r.nn] .*= (w[1:r.nn].>0)
        w[1:r.nn]./=w[1]; w[r.nn+1:end]./=-w[r.nn+1]
        v = (eye(2r.nn) - μ*dt*(r.AA+r.B)) \ w
        v[1:r.nn] .*= (v[1:r.nn].>0)
        v[1:r.nn]./=v[1]; v[r.nn+1:end]./=-v[r.nn+1]
        fig=Plots.plot!(r.y[1:r.nn],v[1:r.nn];label=:none)
    end
end
display(fig)
#ynn=r.y[1:nn]
Σ_ = v[1:r.nn]
Φ_ = v[r.nn+1:end];

Lcrit = r.Ls["Lcrit"]
Lvisc = r.Ls["Lvisc"]


Plots.plot(r.y[1:r.nn],Σ_)
#Plots.vline!([yc,-yc],label="crit",legend=:bottomleft)
#Plots.vline!([yv,-yv],label="visc")
#Plots.vline!([L0,-L0],label="moon")
