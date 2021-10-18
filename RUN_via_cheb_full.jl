include("viaCheb.jl")
using .bycheb
using LinearAlgebra, Plots

method="Explicit"
#method="Implicit"
#method="Mixed"; μ = 0.5


α= 1.0
ν = 30.0
s = 30.0
Lrel = 10

#N=40
N = 15

M = 1

TINY = 0.001

#################

function eye(n::Integer)
    diagm(ones(n))
end

r = bycheb.ring(N,(α,ν,s,Lrel))
ws = sqrt(3*ν/s) # wave speed
v=[ones(N);zeros(N)]
#Σ = view(v,    1:nn )
#Φ = view(v, nn+1:e)
minAd1 = minimum( abs.(diff(r.y)./(abs.(r.u[1:end-1]).+ws)))
minAd2 = minimum( abs.(diff(r.y)./(abs.(r.u[2:end]).+ws)))
minDiff = minimum( diff(r.y).^2 / (3*ν) )
dt = TINY * 0.5 * min(minAd1, minAd2, minDiff)
#dt = 0.5 * minimum( abs.(diff(r.y)./(abs.(r.u[1:end-1]).+ws)))
Nm = trunc(Int,1.0/(dt*TINY))
for i in 1:Nm÷10 #M*Nm
    global v
    local w
    if method=="Explicit"
        w = [abs.(v[1:N]); -abs.(v[N+1:end])]; w[1]=0.0; w[N]=0.0;
        w[N+1]=0.0; w[end]=0.0; w./=sum(w[1:N]);
        w[1:N] = 0.5*(w[1:N]+w[N:-1:1])
        w[N+1:2N] = 0.5 * (w[N+1:2N]-w[2N:-1:N+1])
        v = w + dt*(-r.AAf+r.Bf)*w
    elseif method=="Implicit"
        w = [abs.(v[1:N]); -abs.(v[N+1:end])]; w[1]=0.0; w[N]=0.0;
        w[N+1]=0.0; w[end]=0.0; w./=sum(w[1:nn]);
        v = (eye(2*N) - dt*(-r.AA+r.B))\w
    elseif method=="Mixed"
        w  = [abs.(v[1:N]); -abs.(v[N+1:end])]; w[1]=0.0; w[N]=0.0;
        w[N+1]=0.0; w[end]=0.0; w./=sum(w[1:nn]);
        vv = w + μ*dt*(-r.AA+r.B)*w
        v = (eye(2*N) - μ*dt*(-r.AA+r.B))\vv
    end
end
#ynn=r.y[1:nn]
Σ = v[1:N]
Φ = v[N+1:end];
L0 = r.Ls["L0"]
Lcrit = r.Ls["Lcrit"]
Lvisc = r.Ls["Lvisc"]
yc = L0-Lcrit
yv = L0-Lvisc

Plots.plot(r.y,Σ)
Plots.vline!([yc,-yc],label="crit",legend=:bottomleft)
Plots.vline!([yv,-yv],label="visc")
Plots.vline!([L0,-L0],label="moon")
