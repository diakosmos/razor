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
N = 71

M = 1

TINY = 0.01

#################

function eye(n::Integer)
    diagm(ones(n))
end

nn = N÷2 + N%2
r = bycheb.ring(N,(α,ν,s,Lrel))
ws = sqrt(3*ν/s) # wave speed
v=[ones(nn);zeros(nn)]
#Σ = view(v,    1:nn )
#Φ = view(v, nn+1:e)
minAd1 = minimum( abs.(diff(r.y)./(abs.(r.u[1:end-1]).+ws)))
minAd2 = minimum( abs.(diff(r.y)./(abs.(r.u[2:end]).+ws)))
minDiff = minimum( diff(r.y).^2 / (3*ν) )
dt = TINY * 0.5 * min(minAd1, minAd2, minDiff)

Nm = trunc(Int,1.0/dt)
for i in 1:M*Nm
    global v
    local w
    if method=="Explicit"
        w = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; w./=sum(w[1:nn]);
        #Base.Threads.@threads
        for j in 1:2nn
            v[j] = w[j] + dt*sum((-r.AA[j,:]+r.B[j,:]).*w)
        end
    elseif method=="Implicit"
        w = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; w./=sum(w[1:nn]);
        v = (eye(2*nn) - dt*(-r.AA+r.B))\w
    elseif method=="Mixed"
        w  = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; w./=sum(w[1:nn]);
        vv = w + μ*dt*(-r.AA+r.B)*w
        v = (eye(2*nn) - μ*dt*(-r.AA+r.B))\vv
    end
end
ynn=r.y[1:nn]
Σ = v[1:nn]
Φ = v[nn+1:end];
L0 = r.Ls["L0"]
Lcrit = r.Ls["Lcrit"]
Lvisc = r.Ls["Lvisc"]

Plots.plot(ynn,Σ)
Plots.vline!([L0-Lcrit],label="crit",legend=:bottomleft)
Plots.vline!([L0-Lvisc],label="visc")
Plots.vline!([L0],label="moon")
