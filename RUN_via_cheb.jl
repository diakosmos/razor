"""
This is a mess at the moment. :)
"""

include("viaCheb.jl")
include("classical.jl")
using .bycheb
using LinearAlgebra, Plots
using .Classical
using Plots

method="Explicit"
#method="Implicit"
#method="ImplicitPinv"
#method="Mixed"
μ = 0.5

print("Method: $method\n")

normby = "sum" #true

FACTOR = 0.1

α= 1.0    * FACTOR
ν = 1.0/3   * FACTOR * 9
s = 1.0   * FACTOR * 0.3
Lrel = 5

#N=40
N = floor(Int, 73) #64 / √2) # floor(Int, 73 * √2)

M = 500.0
nMx=0

print("M: $M\n")

TINY = 0.1#1.0#0.1 #0.1

β = 1.0
print("β: $β\n")

#################
rc = Classical.ring(32;xf=4.0)
fig = Plots.plot()

#####

function eye(n::Integer)
    diagm(ones(n))
end

nn = N÷2 + N%2
r = bycheb.ring(N,(α,ν,s,Lrel))
Lcrit = r.Ls["Lcrit"]
Lvisc = r.Ls["Lvisc"]
#Plots.vline!([L0-Lcrit],label="crit",legend=:bottomleft)
Plots.vline!([Lcrit/Lvisc],label="crit",legend=:bottom)#:topleft)
#Plots.vline!([L0-Lvisc],label="visc")
Plots.vline!([Lvisc/Lvisc],label="visc")
#display(fig)
ws = sqrt(3*ν/s) # wave speed
v=[ones(nn);zeros(nn)]
#v[1:nn÷2].=0.0
#Σ = view(v,    1:nn )
#Φ = view(v, nn+1:e)
minAd1 = minimum( abs.(diff(r.y)./(abs.(r.u[1:end-1]).+ws)))
minAd2 = minimum( abs.(diff(r.y)./(abs.(r.u[2:end]).+ws)))
minDiff = minimum( diff(r.y).^2 / (3*ν) )
dt = TINY * 0.5 * min(minAd1, minAd2, minDiff)
print("dt: $dt \n")

IM = eye(2*nn) - dt*(-r.AA+r.B)
MM = eye(2*nn) - μ*dt*(-r.AA+r.B)
#PinvIM = pinv(IM)
Nm = r.Ls["L0"]^5/(5α) / dt# advective time is L0^5/(5α)/dt
NNN = M * Nm
Np = 1002003 #trunc(Int,Nm*0.001)#0

print("Nm: $Nm   NNN: $NNN \n")


rc.β = β
for j in 1:1000000
    Classical.iterate!(rc)
end
Plots.plot!(rc.x_[1:24],rc.σ_[1:24];legend=false)
display(fig)


function normΣ(w,nn,flag)
    if flag=="L1 norm"
        w./=(sum(w[1:nn]) / nn) #not quite right but close?
    elseif flag=="maxnorm"
        w./=maximum(w[1:nn])
    else
        @warn "Error in normΣ.\n"
    end
end
function myplot(r,v;nM="final",first=true)
    ynn=r.y[1:nn]
    Σ = v[1:nn]
    Φ = v[nn+1:end]
    L0 = r.Ls["L0"]
    local Lcrit = r.Ls["Lcrit"]
    local Lvisc = r.Ls["Lvisc"]
    NX = r.parmd["X"] # key dimensionless number

    XXXv = (L0.-ynn)./Lvisc
    XXXc = (L0.-ynn)./Lcrit
    if first
        f = Plots.plot
    else
        f = Plots.plot!
    end
    f(XXXv,Σ,title="X: $NX  α: $α  ν: $ν\ns: $s  N: $N  M: $M  Lrel: $Lrel",label=:none,xlims=(0.0,(L0/Lvisc)))

    #Plots.vline!([L0],label="moon")
end
function myplot!(r,v;nM="final")
    myplot(r,v;nM=nM,first=false)
end

for i in 1:trunc(Int,1000*Np)#NNN)
    global v
    local w
    if method=="Explicit"
        w = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; # w./=sum(w[1:nn]);
        if normby=="max"
            w./=maximum(w[1:nn])
        elseif normby=="sum"
            w./=sum(w[1:nn]/nn);
        elseif normby=="bdy"
            if w[1] <= 0
                w./=sum(w[1:nn]/nn)/nn
            else
                w./= w[1]
                #print("Got here!")
            end
        else
            @warn "Failure\n"
        end#if
        bycheb.setκ!(r,w[1:nn],β)
        v = w + dt*(-r.AA+r.B)*w
    elseif method=="Implicit"
        w = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; # w[1]=0.0; w[nn+1]=0.0;
        if normby=="max"
            w./=maximum(w[1:nn])
        elseif normby=="sum"
            w./=sum(w[1:nn]/nn);
        elseif normby=="bdy"
            if w[1] <= 0
                w./=sum(w[1:nn]/nn)/nn
            else
                w./= w[1]
                #print("Got here!")
            end
        else
            @warn "Failure\n"
        end#if
        bycheb.setκ!(r,w[1:nn],β)
        local IM = eye(2*nn) - dt*(-r.AA+r.B)
        v = IM\w
    elseif method=="ImplicitPinv"
        w = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; w./=sum(w[1:nn]);
        if normbymax
            w./=maximum(w[1:nn])
        else
            w./=sum(w[1:nn]);
        end#if
        bycheb.setκ!(r,w[1:nn],β)
        IM = eye(2*nn) - dt*(-r.AA+r.B)
        local PinvIM = pinv(IM)
        for k in 1:100
            w = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; w./=sum(w[1:nn]);
            if normbymax
                w./=maximum(w[1:nn])
            else
                w./=sum(w[1:nn]);
            end#if
            v = PinvIM * w
        end
    elseif method=="Mixed"
        w  = [abs.(v[1:nn]); -abs.(v[nn+1:end])]; w[1]=0.0; w[nn+1]=0.0; w./=sum(w[1:nn]);
        vv = w + μ*dt*(-r.AA+r.B)*w
        v = MM\vv
    else
        print("WARNING: failure.\n")
    end
    if i%Np==0
        global nMx+=1
        myplot!(r,v;nM="$nMx")
        display(fig)
    end
end
myplot!(r,v)
display(fig)
