
module bycheb
using LinearAlgebra
include("cheb.jl")
using .cheb

export ring

EXTEND = 0.0 #0.5

"""
Σ is an even function of x
Φ is an odd function of x
U is an odd function of x
"""

mutable struct ring
    N::Integer              # number of nodes for full ring
    nn::Integer
    parmd::Dict{String,Float64} # fundamental physical parameter
    Ls::Dict{String,Float64}  # dictionary of fundamental lengthscales
    y::Any
    u::Any
    A::Any
    AA::Any
    B::Any
    D1::Any
    De::Any
    Do::Any
    Af::Any
    AAf::Any
    Bf::Any
#    X_::Array{Float64,1}
#    Σ_::Array{Float64,1}
#    Φ_::Array{Float64,1}
end
function ring(N=20, params=(1.0,1.0,1.0,10.0))
    nn = N÷2 + N%2  # nn: "non-negative" (positions)
    # Get key physical parameters, and lengthscales to non-dimensionalize X-axis.
    α = params[1]; ν = params[2]; s = params[3]; Lrel = params[4] # velocity amplitude, kinematic viscosity, relaxation time, Lrel = L/Lc
    q = 3*ν/s; κ = 3*ν
    NX = α^2/κ^5/s^3
    parmd = Dict("α"=>α, "ν"=>ν, "s"=>s, "q"=>q, "κ"=>κ, "Lrel"=>Lrel, "X"=>NX)
    pa = sort(collect(parmd),by=x->x[2])
    for p in pa
        print("$p\n")
    end
    #
    Lc = (α/√q)^0.25        # position of critical point for hyperbolic transport (could call it "L0" I suppose)
    Lv = (α/κ)^(1.0/3.0) # viscous "critical" point
    Lαs = (α*s)^0.2          # another lengthscale constructed from α, ν, s
    Lνs = √(ν*s)             # another lengthscale      "        "    "   "
    L0 = Lrel*max(Lc,Lv)
    L = L0-Lc
    Ls = Dict("Lcrit"=>Lc, "Lvisc"=>Lv, "Lαs"=>Lαs, "Lνs"=>Lνs, "L"=>L, "L0"=>L0)
    a = sort(collect(Ls),by=x->x[2]) # this is a bit ugly... ugh
    for e in a
        print("$e\n")
    end
    #
    (x,DM,DeM,DoM) = chebdif(N,1)
    D1 = DM[:,:,1]   # full derivative matrix
    De = DeM[:,:,1]  # derivative of even function
    Do = DoM[:,:,1]  # derivative of odd function

    xnn = x[1:nn]
    y = (L+EXTEND*Lc)*x
    ynn = y[1:nn]

    u = (α ./ (y.+L0).^4 ) - (α ./ (y.-L0).^4)
    unn = u[1:nn]

    A = [ Do   0*De  ; 0*Do    De]./(L-Lc)
    AA = A*diagm([unn;unn])
    #DeD = De*diagm()
    B = [ 0*De  -Do  ;  -(κ/s)*De   -(1/s)*diagm(ones(nn))]

    Af = [D1  0*D1 ; 0*D1 D1]./(L-Lc)
    AAf = Af*diagm([u;u])
    Bf =  [0*D1   -D1 ; -(κ/s)*D1   -(1/s)*diagm(ones(N)) ]

    # done.
    return ring(N,nn,parmd,Ls,y,u,A,AA,B,D1,De,Do,Af,AAf,Bf)
end

function setκ!(r::ring,Σ::Array{Float64,1},β::Float64)
    s = r.parmd["s"]; κ = r.parmd["κ"]
    Deκ = r.De * diagm(Σ.^β)
    B = [ 0*r.De  -r.Do  ;  -(κ/s)*Deκ   -(1/s)*diagm(ones(r.nn))]
    r.B = B
end

end#module
