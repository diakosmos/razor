module Classical
using LinearAlgebra
#using StaticArrays

function visco(Σrel,β)
    return Diagonal(Σrel.^β)
end

"""
A planetary particulate disk is called a "ring", but this terminology is confusing.
"""
mutable struct ring
#    α::Float64
#    Δx::Float64
#    Δt::Float64
#    xi, xf ::Float64
    N::Integer
    β::Float64
    ΔtVN::Float64
    # Defined on nodes:
    X_::Array{Float64,1}
    U_::Array{Float64,1}
    # Defined on cells:
    x_::Array{Float64,1}
    σ_::Array{Float64,1}
    Gᵐ::Any # Tridiagonal
#    μ_::Array{Float64,1}
#    ν_::Array{Float64,1}

end
function ring(N=5,β=0,xf=3.0,xi=xf/N)
    xtemp_ = collect(range(xi,xf,length=N+1))
    X_ = xtemp_[2:end]
    U_ = 1.0 ./ X_.^4
    ΔX_ = diff(xtemp_)
    x_ = 0.5*(xtemp_[1:end-1]+xtemp_[2:end]) # midpoints
    σ_ = ones(N)
    νᵐ = visco(σ_,β)
    Mᵐ = Diagonal(ΔX_)
    Mⁱᵐ= Diagonal(1.0./ΔX_)
    Kdiag = [1.0 ./ diff(x_)..., 0.0]
    Ksuper = -Kdiag[1:end-1]
    Ksub = zeros(N-1)
    Kᵐ = Tridiagonal(Ksub,Kdiag,Ksuper)
    Dᵐ = Tridiagonal(ones(N-1),-ones(N),zeros(N-1))
    Gᵐ = Dᵐ * Kᵐ * νᵐ
    μ_ = Mᵐ * σ_ # mass (per azimuthal length) in annulus
    ΔtVN = 0.5 * minimum(diff(X_).^2 ./ diag(νᵐ)[1:end-1] )
    return ring(N,β,ΔtVN,X_,U_,x_,σ_,Gᵐ)#,μ_)
end

function diffuse!(r::ring;Δt=r.ΔtVN,α=0.5)
    Gᵐ = r.Gᵐ
    N = r.N
    Aᵐ = I(N) -    α  * Δt * Gᵐ
    Bᵐ = I(N) + (1-α) * Δt * Gᵐ
    # Solve:
    r.σ_ = Aᵐ \ (Bᵐ * r.σ_)
    r.σ_[1]   = 0.0
    r.σ_[end] = 1.0
end

end#module
