module cheb
using LinearAlgebra

export chebdif, projmat

"""
flippy():
    Utility function that flips matrix up-down and left-right
        (or, equivalently, rotates it by 180deg).
"""
function flippy(M)
    M[end:-1:1,end:-1:1]
end

"""
chebdif():

Julia translation of the MATLAB routine by J.A.C. Weideman & S.C. Reddy (1998).

Original MATAB documentation:

%  The function DM =  chebdif(N,M) computes the differentiation
%  matrices D1, D2, ..., DM on Chebyshev nodes.
%
%  Input:
%  N:        Size of differentiation matrix.
%  M:        Number of derivatives required (integer).
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced
%  accuracy suggested by W. Don and S. Solomonoff in
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric
%  identities to avoid the computation of differences
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%
% J.A.C. Weideman, S.C. Reddy 1998.

"""
function chebdif(N::Integer,M::Integer)
    I  = diagm(ones(N))             # Identity matrix
    L  = I.!=0                      # Logical identity matrix

    n1 = N÷2;    n2 = n1 + N%2      # Indices used for flipping trick.

    k  = collect(0:N-1)             # Compute theta vector.
    th = k * π/(N-1)

    temp = collect(N-1:-2:1-N)
    x  = sin.(π*temp/(2.0*(N-1)))   # Compute Chebyshev points

    T  = repeat(th/2,1,N)
    DX = 2*sin.(T'+T) .* sin.(T'-T) # Trigonometric identity.

    DX = [DX[1:n1,:]; -flippy(DX[1:n2,:])]
    DX[L] .= 1.0

    temp = (-1.0).^k
    C = temp * temp'  # generates a Toeplitz matrix with diagonal of +1 and alternating signs on diags
    C[1,:]   *= 2;      C[end,:] *= 2
    C[:,1]   *= 0.5;    C[:,end] *= 0.5

    Z = 1.0./DX
    Z[L] .= 0.0
    D = diagm(ones(N))

    s = size(D)
    DM = zeros(s...,M)
    De = zeros(n2,n2,M)
    Do = zeros(n2,n2,M)
    (Pe,Po,Qe,Qo) = projmat(N)

    for ell in 1:M
        D = ell * Z .* (C .* repeat(diag(D),1,N) - D)
        D[L] = -sum(D,dims=2)
        DM[:,:,ell] = D
        if ell%2 == 0 # 2nd, 4th, 6th etc derivatives do not change parity
            De[:,:,ell] = Pe * D * Qe
            Do[:,:,ell] = Po * D * Qo
        else # 1st, 3rd, 5th etc derivatives change parity
            De[:,:,ell] = Po * D * Qe
            Do[:,:,ell] = Pe * D * Qo
        end
    end
    return (x,DM,De,Do)
end

"""
projmat():

    Constructs projection matrices to project (and de-project) onto space of
    even and odd functions, from general space of functions.

    I could try to be clever here, dealing with N odd or even, but I forgo
    cleverness in favor of readability. (In all likelihood this could be done
    in far fewer lines, I'm sure.)

    2*Peve' de-projects an even function back onto the full function space.
    2*Podd' de-projects an odd function back onto the full function space.

    Peve is the inverse of dPeve
    Podd is the inverse of dPodd
"""
function projmat(N::Integer)
    n1 = N÷2;  n2 = n1 + N%2
    isNodd = (n2>n1) # is N odd?
    top = diagm(ones(n2))
    if isNodd
        Qodd = [ top[1:n2,:] ; -top[n1:-1:1,:] ]
        Qodd[n2,:] .*= 0
        Podd = Qodd' *0.5
        Qeve = [ top[1:n2,:] ; top[n1:-1:1,:] ]
        Peve = Qeve' *0.5
        Peve[:,n2] .*= 2
    else
        Qodd = [ top[1:n1,:] ; -top[n1:-1:1,:] ]
        Qeve = [ top[1:n1,:] ;  top[n1:-1:1,:] ]
        Podd = Qodd' *0.5
        Peve = Qeve' *0.5
    end
    return (Peve,Podd,Qeve,Qodd)
end

end#module
