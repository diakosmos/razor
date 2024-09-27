# just for testing/debugging purposes
function junky(y,hrel)
    ℓ₀ = 3
    if true
        @assert -1 ≤ y ≤ 1
        @assert hrel ≥ 0.0
    end
    ξL = y * (ℓ₀-1) + ℓ₀ #  1 ≤ ξL = position relative to mass on left  (dim'less)
    ξR = y * (ℓ₀-1) - ℓ₀ # -1 ≥ ξR = position relative to mass on right (dim'less)
    if hrel ≈ 0.0
        return ξL^-4 - ξR^-4
    else
        aa = 0.711557; bb = -7.58607
        function g(hox)
            min(1/2.5,abs(hox))
        end
        function f(hox)
            (1+ aa*g(hox) + bb*g(hox)^2)^(-1)
        end
        function ww(ξ)
            sign(ξ)/ξ^4 * f(hrel/ξ)
        end
        wL = ww(ξL)
        wR = ww(ξR)
        return wL + wR
    end
end
