# Time-reversal defect order via `solve`.
#
# For an "effectively symmetric" method of antisymmetric order p,
#   e_rev(h) = |Phi_{-h}(Phi_h(y0)) - y0| = O(h^{p+1})
# Expected slopes:
#   EES25  -- antisymmetric order 5  ->  slope ~6
#   EES27  -- antisymmetric order 7  ->  slope ~8

using Test
using Statistics
using EffectivelySymmetric
using OrdinaryDiffEqCore: ODEProblem, solve

const REV_THETA = 0.3
f_rev(u, p, t) = sin(u) + REV_THETA * u^2

function _one_step(alg, u0, t0, t1, h)
    sol = solve(ODEProblem(f_rev, u0, (t0, t1)), alg;
                dt = h, adaptive = false,
                save_everystep = false, save_start = false)
    return sol.u[end]
end

function _reversal_error(alg, u0, h)
    u1 = _one_step(alg, u0, 0.0, h, h)
    u_back = _one_step(alg, u1, h, 0.0, h)
    return abs(u_back - u0)
end

function _reversal_slope(alg; u0 = 0.7, h0 = 0.4, m = 10)
    hs = [h0 * 2.0^(-k) for k in 0:m]
    errs = [_reversal_error(alg, u0, h) for h in hs]
    keep = [k for k in eachindex(errs) if isfinite(errs[k]) && errs[k] > 1e-14]
    @assert length(keep) >= 3
    x = log.(hs[keep]); y = log.(errs[keep])
    # use only the asymptotic tail (skip the coarsest two points)
    tail = max(1, length(x) - 5):length(x)
    xs = x[tail]; ys = y[tail]
    xm = mean(xs); ym = mean(ys)
    return sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2)
end

@testset "EES25 reversal slope ≈ 6" begin
    s = _reversal_slope(EES25())
    @test 5.3 < s < 6.7
end

@testset "EES27 reversal slope ≈ 8" begin
    s = _reversal_slope(EES27())
    @test 6.8 < s < 8.7
end
