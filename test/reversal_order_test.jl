# Time-reversal defect order via `solve`.
#
# For an "effectively symmetric" method of antisymmetric order p,
#   e_rev(h) = |Phi_{-h}(Phi_h(y0)) - y0| = O(h^{p+1})
# Expected orders:
#   EES25  -- antisymmetric order 5  ->  defect order ~6
#   EES27  -- antisymmetric order 7  ->  defect order ~8
#
# DiffEqDevTools.test_convergence cannot drive the composition Phi_{-h}∘Phi_h
# (it measures error against an analytic solution), so the defects are computed
# through `solve` and the order is estimated with DiffEqDevTools' shared
# estimator, calc𝒪estimates — the same one ConvergenceSimulation uses.

using Test
using EffectivelySymmetric
using OrdinaryDiffEqCore: ODEProblem, solve
using DiffEqDevTools: calc𝒪estimates

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

# `hs` must halve between entries (calc𝒪estimates averages log2 error ratios)
# and stay coarse enough that the defect sits above round-off saturation.
function _reversal_order(alg, hs; u0 = 0.7)
    errs = [_reversal_error(alg, u0, h) for h in hs]
    return only(last(calc𝒪estimates(:final => errs)))
end

@testset "EES25 reversal order ≈ 6" begin
    hs = [0.4 * 2.0^-k for k in 0:5]  # defects span ~9e-5 .. ~1e-13
    @test _reversal_order(EES25(), hs) ≈ 6 atol = 0.3
end

@testset "EES27 reversal order ≈ 8" begin
    hs = [0.2 * 2.0^-k for k in 0:2]  # defects span ~1.5e-9 .. ~3e-14
    @test _reversal_order(EES27(), hs) ≈ 8 atol = 0.4
end
