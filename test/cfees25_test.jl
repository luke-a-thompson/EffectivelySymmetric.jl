# CFEES25 Lie-group integrator driven through `solve`.
#
# (1) On a state-dependent SO(3) ODE, dR/dt = Ω(R) · R, the trajectory stays
#     in SO(3) to round-off.
# (2) Reversal defect through `solve` on the same problem matches the
#     antisymmetric order of the underlying EES25 tableau (~6).

using Test
using Statistics
using LinearAlgebra
using EffectivelySymmetric
using OrdinaryDiffEqCore: ODEProblem, solve
using OrdinaryDiffEqLinear  # activates ESLinearExt
using SciMLOperators: MatrixOperator
using ExponentialUtilities  # extension dep; load just to activate ESLinearExt

const CF_OMEGA0 = [0.05, 0.30, 1.00]
const CF_B      = [1.80, 0.60, 0.00]

function _hat(v)
    return [  0.0   -v[3]  v[2]
              v[3]   0.0  -v[1]
             -v[2]   v[1]  0.0 ]
end

function _update_A!(A, u, p, t)
    R = reshape(u, 3, 3)
    A .= kron(I(3), _hat(CF_OMEGA0 + R * CF_B))
end

function _so3_problem(tspan)
    R0 = Matrix{Float64}(I, 3, 3)
    u0 = vec(R0)
    A0 = zeros(9, 9)
    op = MatrixOperator(A0; update_func! = _update_A!)
    return ODEProblem(op, u0, tspan)
end

@testset "SO(3) orthogonality preserved" begin
    prob = _so3_problem((0.0, 2.0))
    sol = solve(prob, CFEES25(); dt = 0.02, saveat = 0.04)
    max_err = maximum(sol.u) do u
        R = reshape(u, 3, 3)
        norm(R' * R - I(3))
    end
    @test max_err < 1e-10
end

function _cf_reversal_error(h)
    f_prob = _so3_problem((0.0, h))
    u1 = solve(f_prob, CFEES25(); dt = h, adaptive = false,
               save_everystep = false, save_start = false).u[end]

    # SciMLOperators caches the update_func state on the operator; build a
    # fresh problem so the backward solve doesn't reuse stale coefficients.
    A0 = zeros(9, 9)
    op_b = MatrixOperator(A0; update_func! = _update_A!)
    b_prob = ODEProblem(op_b, u1, (h, 0.0))
    u_back = solve(b_prob, CFEES25(); dt = h, adaptive = false,
                   save_everystep = false, save_start = false).u[end]

    return norm(u_back - vec(Matrix{Float64}(I, 3, 3)))
end

@testset "CFEES25 reversal slope ≈ 6" begin
    hs = [0.2 * 2.0^(-k) for k in 0:6]
    errs = [_cf_reversal_error(h) for h in hs]
    keep = [k for k in eachindex(errs) if isfinite(errs[k]) && errs[k] > 1e-13]
    @assert length(keep) >= 4
    x = log.(hs[keep]); y = log.(errs[keep])
    tail = max(1, length(x) - 4):length(x)
    xs = x[tail]; ys = y[tail]
    xm = mean(xs); ym = mean(ys)
    slope = sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2)
    @test 5.0 < slope < 7.5
end
