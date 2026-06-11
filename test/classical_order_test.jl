# Classical convergence order via DiffEqDevTools.test_convergence.
#
# Drives EES25, EES27, EES25_2N, and EES27_2N on linear and nonlinear scalar
# test problems with known analytic solutions. All four methods have classical
# order 2 (the 5/7 in the names is the antisymmetric order).

using Test
using EffectivelySymmetric
using OrdinaryDiffEqCore: ODEProblem, ODEFunction
using OrdinaryDiffEqLowStorageRK  # activates ESLowStorageRKExt
using DiffEqDevTools: test_convergence

const CL_DTS = 1 .// 2 .^ (8:-1:4)
const CL_TEST_TOL = 0.25

f_linear(u, p, t) = -u
f_linear_analytic(u0, p, t) = u0 * exp(-t)
const PROB_LINEAR = ODEProblem(
    ODEFunction(f_linear; analytic = f_linear_analytic),
    1.0, (0.0, 1.0)
)

f_nonlinear(u, p, t) = u^2
f_nonlinear_analytic(u0, p, t) = u0 / (1 - u0 * t)
const PROB_NONLINEAR = ODEProblem(
    ODEFunction(f_nonlinear; analytic = f_nonlinear_analytic),
    0.5, (0.0, 1.0)
)

@testset "$(nameof(typeof(alg))) classical order 2" for alg in
        (EES25(), EES27(), EES25_2N(), EES27_2N())
    for prob in (PROB_LINEAR, PROB_NONLINEAR)
        sim = test_convergence(CL_DTS, prob, alg)
        @test sim.𝒪est[:final] ≈ 2 atol = CL_TEST_TOL
    end
end
