# Classical convergence order via the actual `solve` interface.
#
# Drives EES25, EES27, and EES25_2N through `solve(prob, alg; dt = h, adaptive = false)`
# on a scalar linear test problem and fits the slope of log(err) vs log(h).

using Test
using Statistics
using EffectivelySymmetric
using OrdinaryDiffEqCore: ODEProblem, solve
using OrdinaryDiffEqLowStorageRK  # activates ESLowStorageRKExt

f_linear(u, p, t) = -u
const CL_U0 = 1.0
const CL_TSPAN = (0.0, 1.0)
const CL_EXACT = exp(-1.0)
const CL_DTS = [1 / 16, 1 / 32, 1 / 64, 1 / 128, 1 / 256]

function _convergence_slope(alg)
    errs = Float64[]
    for dt in CL_DTS
        prob = ODEProblem(f_linear, CL_U0, CL_TSPAN)
        sol = solve(prob, alg; dt = dt, adaptive = false,
                    save_everystep = false, save_start = false)
        push!(errs, abs(sol.u[end] - CL_EXACT))
    end
    x = log.(CL_DTS); y = log.(errs)
    xm = mean(x); ym = mean(y)
    return sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
end

@testset "EES25 classical order" begin
    s = _convergence_slope(EES25())
    @test 1.7 < s < 3.5
end

@testset "EES27 classical order" begin
    s = _convergence_slope(EES27())
    @test 1.7 < s < 4.5
end

@testset "EES25_2N classical order" begin
    s = _convergence_slope(EES25_2N())
    @test 1.7 < s < 3.5
end
