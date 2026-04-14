# Time-reversal defect test for EES25 and CG2.
#
# Antisymmetric order p  =>  e_rev(h) = ||Phi_{-h}(Phi_h(y0)) - y0|| = O(h^{p+1})
# Expected slopes:
#   CG2   -- antisymmetric order 2  ->  slope 3
#   EES25 -- antisymmetric order 5  ->  slope 6
#
# Test ODE: dR/dt = L(R)*R on SO(3), L(R) = hat(omega0 + R*b).
# Standalone step functions (no DiffEq), so negative h works directly.

using Test
using LinearAlgebra
using Statistics

hat(v::AbstractVector) = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]

const ORDER_Omega0 = [0.1, 0.2, 0.3]
const ORDER_B  = [0.5, -0.3, 0.1]
L(R) = hat(ORDER_Omega0 + R * ORDER_B)

function cg2_step(R, h::Float64)
    K1 = h * L(R)
    K2 = h * L(exp(K1) * R)
    return exp((1/2) * K1) * (exp((1/2) * K2) * R)
end

function ees25_step(R, h::Float64)
    K1  = h * L(R)
    Y1  = exp((1/2) * K1) * R
    K2  = h * L(Y1)
    dY2 = -(1/2) * K1 + K2
    Y2  = exp(dY2) * Y1
    K3  = h * L(Y2)
    dY3 = -2 * dY2 + K3
    return exp((1/4) * dY3) * Y2
end

function reversal_order(step::Function, R0, h0, m)
    errors = [norm(step(step(R0, h0 * 2.0^(-k)), -(h0 * 2.0^(-k))) - R0)
              for k in 0:m]
    slopes = [log(errors[k] / errors[k+1]) / log(2) for k in 1:m]

    stable = Float64[]
    for k in 1:m
        errors[k+1] < 1e-12 && break
        errors[k+1] > errors[k] && break
        k < 2 && continue
        push!(stable, slopes[k])
    end
    return isempty(stable) ? mean(slopes[2:max(2, m÷2)]) : mean(stable)
end

@testset "Time-reversal order" begin
    R0 = exp(hat([0.8, 0.5, 0.3]))
    h0 = 0.5
    m  = 20

    p_cg2   = reversal_order(cg2_step, R0, h0, m)
    p_ees25 = reversal_order(ees25_step, R0, h0, m)

    @test p_cg2 > 3.5
    @test p_cg2 < 4.5
    @test p_ees25 > 5.5
    @test p_ees25 < 6.5
end
