# Time-reversal defect test for the supplied EES25 and EES27 tableaux.
#
# Antisymmetric order p  =>  e_rev(h) = |Phi_{-h}(Phi_h(y0)) - y0| = O(h^{p+1})
# Expected slopes:
#   EES25 -- antisymmetric order 5  ->  slope 6
#   EES27 -- antisymmetric order 7  ->  slope 8

using Test
using Statistics

const ORDER_THETA = 0.3

f_order(u) = sin(u) + ORDER_THETA * u^2

function rk_tableau_step(u, h, A, b)
    s = length(b)
    k = Vector{Float64}(undef, s)

    for i in 1:s
        stage = u
        for j in 1:i-1
            stage += h * A[i, j] * k[j]
        end
        k[i] = f_order(stage)
    end

    out = u
    for i in 1:s
        out += h * b[i] * k[i]
    end
    return out
end

function reversal_order(step::Function, u0, h0, m)
    errors = [abs(step(step(u0, h0 * 2.0^(-k)), -(h0 * 2.0^(-k))) - u0)
              for k in 0:m]
    slopes = [log(errors[k] / errors[k + 1]) / log(2) for k in 1:m if errors[k + 1] > 0]

    stable = Float64[]
    for k in 3:length(slopes)
        !isfinite(slopes[k]) && break
        push!(stable, slopes[k])
    end
    return mean(stable)
end

const ORDER_A25 = [
    0.0 0.0 0.0
    1 / 3 0.0 0.0
    -5 / 48 15 / 16 0.0
]
const ORDER_B25 = (1 / 10, 1 / 2, 2 / 5)

const ORDER_R = sqrt(2.0)
const ORDER_A27 = [
    0.0 0.0 0.0 0.0
    (2 - ORDER_R) / 3 0.0 0.0 0.0
    (-4 + ORDER_R) / 24 (4 + ORDER_R) / 8 0.0 0.0
    (-176 + 145 * ORDER_R) / 168 3 * (8 - 5 * ORDER_R) / 56 3 * (3 - ORDER_R) / 7 0.0
]
const ORDER_B27 = (
    (5 - 3 * ORDER_R) / 14,
    (3 + ORDER_R) / 14,
    3 * (-1 + 2 * ORDER_R) / 14,
    (9 - 4 * ORDER_R) / 14,
)

@testset "Time-reversal order" begin
    u0 = 0.7
    h0 = 0.4
    m  = 12

    p_ees25 = reversal_order((u, h) -> rk_tableau_step(u, h, ORDER_A25, ORDER_B25), u0, h0, m)
    p_ees27 = reversal_order((u, h) -> rk_tableau_step(u, h, ORDER_A27, ORDER_B27), u0, h0, m)

    @test p_ees25 > 5.5
    @test p_ees25 < 6.5
    @test p_ees27 > 7.0
    @test p_ees27 < 8.5
end
