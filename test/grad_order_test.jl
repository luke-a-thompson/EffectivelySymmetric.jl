# Checkpointed vs reconstructed discrete-adjoint test for the supplied EES25
# and EES27 tableaux on a scalar nonlinear ODE.

using Test
using ForwardDiff

const GRAD_T_FINAL = 5.0
const GRAD_THETA0 = 0.9
const GRAD_H0 = 5.0 / 8.0
const GRAD_M = 6
const GRAD_U0 = 0.7
const GRAD_TARGET = -0.15

f_grad(u, theta) = theta * sin(u) + 0.1 * u^2

function rk_tableau_step(u, h, theta, A, b)
    s = length(b)
    k = Vector{typeof(u)}(undef, s)

    for i in 1:s
        stage = u
        for j in 1:i-1
            stage += h * A[i, j] * k[j]
        end
        k[i] = f_grad(stage, theta)
    end

    out = u
    for i in 1:s
        out += h * b[i] * k[i]
    end
    return out
end

const GRAD_A25 = [
    0.0 0.0 0.0
    1 / 3 0.0 0.0
    -5 / 48 15 / 16 0.0
]
const GRAD_B25 = (1 / 10, 1 / 2, 2 / 5)

const GRAD_R = sqrt(2.0)
const GRAD_A27 = [
    0.0 0.0 0.0 0.0
    (2 - GRAD_R) / 3 0.0 0.0 0.0
    (-4 + GRAD_R) / 24 (4 + GRAD_R) / 8 0.0 0.0
    (-176 + 145 * GRAD_R) / 168 3 * (8 - 5 * GRAD_R) / 56 3 * (3 - GRAD_R) / 7 0.0
]
const GRAD_B27 = (
    (5 - 3 * GRAD_R) / 14,
    (3 + GRAD_R) / 14,
    3 * (-1 + 2 * GRAD_R) / 14,
    (9 - 4 * GRAD_R) / 14,
)

function _forward_solve(step::Function, h::Float64, theta)
    nsteps = round(Int, GRAD_T_FINAL / h)
    @assert isapprox(nsteps * h, GRAD_T_FINAL; atol = 1e-12, rtol = 0.0)

    u = theta * 0 + GRAD_U0
    for _ in 1:nsteps
        u = step(u, h, theta)
    end
    return u
end

_terminal_cost_gradient(u) = u - GRAD_TARGET

function _local_pullback(
    step::Function,
    u_prev::Float64,
    lambda_next::Float64,
    h::Float64,
    theta::Float64,
)
    scalar_map(z) = lambda_next * step(z[1], h, z[2])
    grad = ForwardDiff.gradient(scalar_map, [u_prev, theta])
    return grad[1], grad[2]
end

function _forward_trajectory(step::Function, h::Float64, theta::Float64)
    nsteps = round(Int, GRAD_T_FINAL / h)
    trajectory = Vector{Float64}(undef, nsteps + 1)
    trajectory[1] = GRAD_U0

    for n in 1:nsteps
        trajectory[n + 1] = step(trajectory[n], h, theta)
    end

    return trajectory
end

function _adjoint_gradient_checkpointed(step::Function, h::Float64, theta::Float64)
    trajectory = _forward_trajectory(step, h, theta)
    lambda = _terminal_cost_gradient(trajectory[end])
    grad_theta = 0.0

    for n in length(trajectory)-1:-1:1
        lambda, local_grad_theta = _local_pullback(step, trajectory[n], lambda, h, theta)
        grad_theta += local_grad_theta
    end

    return grad_theta
end

function _adjoint_gradient_reconstructed(step::Function, h::Float64, theta::Float64)
    nsteps = round(Int, GRAD_T_FINAL / h)
    u = _forward_solve(step, h, theta)
    lambda = _terminal_cost_gradient(u)
    grad_theta = 0.0

    for _ in nsteps:-1:1
        u_prev = step(u, -h, theta)
        lambda, local_grad_theta = _local_pullback(step, u_prev, lambda, h, theta)
        grad_theta += local_grad_theta
        u = u_prev
    end

    return grad_theta
end

function _fit_slope(h_values::Vector{Float64}, errors::Vector{Float64})
    valid = [k for k in eachindex(errors) if errors[k] > 1e-14]
    @assert length(valid) >= 2

    x = log.(h_values[valid])
    y = log.(errors[valid])
    x_mean = sum(x) / length(x)
    y_mean = sum(y) / length(y)
    return sum((x .- x_mean) .* (y .- y_mean)) / sum((x .- x_mean) .^ 2)
end

function run_grad_study(step::Function)
    h_values = [GRAD_H0 * 2.0^(-k) for k in 0:GRAD_M]
    recon_errors = Float64[]

    for h in h_values
        g_chk = _adjoint_gradient_checkpointed(step, h, GRAD_THETA0)
        g_recon = _adjoint_gradient_reconstructed(step, h, GRAD_THETA0)
        push!(recon_errors, abs(g_recon - g_chk))
    end

    return _fit_slope(h_values, recon_errors)
end

@testset "Discrete adjoint reconstruction order" begin
    slope_ees25 = run_grad_study((u, h, theta) -> rk_tableau_step(u, h, theta, GRAD_A25, GRAD_B25))
    slope_ees27 = run_grad_study((u, h, theta) -> rk_tableau_step(u, h, theta, GRAD_A27, GRAD_B27))

    @test slope_ees25 > 4.5
    @test slope_ees25 < 7.0
    @test slope_ees27 > 6.0
    @test slope_ees27 < 9.0
end
