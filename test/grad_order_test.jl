# Checkpointed vs reconstructed discrete-adjoint test.
#
# Scalar parameter theta in SO(3) dynamics:
#     dR/dt = hat(omega0 + theta * R * b) * R
# Terminal loss:
#     J(R(T)) = 0.5 * ||R(T) * e3 - target||^2
#
# Compares g_disc,chk(h) (stored forward states) with g_disc,recon(h)
# (one backward step to recover each previous state).
# The difference e_recon(h) = |g_recon - g_chk| isolates reconstruction error.

using Test
using ForwardDiff
using LinearAlgebra

const GRAD_T_FINAL = 5.0
const GRAD_THETA0 = 1.0
const GRAD_H0 = 0.5
const GRAD_M = 7

const GRAD_OMEGA0 = [0.05, 0.30, 1.00]
const GRAD_B = [1.80, 0.60, 0.00]
const GRAD_E3 = [0.0, 0.0, 1.0]
const GRAD_TARGET = normalize([0.25, -0.55, 0.80])
const GRAD_R0 = Matrix{Float64}(I, 3, 3)

function _primal(x)
    if x isa ForwardDiff.Dual
        return ForwardDiff.value(x)
    end
    return x
end

function _hat(v::AbstractVector{T})::Matrix{T} where {T <: Real}
    return T[
         0.0   -v[3]   v[2]
         v[3]   0.0   -v[1]
        -v[2]   v[1]   0.0
    ]
end

function _so3_exp(xi::AbstractVector{T})::Matrix{T} where {T <: Real}
    K = _hat(xi)
    theta2 = sum(abs2, xi)
    theta = sqrt(theta2)

    if abs(_primal(theta)) < 1e-10
        a = one(T) - theta2 / 6 + theta2^2 / 120
        b = one(T) / 2 - theta2 / 24 + theta2^2 / 720
    else
        a = sin(theta) / theta
        b = (one(T) - cos(theta)) / theta2
    end

    return Matrix{T}(I, 3, 3) + a * K + b * (K * K)
end

function _angular_velocity(R::AbstractMatrix{T}, theta::T)::Vector{T} where {T <: Real}
    return T.(GRAD_OMEGA0) .+ theta .* (R * T.(GRAD_B))
end

function _cg2_step(R::AbstractMatrix{T}, h::Float64, theta::T)::Matrix{T} where {T <: Real}
    k1 = h .* _angular_velocity(R, theta)
    y1 = _so3_exp(k1) * R
    k2 = h .* _angular_velocity(y1, theta)
    return _so3_exp(0.5 .* k1) * (_so3_exp(0.5 .* k2) * R)
end

function _ees25_step(R::AbstractMatrix{T}, h::Float64, theta::T)::Matrix{T} where {T <: Real}
    k1 = h .* _angular_velocity(R, theta)
    y1 = _so3_exp(0.5 .* k1) * R

    k2 = h .* _angular_velocity(y1, theta)
    dy2 = -0.5 .* k1 .+ k2
    y2 = _so3_exp(dy2) * y1

    k3 = h .* _angular_velocity(y2, theta)
    dy3 = -2.0 .* dy2 .+ k3
    return _so3_exp(0.25 .* dy3) * y2
end

function _forward_solve(step::Function, h::Float64, theta::T)::Matrix{T} where {T <: Real}
    nsteps = round(Int, GRAD_T_FINAL / h)
    @assert isapprox(nsteps * h, GRAD_T_FINAL; atol = 1e-12, rtol = 0.0)

    R = T.(GRAD_R0)
    for _ in 1:nsteps
        R = step(R, h, theta)
    end
    return R
end

function _terminal_cost_gradient(R::AbstractMatrix{Float64})::Vector{Float64}
    diff = R * GRAD_E3 - GRAD_TARGET
    return vec(diff * transpose(GRAD_E3))
end

function _local_pullback(
    step::Function,
    R_prev::AbstractMatrix{Float64},
    lambda_next::Vector{Float64},
    h::Float64,
    theta::Float64,
)::Tuple{Vector{Float64}, Float64}
    function scalar_map(z::AbstractVector{T}) where {T <: Real}
        R = reshape(z[1:9], 3, 3)
        theta_local = z[10]
        return dot(lambda_next, vec(step(R, h, theta_local)))
    end

    z = vcat(vec(R_prev), theta)
    grad = ForwardDiff.gradient(scalar_map, z)
    return Vector{Float64}(grad[1:9]), grad[10]
end

function _forward_trajectory(step::Function, h::Float64, theta::Float64)::Vector{Matrix{Float64}}
    nsteps = round(Int, GRAD_T_FINAL / h)
    trajectory = Vector{Matrix{Float64}}(undef, nsteps + 1)
    trajectory[1] = copy(GRAD_R0)

    for n in 1:nsteps
        trajectory[n + 1] = Matrix{Float64}(step(trajectory[n], h, theta))
    end

    return trajectory
end

function _adjoint_gradient_checkpointed(step::Function, h::Float64, theta::Float64)::Float64
    trajectory = _forward_trajectory(step, h, theta)
    lambda = _terminal_cost_gradient(trajectory[end])
    grad_theta = 0.0

    for n in length(trajectory)-1:-1:1
        lambda, local_grad_theta = _local_pullback(step, trajectory[n], lambda, h, theta)
        grad_theta += local_grad_theta
    end

    return grad_theta
end

function _adjoint_gradient_reconstructed(step::Function, h::Float64, theta::Float64)::Float64
    nsteps = round(Int, GRAD_T_FINAL / h)
    R = Matrix{Float64}(_forward_solve(step, h, theta))
    lambda = _terminal_cost_gradient(R)
    grad_theta = 0.0

    for _ in nsteps:-1:1
        R_prev = Matrix{Float64}(step(R, -h, theta))
        lambda, local_grad_theta = _local_pullback(step, R_prev, lambda, h, theta)
        grad_theta += local_grad_theta
        R = R_prev
    end

    return grad_theta
end

function _fit_slope(h_values::Vector{Float64}, errors::Vector{Float64})::Float64
    valid = [k for k in eachindex(errors) if errors[k] > 1e-14]
    @assert length(valid) >= 2

    start_idx = max(first(valid), last(valid) - 3)
    sel = start_idx:last(valid)

    x = log.(h_values[sel])
    y = log.(errors[sel])
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
    slope_cg2   = run_grad_study(_cg2_step)
    slope_ees25 = run_grad_study(_ees25_step)

    @test slope_cg2 > 2.5
    @test slope_cg2 < 4.5
    @test slope_ees25 > 5.0
    @test slope_ees25 < 7.0
end
