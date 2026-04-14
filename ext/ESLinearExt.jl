module ESLinearExt

using EffectivelySymmetric: CFEES25

import OrdinaryDiffEqCore:
    OrdinaryDiffEqConstantCache,
    @cache, alg_cache, initialize!, perform_step!,
    _vec, increment_nf!

using OrdinaryDiffEqLinear: LinearMutableCache
using ExponentialUtilities: exponential!, ExpMethodHigham2005
using SciMLOperators: update_coefficients!
using LinearAlgebra: mul!, rmul!

# ── Caches ──────────────────────────────────────────────────────────────────

@cache struct CFEES25Cache{uType, rateType, WType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    Wtmp::WType
    k::rateType
end

struct CFEES25ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(
        alg::CFEES25, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W         = false .* _vec(rate_prototype) .* _vec(rate_prototype)'
    Wtmp      = similar(W)
    k         = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    return CFEES25Cache(u, uprev, uprev2, zero(u), fsalfirst, W, Wtmp, k)
end

function alg_cache(
        alg::CFEES25, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return CFEES25ConstantCache()
end

# ── initialize! ─────────────────────────────────────────────────────────────

function initialize!(integrator, cache::CFEES25Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return increment_nf!(integrator.stats, 1)
end

# ── perform_step! ───────────────────────────────────────────────────────────
#
# Stages (Kᵢ = h·A evaluated at the given point; ΔYᵢ are Lie-algebra elements):
#
#   K₁  = h·F(t,        Y₀)          ΔY₁ = K₁
#   Y₁  = exp(½ΔY₁) · Y₀
#   K₂  = h·F(t+½h,     Y₁)          ΔY₂ = -½ΔY₁ + K₂
#   Y₂  = exp(ΔY₂)  · Y₁
#   K₃  = h·F(t+h,      Y₂)          ΔY₃ = -2ΔY₂  + K₃
#   Y_{t+h} = exp(¼ΔY₃) · Y₂

function perform_step!(integrator, cache::CFEES25Cache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    (; tmp, k, W, Wtmp) = cache
    exp_method = ExpMethodHigham2005()
    L = integrator.f.f

    # Stage 1: K₁ = h·A(t, Y₀)
    update_coefficients!(L, uprev, p, t)
    copyto!(W, convert(AbstractMatrix, L))
    rmul!(W, dt)

    # Stage 2: Y₁ = exp(½K₁)·Y₀
    copyto!(Wtmp, W)
    rmul!(Wtmp, 1 / 2)
    exponential!(Wtmp, exp_method)
    mul!(tmp, Wtmp, uprev)

    # K₂ = h·A(t+½h, Y₁),  ΔY₂ = -½ΔY₁ + K₂
    update_coefficients!(L, tmp, p, t + dt / 2)
    copyto!(Wtmp, convert(AbstractMatrix, L))
    rmul!(Wtmp, dt)
    @. W = Wtmp - (1 / 2) * W

    # Stage 3: Y₂ = exp(ΔY₂)·Y₁
    copyto!(Wtmp, W)
    exponential!(Wtmp, exp_method)
    mul!(k, Wtmp, tmp)

    # K₃ = h·A(t+h, Y₂),  ΔY₃ = -2ΔY₂ + K₃
    update_coefficients!(L, k, p, t + dt)
    copyto!(Wtmp, convert(AbstractMatrix, L))
    rmul!(Wtmp, dt)
    @. W = Wtmp - 2 * W

    # Update: Y_{t+h} = exp(¼ΔY₃)·Y₂
    copyto!(Wtmp, W)
    rmul!(Wtmp, 1 / 4)
    exponential!(Wtmp, exp_method)
    mul!(u, Wtmp, k)

    integrator.f(integrator.fsallast, u, p, t + dt)
    return increment_nf!(integrator.stats, 1)
end

end
