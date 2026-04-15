# Standard (non-low-storage, non-Lie-group) effectively-symmetric implementation.
# Uses the supplied explicit RK tableaux directly.

import OrdinaryDiffEqCore: OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    @cache, get_fsalfirstlast

const ExplicitEESAlgorithm = Union{EES25, EES27}

# ── Caches ──────────────────────────────────────────────────────────────────

@cache struct EESCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    fsalfirst::rateType
end

get_fsalfirstlast(cache::EESCache, u) = (cache.fsalfirst, cache.k1)

struct EESConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(
        alg::ExplicitEESAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1        = zero(rate_prototype)
    k2        = zero(rate_prototype)
    k3        = zero(rate_prototype)
    k4        = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    return EESCache(u, uprev, uprev2, zero(u), k1, k2, k3, k4, fsalfirst)
end

function alg_cache(
        alg::ExplicitEESAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return EESConstantCache()
end

# ── initialize! ─────────────────────────────────────────────────────────────

function initialize!(integrator, cache::EESCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.fsalfirst
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::EESConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

# ── perform_step! (in-place) ────────────────────────────────────────────────

function perform_step!(integrator, cache::EESCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, k1, k2, k3, k4, fsalfirst) = cache

    if integrator.alg isa EES25
        copyto!(k1, fsalfirst)
        @. tmp = uprev + dt * EES25_A21 * k1

        f(k2, tmp, p, t + EES25_C[2] * dt)
        increment_nf!(integrator.stats, 1)
        @. tmp = uprev + dt * (EES25_A31 * k1 + EES25_A32 * k2)

        f(k3, tmp, p, t + EES25_C[3] * dt)
        increment_nf!(integrator.stats, 1)
        @. u = uprev + dt * (EES25_B[1] * k1 + EES25_B[2] * k2 + EES25_B[3] * k3)
    else
        copyto!(k1, fsalfirst)
        @. tmp = uprev + dt * EES27_A21 * k1

        f(k2, tmp, p, t + EES27_C[2] * dt)
        increment_nf!(integrator.stats, 1)
        @. tmp = uprev + dt * (EES27_A31 * k1 + EES27_A32 * k2)

        f(k3, tmp, p, t + EES27_C[3] * dt)
        increment_nf!(integrator.stats, 1)
        @. tmp = uprev + dt * (EES27_A41 * k1 + EES27_A42 * k2 + EES27_A43 * k3)

        f(k4, tmp, p, t + EES27_C[4] * dt)
        increment_nf!(integrator.stats, 1)
        @. u = uprev + dt * (EES27_B[1] * k1 + EES27_B[2] * k2 + EES27_B[3] * k3 + EES27_B[4] * k4)
    end

    f(fsalfirst, u, p, t + dt)
    increment_nf!(integrator.stats, 1)
end

# ── perform_step! (out-of-place) ────────────────────────────────────────────

function perform_step!(integrator, cache::EESConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    k1 = integrator.fsalfirst

    if integrator.alg isa EES25
        tmp = uprev + dt * EES25_A21 * k1
        k2  = f(tmp, p, t + EES25_C[2] * dt)
        increment_nf!(integrator.stats, 1)
        tmp = uprev + dt * (EES25_A31 * k1 + EES25_A32 * k2)
        k3  = f(tmp, p, t + EES25_C[3] * dt)
        increment_nf!(integrator.stats, 1)
        u   = uprev + dt * (EES25_B[1] * k1 + EES25_B[2] * k2 + EES25_B[3] * k3)
    else
        tmp = uprev + dt * EES27_A21 * k1
        k2  = f(tmp, p, t + EES27_C[2] * dt)
        increment_nf!(integrator.stats, 1)
        tmp = uprev + dt * (EES27_A31 * k1 + EES27_A32 * k2)
        k3  = f(tmp, p, t + EES27_C[3] * dt)
        increment_nf!(integrator.stats, 1)
        tmp = uprev + dt * (EES27_A41 * k1 + EES27_A42 * k2 + EES27_A43 * k3)
        k4  = f(tmp, p, t + EES27_C[4] * dt)
        increment_nf!(integrator.stats, 1)
        u   = uprev + dt * (EES27_B[1] * k1 + EES27_B[2] * k2 + EES27_B[3] * k3 + EES27_B[4] * k4)
    end

    integrator.k[1] = k1
    integrator.fsalfirst = f(u, p, t + dt)
    increment_nf!(integrator.stats, 1)
    integrator.u = u
end
