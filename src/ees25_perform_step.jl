# Standard (non-low-storage, non-Lie-group) EES25 implementation.
# Uses the 2N recurrence as additions: tmp = A*tmp + dt*k, u = u + B*tmp.

import OrdinaryDiffEqCore: OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    @cache, get_fsalfirstlast

# ── Caches ──────────────────────────────────────────────────────────────────

@cache struct EES25Cache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    k::rateType
    fsalfirst::rateType
end

get_fsalfirstlast(cache::EES25Cache, u) = (cache.fsalfirst, cache.k)

struct EES25ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(
        alg::EES25, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k         = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    return EES25Cache(u, uprev, uprev2, zero(u), k, fsalfirst)
end

function alg_cache(
        alg::EES25, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return EES25ConstantCache()
end

# ── initialize! ─────────────────────────────────────────────────────────────

function initialize!(integrator, cache::EES25Cache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.fsalfirst
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::EES25ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

# ── perform_step! (in-place) ────────────────────────────────────────────────

function perform_step!(integrator, cache::EES25Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, k, fsalfirst) = cache

    # Stage 1: tmp = dt * f(uprev), u = uprev + B1 * tmp
    f(fsalfirst, uprev, p, t)
    increment_nf!(integrator.stats, 1)
    @. tmp = dt * fsalfirst
    @. u = uprev + EES25_B1 * tmp

    # Stage 2: tmp = A[1]*tmp + dt*f(u), u = u + B2end[1]*tmp
    f(k, u, p, t + EES25_C2end[1] * dt)
    increment_nf!(integrator.stats, 1)
    @. tmp = EES25_A2end[1] * tmp + dt * k
    @. u = u + EES25_B2end[1] * tmp

    # Stage 3: tmp = A[2]*tmp + dt*f(u), u = u + B2end[2]*tmp
    f(k, u, p, t + EES25_C2end[2] * dt)
    increment_nf!(integrator.stats, 1)
    @. tmp = EES25_A2end[2] * tmp + dt * k
    @. u = u + EES25_B2end[2] * tmp
end

# ── perform_step! (out-of-place) ────────────────────────────────────────────

function perform_step!(integrator, cache::EES25ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator

    # Stage 1
    tmp = dt * integrator.fsalfirst
    u   = uprev + EES25_B1 * tmp

    # Stage 2
    k   = f(u, p, t + EES25_C2end[1] * dt)
    increment_nf!(integrator.stats, 1)
    tmp = EES25_A2end[1] * tmp + dt * k
    u   = u + EES25_B2end[1] * tmp

    # Stage 3
    k   = f(u, p, t + EES25_C2end[2] * dt)
    increment_nf!(integrator.stats, 1)
    tmp = EES25_A2end[2] * tmp + dt * k
    u   = u + EES25_B2end[2] * tmp

    integrator.k[1] = integrator.fsalfirst
    integrator.fsalfirst = f(u, p, t + dt)
    increment_nf!(integrator.stats, 1)
    integrator.u = u
end
