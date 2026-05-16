module ESLinearExt

using EffectivelySymmetric:
    CFEES25, CFEES27,
    EES25_A2end, EES25_B1, EES25_B2end, EES25_C2end,
    EES27_A2end, EES27_B1, EES27_B2end, EES27_C2end

import OrdinaryDiffEqCore:
    OrdinaryDiffEqConstantCache,
    @cache, alg_cache, initialize!, perform_step!,
    _vec, increment_nf!, constvalue

using OrdinaryDiffEqLinear: LinearMutableCache
using ExponentialUtilities: exponential!, ExpMethodHigham2005, alloc_mem
using SciMLOperators: update_coefficients!
using LinearAlgebra: mul!, rmul!

const CFEESAlgorithm = Union{CFEES25, CFEES27}

# ── Caches ──────────────────────────────────────────────────────────────────
#
# Commutator-free EES methods share the same 2N-style recurrence:
#
#   ΔY₁ = K₁ = h·L(Y₀, t)
#   Y₁  = exp(B₁·ΔY₁) · Y₀
#   For i = 2..s:
#       Kᵢ  = h·L(Yᵢ₋₁, t + cᵢ·h)
#       ΔYᵢ = Aᵢ·ΔYᵢ₋₁ + Kᵢ
#       Yᵢ  = exp(Bᵢ·ΔYᵢ) · Yᵢ₋₁
#
# CFEES25 and CFEES27 differ only in the (Aᵢ, Bᵢ, cᵢ) coefficients and the
# total stage count, so one parameterised cache + one loop-based perform_step!
# handles both — matching the LowStorageRK2N idiom.

struct CFEESConstantCache{N, T, T2} <: OrdinaryDiffEqConstantCache
    A2end::NTuple{N, T}
    B1::T
    B2end::NTuple{N, T}
    C2end::NTuple{N, T2}
end

function CFEES25ConstantCache(T, T2)
    A2end = (convert(T, EES25_A2end[1]), convert(T, EES25_A2end[2]))
    B1    = convert(T, EES25_B1)
    B2end = (convert(T, EES25_B2end[1]), convert(T, EES25_B2end[2]))
    C2end = (convert(T2, EES25_C2end[1]), convert(T2, EES25_C2end[2]))
    return CFEESConstantCache{2, T, T2}(A2end, B1, B2end, C2end)
end

function CFEES27ConstantCache(T, T2)
    A2end = (
        convert(T, EES27_A2end[1]),
        convert(T, EES27_A2end[2]),
        convert(T, EES27_A2end[3]),
    )
    B1    = convert(T, EES27_B1)
    B2end = (
        convert(T, EES27_B2end[1]),
        convert(T, EES27_B2end[2]),
        convert(T, EES27_B2end[3]),
    )
    C2end = (
        convert(T2, EES27_C2end[1]),
        convert(T2, EES27_C2end[2]),
        convert(T2, EES27_C2end[3]),
    )
    return CFEESConstantCache{3, T, T2}(A2end, B1, B2end, C2end)
end

_cfees_tab(::CFEES25, T, T2) = CFEES25ConstantCache(T, T2)
_cfees_tab(::CFEES27, T, T2) = CFEES27ConstantCache(T, T2)

@cache struct CFEESCache{uType, rateType, WType, expType, tabType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    Wtmp::WType
    k::rateType
    exp_cache::expType
    tab::tabType
end

function alg_cache(
        alg::CFEESAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W         = false .* _vec(rate_prototype) .* _vec(rate_prototype)'
    Wtmp      = similar(W)
    exp_cache = alloc_mem(W, ExpMethodHigham2005())
    k         = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    tab       = _cfees_tab(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return CFEESCache(u, uprev, uprev2, zero(u), fsalfirst, W, Wtmp, k, exp_cache, tab)
end

function alg_cache(
        alg::CFEESAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return _cfees_tab(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

# ── initialize! ─────────────────────────────────────────────────────────────

function initialize!(integrator, cache::CFEESCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return increment_nf!(integrator.stats, 1)
end

# ── perform_step! ───────────────────────────────────────────────────────────

function perform_step!(integrator, cache::CFEESCache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    (; tmp, k, W, Wtmp, exp_cache, tab) = cache
    (; A2end, B1, B2end, C2end) = tab
    exp_method = ExpMethodHigham2005()
    L = integrator.f.f
    s = length(B2end) + 1

    # Stage 1: ΔY₁ = K₁ = h·L(Y₀, t), Y₁ = exp(B₁·ΔY₁)·Y₀
    update_coefficients!(L, uprev, p, t)
    copyto!(W, convert(AbstractMatrix, L))
    rmul!(W, dt)

    copyto!(Wtmp, W)
    rmul!(Wtmp, B1)
    exponential!(Wtmp, exp_method, exp_cache)
    mul!(tmp, Wtmp, uprev)

    # Stages 2..s: ping-pong Y between `tmp` and `k`; write final into `u`.
    Y_prev = tmp
    Y_curr = k
    @inbounds for i in 2:s
        update_coefficients!(L, Y_prev, p, t + C2end[i - 1] * dt)
        copyto!(Wtmp, convert(AbstractMatrix, L))
        rmul!(Wtmp, dt)
        @. W = A2end[i - 1] * W + Wtmp

        copyto!(Wtmp, W)
        rmul!(Wtmp, B2end[i - 1])
        exponential!(Wtmp, exp_method, exp_cache)

        target = (i == s) ? u : Y_curr
        mul!(target, Wtmp, Y_prev)

        if i < s
            Y_prev, Y_curr = Y_curr, Y_prev
        end
    end

    integrator.f(integrator.fsallast, u, p, t + dt)
    return increment_nf!(integrator.stats, 1)
end

end
