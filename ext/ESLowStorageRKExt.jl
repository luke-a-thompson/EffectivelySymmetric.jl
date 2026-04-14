module ESLowStorageRKExt

using EffectivelySymmetric:
    EES25_2N, EES25_A2end, EES25_B1, EES25_B2end, EES25_C2end

import OrdinaryDiffEqCore: alg_cache, constvalue

using OrdinaryDiffEqLowStorageRK:
    LowStorageRK2NCache, LowStorageRK2NConstantCache

function EES25_2NConstantCache(T, T2)
    A2end = (convert(T, EES25_A2end[1]), convert(T, EES25_A2end[2]))
    B1    = convert(T, EES25_B1)
    B2end = (convert(T, EES25_B2end[1]), convert(T, EES25_B2end[2]))
    c2end = (convert(T2, EES25_C2end[1]), convert(T2, EES25_C2end[2]))
    return LowStorageRK2NConstantCache{2, T, T2}(A2end, B1, B2end, c2end)
end

function alg_cache(
        alg::EES25_2N, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = EES25_2NConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    tmp = zero(u)
    williamson_condition = alg.williamson_condition
    if calck
        k = zero(rate_prototype)
        williamson_condition = false
    else
        if williamson_condition
            k = tmp
        else
            k = zero(rate_prototype)
        end
    end
    return LowStorageRK2NCache(
        u, uprev, k, tmp, tab, williamson_condition, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::EES25_2N, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return EES25_2NConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

end
