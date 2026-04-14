"""
    EES25()

3-stage effectively-symmetric integrator (standard explicit RK form).
Uses `OrdinaryDiffEqCore` only -- always available when the package is loaded.
"""
struct EES25 <: OrdinaryDiffEqAlgorithm end

"""
    EES25_2N(; williamson_condition = true)

3-stage effectively-symmetric integrator in 2N low-storage form.
Requires `OrdinaryDiffEqLowStorageRK` to be loaded.
"""
Base.@kwdef struct EES25_2N{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end

"""
    CFEES25()

3-stage commutator-free Lie-group integrator on SO(n).
Requires `OrdinaryDiffEqLinear`, `ExponentialUtilities`, and `SciMLOperators` to be loaded.

The update rule replaces additions with matrix exponentials:

    ΔY₁ = hK₁,           Y₁ = exp(½ΔY₁)·Y₀
    ΔY₂ = -½ΔY₁ + hK₂,  Y₂ = exp(ΔY₂)·Y₁
    ΔY₃ = -2ΔY₂  + hK₃,  Y_{t+h} = exp(¼ΔY₃)·Y₂
"""
struct CFEES25 <: OrdinaryDiffEqAlgorithm end
