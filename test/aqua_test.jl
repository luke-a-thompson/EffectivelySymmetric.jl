# Aqua.jl quality assurance. `deps_compat` and `stale_deps` are disabled
# until the package gains full `[compat]` coverage in a follow-up.

using Test
using Aqua
using EffectivelySymmetric

Aqua.test_all(
    EffectivelySymmetric;
    ambiguities = false,
    deps_compat = false,
    stale_deps = false,
    piracies = false,
)
