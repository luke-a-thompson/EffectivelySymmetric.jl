# EffectivelySymmetric.jl

Effectively-symmetric fixed-step integrators for the [SciML](https://sciml.ai) ecosystem.

The package currently provides:

- `EES25`: revised 3-stage explicit RK tableau.
- `EES27`: new 4-stage explicit RK tableau.
- `EES25_2N` and `CFEES25`: the existing low-storage and commutator-free EES25 variants.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/<org>/EffectivelySymmetric.jl")
```

## Algorithms

### `EES25()` and `EES27()`: Standard ODE

Always available. Works with any `ODEProblem`.

```julia
using EffectivelySymmetric, OrdinaryDiffEqCore

prob = ODEProblem(f, u0, tspan)
sol25 = solve(prob, EES25(); dt = 0.01)
sol27 = solve(prob, EES27(); dt = 0.01)
```

### `EES25_2N()`: 2N Low-Storage

Activated by loading `OrdinaryDiffEqLowStorageRK`.

```julia
using EffectivelySymmetric, OrdinaryDiffEqLowStorageRK

sol = solve(prob, EES25_2N(); dt = 0.01)
```

### `CFEES25()`: Commutator-Free Lie-Group

Activated by loading `OrdinaryDiffEqLinear`, `ExponentialUtilities`, and `SciMLOperators`.

```julia
using EffectivelySymmetric
using OrdinaryDiffEqLinear, SciMLOperators, ExponentialUtilities

op   = MatrixOperator(A0; update_func! = update_A!)
prob = ODEProblem(op, u0, tspan)
sol  = solve(prob, CFEES25(); dt = 0.01)
```

## Method

The revised `EES25` tableau is

```
0   |
1/3 | 1/3
5/6 | -5/48   15/16
----+----------------
    | 1/10    1/2     2/5
```

The `EES27` tableau is

```
0 |
(2-√2)/3 | (2-√2)/3
(2+√2)/6 | (-4+√2)/24   (4+√2)/8
(4+√2)/6 | (-176+145√2)/168   3(8-5√2)/56   3(3-√2)/7
---------+----------------------------------------------------
         | (5-3√2)/14   (3+√2)/14   3(-1+2√2)/14   (9-4√2)/14
```

## License

Apache 2.0
