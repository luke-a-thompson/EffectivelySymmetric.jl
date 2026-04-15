# EffectivelySymmetric.jl

Effectively-symmetric fixed-step integrators for the [SciML](https://sciml.ai) ecosystem.

The package currently provides:

- `EES25`: revised 3-stage explicit RK tableau.
- `EES27`: new 4-stage explicit RK tableau.
- `EES25_2N` and `CFEES25`: the existing low-storage and commutator-free EES25 variants.

EES25 is a fixed-step explicit method whose antisymmetric defect is $O(h^6)$. This makes it ideal for long-time structure-preserving integration and memory-efficient discrete adjoints via backward reconstruction (no checkpointing needed).

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

Activated by loading `OrdinaryDiffEqLinear`, `ExponentialUtilities`, and `SciMLOperators`. Solves linear ODEs of the form `du/dt = A(u,t) u` using matrix exponentials, preserving Lie-group structure (e.g. $\mathrm{SO}(n)$ orthogonality).

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

## Examples

`examples/so3_sphere.jl` solves a state-dependent $\mathrm{SO}(3)$ ODE and compares CG2 with CFEES25. Both traces agree on the forward trajectory; the difference is in the antisymmetric error order (CG2: 3, EES25: 6).

![SO(3) trajectory comparison: CG2 vs CFEES25](so3_cg2.png)

## Citation
If you use this package, please cite the paper _Explicit and Effectively Symmetric Runge-Kutta Methods_:
```bibtex
@misc{https://doi.org/10.48550/arxiv.2507.21006,
  doi = {10.48550/ARXIV.2507.21006},
  url = {https://arxiv.org/abs/2507.21006},
  author = {Shmelev,  Daniil and Ebrahimi-Fard,  Kurusch and Tapia,  Nikolas and Salvi,  Cristopher},
  keywords = {Numerical Analysis (math.NA),  Classical Analysis and ODEs (math.CA),  Rings and Algebras (math.RA),  FOS: Mathematics,  FOS: Mathematics,  16T05,  65L05,  65L06,  05C05},
  title = {Explicit and Effectively Symmetric Runge-Kutta Methods},
  publisher = {arXiv},
  year = {2025},
  copyright = {arXiv.org perpetual,  non-exclusive license}
}
```

## License

Apache 2.0
