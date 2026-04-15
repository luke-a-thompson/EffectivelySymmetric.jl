# SO(3) ODE — CG2 vs EffectivelySymmetric Lie-group integrators, plotted on the unit sphere.
#
# State  : R ∈ SO(3), stored as u = vec(R) ∈ ℝ⁹
# ODE    : dR/dt = Ω(R) · R,  ω(R) = ω₀ + R·b  (state-dependent)
# Methods: CG2 (built-in), CFEES25, and CFEES27 when available
#
# Run from the repo root:
#   julia --project=examples examples/so3_sphere.jl

using OrdinaryDiffEqLinear
using SciMLOperators
using LinearAlgebra
using CairoMakie
using Printf
using EffectivelySymmetric

function hat(v::Vector{Float64})::Matrix{Float64}
    return [  0.0   -v[3]  v[2]
              v[3]   0.0  -v[1]
             -v[2]   v[1]  0.0 ]
end

const ω₀ = [0.05, 0.30, 1.00]
const b   = [1.80, 0.60, 0.00]

function update_A!(A::Matrix{Float64}, u::Vector{Float64}, p, t)
    R = reshape(u, 3, 3)
    A .= kron(I(3), hat(ω₀ + R * b))
end

R₀    = Matrix{Float64}(I, 3, 3)
u₀    = vec(R₀)
A₀    = zeros(9, 9)
op    = MatrixOperator(A₀; update_func! = update_A!)
tspan = (0.0, 5.0)
prob  = ODEProblem(op, u₀, tspan)

function run_solver(alg, label)
    t0  = time()
    sol = solve(prob, alg; dt = 0.02, saveat = 0.04)
    elapsed = time() - t0

    so3_err = maximum(norm(reshape(u, 3, 3)' * reshape(u, 3, 3) - I(3)) for u in sol.u)

    @printf("%-8s  steps=%-4d  SO(3) err=%.2e  time=%.3f s\n",
            label, length(sol.t), so3_err, elapsed)
    return sol
end

spin_axis(sol) = reduce(hcat, [Float32.(reshape(u, 3, 3)[:, 3]) for u in sol.u])'

θ = range(0.0, π,   60)
φ = range(0.0, 2π, 120)
sx = Float32[sin(t) * cos(p) for t in θ, p in φ]
sy = Float32[sin(t) * sin(p) for t in θ, p in φ]
sz = Float32[cos(t)          for t in θ, p in φ]

BG = RGBf(0.07, 0.07, 0.10)

fig = Figure(size = (1100, 950), backgroundcolor = BG)

method_specs = [
    (label = "CG2", alg = CG2(), color = RGBAf(1.0, 0.45, 0.25, 0.9)),
    (label = "CFEES25", alg = CFEES25(), color = RGBAf(0.20, 0.85, 0.75, 0.9)),
]

if isdefined(EffectivelySymmetric, :CFEES27)
    push!(method_specs, (
        label = "CFEES27",
        alg = getfield(EffectivelySymmetric, :CFEES27)(),
        color = RGBAf(0.95, 0.78, 0.22, 0.9),
    ))
else
    println("CFEES27 is not available in the current EffectivelySymmetric build; skipping it.")
end

results = [(; spec..., sol = run_solver(spec.alg, spec.label)) for spec in method_specs]
trajectories = [spin_axis(result.sol) for result in results]
title_methods = join([result.label for result in results], " vs ")

ax = Axis3(
    fig[1, 1];
    aspect          = (1, 1, 1),
    title           = "SO(3) e₃ trajectory  ·  " * title_methods,
    titlecolor      = :white,
    backgroundcolor = BG,
    xgridcolor      = RGBAf(1, 1, 1, 0.08),
    ygridcolor      = RGBAf(1, 1, 1, 0.08),
    zgridcolor      = RGBAf(1, 1, 1, 0.08),
    xticklabelcolor = :gray55,
    yticklabelcolor = :gray55,
    zticklabelcolor = :gray55,
)

surface!(ax, sx, sy, sz;
    color        = fill(RGBAf(0.30, 0.50, 0.85, 0.15), size(sx)...),
    shading      = NoShading,
    transparency = true,
)
wireframe!(ax, sx, sy, sz;
    color = RGBAf(0.4, 0.6, 1.0, 0.06), linewidth = 0.3)

for (result, traj) in zip(results, trajectories)
    lines!(ax, traj[:, 1], traj[:, 2], traj[:, 3];
        color = result.color, linewidth = 2.2)
end

start_traj = first(trajectories)
scatter!(ax, [start_traj[1, 1]], [start_traj[1, 2]], [start_traj[1, 3]];
    color = :white, markersize = 14, strokewidth = 0)

for (result, traj) in zip(results, trajectories)
    scatter!(ax, [traj[end, 1]], [traj[end, 2]], [traj[end, 3]];
        color = RGBAf(result.color.r, result.color.g, result.color.b, 1.0),
        markersize = 12, strokewidth = 1, strokecolor = :white)
end

legend_elements = Any[
    LineElement(color = RGBAf(r.color.r, r.color.g, r.color.b, 1.0), linewidth = 3)
    for r in results
]
push!(legend_elements, MarkerElement(color = :white, marker = :circle, markersize = 10))

legend_labels = [r.label for r in results]
push!(legend_labels, "start (shared)")

Legend(fig[1, 2],
    legend_elements,
    legend_labels;
    backgroundcolor = RGBf(0.13, 0.13, 0.16),
    labelcolor = :white,
    framecolor = RGBAf(1, 1, 1, 0.1),
)

out = joinpath(@__DIR__, "..", "figures", "so3_cg2.png")
mkpath(dirname(out))
save(out, fig; px_per_unit = 2)
println("Saved → $out")
