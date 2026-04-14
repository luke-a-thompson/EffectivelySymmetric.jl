# bench.jl — BenchmarkTools harness for CG2 vs CFEES25
#
# Run from the repo root:
#   julia --project=benchmark benchmark/bench.jl
#
# Requires: OrdinaryDiffEqLinear, SciMLOperators, ExponentialUtilities,
#           BenchmarkTools (install into the benchmark environment).

using OrdinaryDiffEqLinear
using SciMLOperators
using LinearAlgebra
using BenchmarkTools
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

function make_prob()
    R₀ = Matrix{Float64}(I, 3, 3)
    u₀ = vec(R₀)
    A₀ = zeros(9, 9)
    op = MatrixOperator(A₀; update_func! = update_A!)
    return ODEProblem(op, u₀, (0.0, 5.0))
end

const PROB = make_prob()

solve_cg2()     = solve(PROB, CG2();     dt = 0.02, save_everystep = false, save_start = false, save_end = false)
solve_cfees25() = solve(PROB, CFEES25(); dt = 0.02, save_everystep = false, save_start = false, save_end = false)

function peak_live_delta_kib(f)::Int
    GC.gc(true)
    before = Base.gc_live_bytes()
    f()
    after = Base.gc_live_bytes()
    return div(max(after - before, 0), 1024)
end

print("Warming up CG2...     "); flush(stdout); solve_cg2();     println("done")
print("Warming up CFEES25... "); flush(stdout); solve_cfees25(); println("done")

println("\nRunning benchmarks (this may take ~30 s)…")

suite = BenchmarkGroup()
suite["CG2"]     = @benchmarkable solve_cg2()
suite["CFEES25"] = @benchmarkable solve_cfees25()

tune!(suite)
results = run(suite; verbose = false)

live_cg2     = peak_live_delta_kib(solve_cg2)
live_cfees25 = peak_live_delta_kib(solve_cfees25)

println("\n── Results ─────────────────────────────────────────────────────────")
for (name, t, live) in (("CG2", results["CG2"], live_cg2), ("CFEES25", results["CFEES25"], live_cfees25))
    med = median(t)
    mn  = minimum(t)
    mx  = maximum(t)
    @printf("%-8s  median=%6.2f ms  min=%6.2f ms  max=%6.2f ms  allocs=%d  peak_live=%+d KiB\n",
            name,
            time(med) / 1e6,
            time(mn)  / 1e6,
            time(mx)  / 1e6,
            allocs(med),
            live)
end

r = BenchmarkTools.ratio(median(results["CFEES25"]), median(results["CG2"]))
@printf("\nCFEES25 / CG2  time ratio = %.2f×\n", time(r))
