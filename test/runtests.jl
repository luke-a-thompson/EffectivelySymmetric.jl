using Test

@testset "EffectivelySymmetric.jl" begin
    @testset "Classical convergence (solve)" begin
        include("classical_order_test.jl")
    end
    @testset "Time-reversal defect (solve)" begin
        include("reversal_order_test.jl")
    end
    @testset "CFEES25 Lie-group (solve)" begin
        include("cfees25_test.jl")
    end
    @testset "Aqua quality checks" begin
        include("aqua_test.jl")
    end
end
