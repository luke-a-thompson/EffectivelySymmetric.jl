using Test

@testset "EffectivelySymmetric.jl" begin
    include("order_test.jl")
    include("grad_order_test.jl")
end
