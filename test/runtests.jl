using Test
using GastroFlux

println("="^75)
println("Running GastroFlux Test Suite...")
println("="^75)

@testset "GastroFlux.jl Full Verification Suite" begin
    include("test_solvers.jl")
    include("test_solubility.jl")
    include("test_dissolution.jl")
    include("test_physiology.jl")
    include("test_acat_pbpk.jl")
    include("test_bcs_nca.jl")
end

println("="^75)
println("All GastroFlux tests passed successfully!")
println("="^75)
