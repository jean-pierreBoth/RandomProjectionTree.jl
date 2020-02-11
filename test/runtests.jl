using Test
using RandomProjectionTree


using Logging
using Base.CoreLogging

logger = ConsoleLogger(stdout, CoreLogging.Debug)
global_logger(logger)


include("testrandom.jl")
@testset "random" begin
    @test check_for_random_data()
end



include("testpli.jl")