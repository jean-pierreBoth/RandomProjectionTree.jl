using Printf
using Distances

function check_for_random_data()
    nbspectra = 20000
    nbval = 5000
    @debug "check_for_random_data: dispatching nbspecta : " nbspectra
    spectra = Array{Array{Float64,1},1}(undef, nbspectra)
    for i in 1:nbspectra
        spectra[i] = rand(nbval)
    end
    D=Jaccard()
    depth = 8
    argument = RPTreeArg(D, 8, 1.5)    
    rptree=RPTree(argument, spectra)
    @time leafCenters = randomProjection(rptree)
    leaves, diameters = analyzeSplittingInfo(rptree)
    #
    if length(leaves) != 2^depth
        @printf stdout "\n bad number of leaves , shoud be %d, got %d " 2^8  length(leaves)
        return false
    end
    #
    return true
end
