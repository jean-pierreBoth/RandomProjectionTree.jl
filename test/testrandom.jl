
function check_for_random_data()
    nbspectra = 10000
    nbval = 5000
    spectra = Array{Array{Float64,1},1}(undef, nbspectra)
    for i in 1:nbpsectra
        spectra[i] = rand(nbval)
    end
    D=Jaccard()
    depth = 8
    argument = RPTreeArg(D, 8, 1.5)    
    rptree=RPTree(argument, spectra)
    leaves, diameters = analyzeSplittingInfo(rptree)
    #
    if size(leaves) != 2^depth
        @printf stdout "\n bad number of leaves , shoud be %d, got %d " 2^8  size(leaves)
        return false
    end
    #
    return true
end
