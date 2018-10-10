filename="/home.2/Data/Spectro/LILLE/OrbiTrap/20120502_MASDA_LC20C37/recalibrated-20120503_MASDA_L20C37.imzML.pli"

include("/home/jpboth/Julia/src/pli.jl")

using MsPli

function check_for_pli(filename)
    exist = isfile(filename)

    if !exist
        @printf stdout "\n pli file not found %s" filename 
        return false
    end
    #
    spectra = MsPli.getAreasAsVector(pli)
    argument = RPTreeArg(D, 8, 1.3)
end



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
