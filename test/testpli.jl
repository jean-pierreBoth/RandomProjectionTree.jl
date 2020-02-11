pliname= "LILLE/OrbiTrap/20120502_MASDA_LC20C37/recalibrated-20120503_MASDA_L20C37.imzML.pli"

dirname = "/home.2/Data/Spectro/"

filename = dirname*dirname

include("/home/jpboth/Julia/src/pli.jl")

using MsPli

function check_for_pli(filename)
    exist = isfile(filename)

    if !exist
        @printf stdout "\n pli file not found %s" filename 
        return true
    end
    #
    spectra = MsPli.getAreasAsVector(pli)
    argument = RPTreeArg(D, 8, 1.3)
    true
end
