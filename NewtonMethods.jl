# Precompile the module
__precompile__(true)

module NewtonMethods

# Include common files
include("Jacobian.jl")
include("Options.jl")
include("Constants.jl")
include("Error.jl")
include("Helper.jl")
include("CheckOptionsNLEQ1.jl")

# Create a global container for the options since they need to be stored for
# further use
#commonWk = Dict{ASCIIString,OptionsNLEQ}()
#commonWk["NLEQ1"] = OptionsNLEQ()

global wk = OptionsNLEQ()

# Include the solver specific files
include("NLEQ1.jl")

# Export the required methods
export nleq1,OptionsNLEQ
export OPT_RTOL, OPT_QSUCC, OPT_MODE, OPT_JACGEN, OPT_JACFCN, OPT_MSTOR, OPT_ML,
        OPT_MU, OPT_ISCAL, OPT_PRINTWARNING, OPT_PRINTITERATION,
        OPT_PRINTIOWARN, OPT_PRINTIOMON, OPT_PRINTIOSOL, OPT_PRINTSOLUTION,
        OPT_NONLIN, OPT_QRANK1, OPT_QORDI, OPT_QSIMPL, OPT_NOROWSCAL,
        OPT_BOUNDEDDAMP, OPT_IORMON, OPT_NITMAX, OPT_FCBAND, OPT_SIGMA,
        OPT_SIGMA2, OPT_AJDEL, OPT_AJMIN, OPT_ETADIF, OPT_ETAINI, OPT_NBROY,
        OPT_FCSTART, OPT_FCMIN
end
