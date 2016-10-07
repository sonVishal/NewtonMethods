# TODO: Make everything work with Float64 as well as BigFloat
# TODO: Test all the Jacobian functions
# Precompile the module
__precompile__(true)
"""
## References:
 1. P. Deuflhard:
     Newton Methods for Nonlinear Problems. -
     Affine Invariance and Adaptive Algorithms.
     Series Computational Mathematics 35, Springer (2004)
 2. U. Nowak, L. Weimann:
     A Family of Newton Codes for Systems of Highly Nonlinear
     Equations - Algorithm, Implementation, Application.
     ZIB, Technical Report TR 90-10 (December 1990)

 Please generate the documentation using the following steps
"""
module NewtonMethods

using ForwardDiff

# Export the required methods
export nleq1, nleq2, n1int, n2int, OptionsNLEQ, clearWorkspace, clearAllWorkspaces
export OPT_RTOL, OPT_QSUCC, OPT_MODE, OPT_JACGEN, OPT_JACFCN, OPT_MSTOR, OPT_ML,
        OPT_MU, OPT_ISCAL, OPT_PRINTWARNING, OPT_PRINTITERATION,
        OPT_PRINTIOWARN, OPT_PRINTIOMON, OPT_PRINTIOSOL, OPT_PRINTSOLUTION,
        OPT_NONLIN, OPT_QRANK1, OPT_QORDI, OPT_QSIMPL, OPT_NOROWSCAL,
        OPT_BOUNDEDDAMP, OPT_IORMON, OPT_NITMAX, OPT_FCBAND, OPT_SIGMA,
        OPT_SIGMA2, OPT_AJDEL, OPT_AJMIN, OPT_ETADIF, OPT_ETAINI, OPT_NBROY,
        OPT_FCSTART, OPT_FCMIN, OPT_IRANK, OPT_COND, STATS_NITER, STATS_NCORR,
        STATS_NFCN, STATS_NFCNJ, STATS_NJAC, STATS_NREJR1, STATS_XSCAL,
        STATS_RTOL, STATS_XITER, STATS_NATLEVEL, STATS_SIMLEVEL, STATS_STDLEVEL,
        STATS_PRECISION, STATS_DAMPINGFC, STATS_NEW, STATS_ICONV, STATS_CONV,
        STATS_SUMX, STATS_DLEVF

# Include common files
include("Jacobian.jl")
include("Options.jl")
include("Constants.jl")
include("SolverSpecific.jl")
include("Common.jl")

# Machine Constants
global epMach = 1e-17
global small  = 1e-150
global great  = 1.0/small

global wkNLEQ1 = OptionsNLEQ()
global wkNLEQ2 = OptionsNLEQ()

function clearWorkspace(name::AbstractString)
    if name == "NLEQ1" || name == "nleq1" || name == "Nleq1"
        empty!(wkNLEQ1.options)
    elseif name == "NLEQ2" || name == "nleq2" || name == "Nleq2"
        empty!(wkNLEQ2.options)
    else
        println("Invalid option. Please specify the correct argument.")
        println("\"Nleq1\" or \"NLEQ1\" or \"nleq1\" - to clear the workspace for nleq1 function")
        println("\"Nleq2\" or \"NLEQ2\" or \"nleq2\" - to clear the workspace for nleq2 function")
    end
end

function clearAllWorkspaces()
    empty!(wkNLEQ1.options)
    empty!(wkNLEQ2.options)
end

# Include the solver specific files
include("NLEQ1.jl")
include("NLEQ2.jl")

end
