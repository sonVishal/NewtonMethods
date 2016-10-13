# TODO: Make everything work with Float64 as well as BigFloat
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

Currently the following two solvers are implemented:
- NLEQ1: Damped Newton-algorithm for systems of highly nonlinear
    equations - damping strategy due to Ref. (1).
- NLEQ2: Damped Newton-algorithm with rank strategy for systems of
    highly nonlinear equations - damping strategy due to Ref.(1).

## Licence
The MIT License (MIT)

Copyright (c) 2016 Vishal Sontakke

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Warranty
This code has been tested up to a certain level. Defects and
weaknesses, which may be included in the code, do not establish
any warranties by ZIB. ZIB does not take over any liabilities
which may follow from acquisition or application of this code.
"""
module NewtonMethods

using ForwardDiff

# Export the required methods
export nleq1, nleq2, OptionsNLEQ, clearWorkspace, clearAllWorkspaces

# Export the required options
export OPT_RTOL, OPT_QSUCC, OPT_MODE, OPT_JACGEN, OPT_JACFCN, OPT_MSTOR, OPT_ML,
    OPT_MU, OPT_ISCAL, OPT_PRINTWARNING, OPT_PRINTITERATION, OPT_PRINTIOWARN,
    OPT_PRINTIOMON, OPT_PRINTIOSOL, OPT_PRINTSOLUTION, OPT_NONLIN, OPT_QRANK1,
    OPT_QORDI, OPT_QSIMPL, OPT_NOROWSCAL, OPT_BOUNDEDDAMP, OPT_IORMON, OPT_NITMAX,
    OPT_FCBAND, OPT_SIGMA, OPT_SIGMA2, OPT_AJDEL, OPT_AJMIN, OPT_ETADIF, OPT_ETAINI,
    OPT_NBROY, OPT_FCSTART, OPT_FCMIN, OPT_IRANK, OPT_COND, OPT_STORE

# Export the required statistics
export STATS_NITER, STATS_NCORR, STATS_NFCN, STATS_NFCNJ, STATS_NJAC, STATS_NREJR1,
    STATS_XSCAL, STATS_RTOL, STATS_XITER, STATS_NATLEVEL, STATS_SIMLEVEL,
    STATS_STDLEVEL, STATS_PRECISION, STATS_DAMPINGFC, STATS_NEW, STATS_ICONV,
    STATS_CONV, STATS_SUMX, STATS_DLEVF, STATS_SUBCOND, STATS_SENS, STATS_IFAIL
# Include common files
include("Options.jl")
include("Constants.jl")
include("Common.jl")
include("Jacobian.jl")
include("SolverSpecific.jl")

global wkNLEQ1 = OptionsNLEQ()
global wkNLEQ2 = OptionsNLEQ()

"""
function clearWorkspace(name::AbstractString)

Function used to clear a workspace.

## Input parameter
`name` is the string which signifies which workspace is to be cleared.

For clearing the workspace related to NLEQ1 the following three strings are valid:
`NLEQ1`, `nleq1`, and `Nleq1`.

For clearing the workspace related to NLEQ2 the following three strings are valid:
`NLEQ2`, `nleq2`, and `Nleq2`
"""
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

"""
function clearWorkspace()

Function used to clear all workspaces.
"""
function clearWorkspace()
    empty!(wkNLEQ1.options)
    empty!(wkNLEQ2.options)
end

# Include the solver specific files
include("NLEQ1.jl")
include("NLEQ2.jl")

end
