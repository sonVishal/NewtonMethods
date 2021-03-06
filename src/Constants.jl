# AbstractString constants for NLEQ options
const OPT_RTOL                  = "relativeTolerance"
const OPT_QSUCC                 = "successiveCall"
const OPT_MODE                  = "mode"
const OPT_JACGEN                = "jacobianMethod"
const OPT_JACFCN                = "jacobianFunction"
const OPT_MSTOR                 = "jacobiStorage"
const OPT_ML                    = "lowerBandwidth"
const OPT_MU                    = "upperBandwidth"
const OPT_ISCAL                 = "scaling"
const OPT_PRINTWARNING          = "printWarning"
const OPT_PRINTITERATION        = "printIterationMonitor"
const OPT_PRINTIOWARN           = "printIOwarning"
const OPT_PRINTIOMON            = "printIOmonitor"
const OPT_PRINTIOSOL            = "printIOsolution"
const OPT_PRINTSOLUTION         = "printSolution"
const OPT_NONLIN                = "nonlinearType"
const OPT_QRANK1                = "rank1BroydenUpdates"
const OPT_QORDI                 = "ordinaryNewton"
const OPT_QSIMPL                = "simpleNewton"
const OPT_NOROWSCAL             = "automaticRowScaling"
const OPT_BOUNDEDDAMP           = "boundedDampingStrategy"
const OPT_IORMON                = "convergenceOrderMonitor"
const OPT_NITMAX                = "maxIter"
const OPT_FCBAND                = "boundedDampingFactor"
const OPT_SIGMA                 = "broydenDecisionParameter"
const OPT_SIGMA2                = "dampingFactorParameter"
const OPT_AJDEL                 = "relativePerturbation"
const OPT_AJMIN                 = "thresholdValue"
const OPT_ETADIF                = "targetRelativePerturbation"
const OPT_ETAINI                = "initialDenominatorDifference"
const OPT_NBROY                 = "maxBroydenSteps"
const OPT_FCSTART               = "initialDampingFactor"
const OPT_FCMIN                 = "minimumDampingFactor"
const OPT_IRANK                 = "initialRank"
const OPT_COND                  = "subConditionNumber"
const OPT_STORE                 = "storeIntermediateValues"

# Statistics
const STATS_NITER               = "numberOfNewtonIterations"
const STATS_NCORR               = "numberOfCorrectorSteps"
const STATS_NFCN                = "numberOfFunctionEvals"
const STATS_NFCNJ               = "numberOfFunctionEvalsForJac"
const STATS_NJAC                = "numberOfJacobianCalls"
const STATS_NREJR1              = "numRejectedRank1NewtonSteps"
const STATS_XSCAL               = "lastScalingVector"
const STATS_RTOL                = "finalRelativeAccuracy"
const STATS_XITER               = "allNewtonIterates"
const STATS_NATLEVEL            = "seqNaturalLevels"
const STATS_SIMLEVEL            = "seqNaturalLevelsSimpl"
const STATS_STDLEVEL            = "seqStandardLevels"
const STATS_PRECISION           = "achievedPrecision"
const STATS_DAMPINGFC           = "seqDampingFactors"
const STATS_NEW                 = "consecutiveRank1Updates"
const STATS_ICONV               = "convergenceMonitorState"
const STATS_CONV                = "conv"
const STATS_SUMX                = "sumx"
const STATS_DLEVF               = "dlevf"
const STATS_SUBCOND             = "subcondition"
const STATS_SENS                = "sensitivity"
const STATS_IFAIL               = "failureIndicator"

# Temporary Workspace variables
const WK_A                      = "a"
const WK_DXSAVE                 = "dxsave"
const WK_QA_DXSAVE              = "qa_dxsave"
const WK_DX                     = "dx"
const WK_DXQ                    = "dxq"
const WK_QA                     = "qa"
const WK_XA                     = "xa"
const WK_XWA                    = "xwa"
const WK_F                      = "f"
const WK_FA                     = "fa"
const WK_ETA                    = "eta"
const WK_XW                     = "xw"
const WK_FW                     = "fw"
const WK_DXQA                   = "dxqa"
const WK_SUMXA0                 = "sumxa0"
const WK_SUMXA1                 = "sumxa1"
const WK_QU                     = "qu"
const WK_T1                     = "t1"
const WK_T2                     = "t2"
const WK_T3                     = "t3"
const WK_FCMON                  = "fcmon"
const WK_FCA                    = "fca"
const WK_FCKEEP                 = "fckeep"
const WK_FCPRI                  = "fcpri"
const WK_DMYCOR                 = "dmycor"
const WK_SUMXS                  = "sumxs"
const WK_SENS1                  = "sens1"

"""
function getMachineConstants(T::DataType)

Get the machine constants depending on the DataType T.

Supported DataTypes are Float64 and BigFloat.
"""
function getMachineConstants(T::DataType)
    if T != Float64 && T != DataType
        println("ERROR: Wrong DataType provided. You provided $T, expected Float64 or BigFloat.")
        error("getMachineConstants - wrong DataType provided");
    end

    epMach = eps(T)
    small  = sqrt(realmin(T)/eps(T))
    great  = sqrt(realmax(T)/10.0)

    return (epMach, small, great)
end
