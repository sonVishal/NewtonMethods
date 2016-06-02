# String constants for NLEQ options

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
const OPT_PRINTIO               = "printIO"
const OPT_PRINTFILENAME         = "fileNameToPrintTo"
const OPT_PRINTFILEMODE         = "filePrintMode"
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

# Temporary Workspace variables
const WK_A                      = "a"
const WK_DXSAVE                 = "dxsave"
const WK_DX                     = "dx"
const WK_DXQ                    = "dxq"
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
const WK_FCMON                  = "fcmon"
const WK_FCA                    = "fca"
const WK_FCKEEP                 = "fckeep"
const WK_FCPRI                  = "fcpri"
const WK_DMYCOR                 = "dmycor"
const WK_SUMXS                  = "sumxs"

function getMachineConstants(i)
# 1 = realmin
# 2 = realmax
# 3 = eps
# 4 = 2.0*eps
# 5 = log10(2.0)
# 6 = sqrt(realmin/eps)
    if i < 1 || i > 6
        error("getMachineConstants - i out of bounds");
    end

    if i == 1
        return realmin(Float64);
    elseif i == 2
        return realmax(Float64);
    elseif i == 3
        return eps(Float64);
    elseif i == 4
        return 2.0*eps(Float64);
    elseif i == 5
        return log10(2.0);
    elseif i == 6
        return sqrt(getMachineConstants(1)/getMachineConstants(3));
    end
end
