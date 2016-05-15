# String constants for NLEQ options

const OPT_RTOL                  = "relativeTolerance"
const OPT_QSUCC                 = "successiveCall"
const OPT_MODE                  = "mode"
const OPT_JACOBIFCN             = "jacobianMethod"
const OPT_MSTOR                 = "jacobiStorage"
const OPT_ML                    = "lowerBandwidth"
const OPT_MU                    = "upperBandwidth"
const OPT_ISCAL                 = "scaling"
const OPT_PRINTWARNING          = "printWarning"
const OPT_PRINTIO               = "printIO"
const OPT_NONLIN                = "nonlinearType"
const OPT_QRANK1                = "rank1BroydenUpdates"
const OPT_QORDI                 = "ordinaryNewton"
const OPT_QSIMPL                = "simpleNewton"
const OPT_QNSCAL                = "automaticRowScaling"
const OPT_IBDAMP                = "boundedDampingStrategy"
const OPT_CONVMON               = "convergenceMonitor"
const OPT_NITMAX                = "maxIter"
const OPT_FCBAND                = "boundedDampingFactor"
const OPT_SIGMA                 = "broydenDecisionParameter"
const OPT_SIGMA2                = "dampingFactorParameter"
const OPT_AJDEL                 = "relativePerturbation"
const OPT_AJMIN                 = "thresholdValue"
const OPT_ETADIF                = "targetRelativePerturbation"
const OPT_ETAINI                = "initialDenominatorDifference"
const OPT_NBROY                 = "maxBroydenSteps"
const OPT_NOROWSCAL             = "rowScaling"
const OPT_FCSTART               = "initialDampingFactor"
const OPT_FCMIN                 = "minimumDampingFactor"
const OPT_ICONV                 = "currentConvMonitor"

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
