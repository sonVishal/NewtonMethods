function n1int(n,fcn,jacFcn,x,xScal,opt,retCode,m1,m2,nBroy,fc,fcMin,mStor)
# Begin
# ------------------------------------------------------------------------------
# 1 Initialization
# ------------------------------------------------------------------------------
# 1.1 Control variables
qSucc       = getOption(opt,QPT_QSUCC,0) == 1
qScale      = ~(getOption(opt,OPT_QNSCAL,0) == 1)
qOrdi       = getOption(opt,OPT_QORDI,0) == 1
qSimpl      = getOption(opt,OPT_QSIMPL,0) == 1
qRank1      = getOption(opt,OPT_QRANK1,0) == 1
convMon     = getOption(opt,OPT_CONVMON,0)          # IORMON
iScal       = getOption(opt,OPT_ISCAL,0)
mode        = getOption(opt,OPT_MODE,0)
# jacFcn      = getOption(opt,OPT_JACOBIFCN)
qLU         = ~qSimpl

if convMon == 0
    convMon = 2
end
# ------------------------------------------------------------------------------
# 1.2 Derived dimensional parameters
if mStor == 0
    setOption!(opt,OPT_ML,0)
elseif mStor == 1
    setOption!(opt,OPT_ML,m1-m2)
    setOption!(opt,OPT_MU,m2-1-opt.options[OPT_ML])
end
# ------------------------------------------------------------------------------
# 1.3 Derived internal parameters
fcMin2  = fcMin*fcMin
rSmall  = sqrt(10.0*opt.options[OPT_RTOL])
# ------------------------------------------------------------------------------
# 1.4 Adaptation of input parameters, if necessary
if fc < fcMin
    fc = fcMin
end
if fc > 1.0
    fc = 1.0
end
# ------------------------------------------------------------------------------
# 1.5 Initial preparations
jacobiRefreshFlag   = false     #QJCRFR
stateLinSolver      = false     #QLINIT
internalRetCode     = 0         #IFAIL
fcBand              = 0.0
qIter               = true
if qBDamp == 1
    fcBand = inFcBand           #input value
end
# ------------------------------------------------------------------------------
# 1.5.1 Numerical differentiation related initializations
if jacFcn == numDiffJac
    aJdel = getOption(opt,OPT_AJDEL,0)
    if aJdel < small
        aJdel = sqrt(epMach*10.0)
    end
    aJmin = getOption(opt,OPT_AJMIN,0)
elseif jacFcn == numDiffJacFB
    etaDif = getOption(opt,OPT_ETADIF,0)
    if etaDif <= small
        etaDif = 1.0e-6
    end
    etaIni = getOption(opt,OPT_ETAINI,0)
    if etaIni <= small
        etaIni = 1.0e-6
    end
    epDiff = sqrt(epMach*10.0)
    etaMax = sqrt(epDiff)
    etaMin = epDiff*etaMax
end
# ------------------------------------------------------------------------------
# 1.5.2 Miscellaneous preparations of the first iteration step
if ~qSucc
    nIter   = 0
    nCorr   = 0
    nRejN   = 0     #NREJR1
    nFcn    = 0
    nJac    = 0
    nFcnJ   = 0
    qGenJ   = true
    qIniSc  = true
    fcKeep  = fc
    fcA     = fc
    fcPri   = fc
    fcK2    = fc
    fcMon   = fc
    conv    = 0.0

    if jacFcn == numDiffJacFB
        eta = etaIni*ones(n)
    end

    xa[:] = x[:]

    iConv   = 0
    alphaE  = 0.0
    sumAx1  = 0.0
    sumAx0  = 0.0
    cLin0   = 0.0
    qmStop  = false
# ------------------------------------------------------------------------------
# 1.6 Print monitor header
# ------------------------------------------------------------------------------
# 1.7 Startup step
# ------------------------------------------------------------------------------
# 1.7.1 Computation of residual vector
    res = zeros(x);
    fcn(x,res)
    nFcn += 1
    if length(res) ~= n
        throw(DimensionMismatch("Dimension of the function output does not
match the input dimension. Please check the function $fcn again."))
    end
else
    qIniSc = false
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Main iteration loop

# Repeat
while qIter
# ------------------------------------------------------------------------------
# 2 Startup of iteration step
    if ~jacobiRefreshFlag
# ------------------------------------------------------------------------------
# 2.1 Scaling of variables x(n)
    end

end

end
