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
if qBDamp == 1
    fcBand = inFcBand           #input value
end
# ------------------------------------------------------------------------------
# 1.5.1 Numerical differentiation related initializations
if jacGen == 

end
