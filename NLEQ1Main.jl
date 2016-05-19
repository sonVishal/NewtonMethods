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
    res = zeros(x)
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
        xW = n1scal(n,x,xa,xScal,iScal,qIniSc,opt)
        qIniSc = 0
        if nIter ~= 0
# ------------------------------------------------------------------------------
# 2.2 Aposteriori estimate of damping factor
            dxqa = dxq
            if ~qOrdi
                fcNump = sum((dx./xW).^2)
                th = fc - 1.0
                fcDnm = sum(((dxqa+th*dx)./xW).^2)
# ------------------------------------------------------------------------------
# 2.2.2 Decision criterion for Jacobian update technique
                qGenJ = true
                if fc == fcPri
                    qGenJ = fc < 1.0 || fcA < 1.0 || dMyCor <= fc*sigma || ~qRank1 || new + 2 > nBroy
                    fca = fc
                else
                    dMyCor = fcA*fcA*0.5*sqrt(fcNump/fcDnm)
                    if nonLin <= 3
                        fcCor = min(1.0,dMyCor)
                    else
                        fcCor = min(1.0,0.5*dMyCor)
                    end
                    fcA = max(min(fc,fcCor),fcMin)

                    #TODO: print something

                end
                fck2 = fcA
# ------------------------------------------------------------------------------
# 2.2.1 Computation of the numerator of damping factor predictor
                fcNmp2 = sum((dxqa./xW).^2)
                fcNump = fcNump*fcNmp2
            end
        end
    end
    jacobiRefreshFlag = false
# ------------------------------------------------------------------------------
# 2.3 Jacobian matrix
# ------------------------------------------------------------------------------
# 2.3.1 Jacobian generation by routine jac or difference approximation
# if qGenJ == true
# - or -
# Rank-1 update of Jacobian if qGenJ == false
    if qGenJ && (~qSimpl || nIter == 0)
        new = 0
        if jacGen == 1
            jac(x,J)
            a[:,:] = J[:,:]
        else
            if mStor == 0
                if jacGen == 3
                    #TODO: Jacobian with num diff feedback
                end
                if jacGen == 2
                    #TODO: Jacobian with num diff without feedback
                end
            else if mStor == 1
                if jacGen == 3
                    #TODO: banded num diff with feedback
                end
                if jacGen == 2
                    #TODO: banded num diff without feedback
                end
            end
        end
        nJac += 1

        #TODO: break if jacobian computation fails
        # Use try catch method

    else if ~qSimpl
        new += 1
    end

    if new == 0 && (qLU || nIter == 0)
# ------------------------------------------------------------------------------
# 2.3.2.1 Save scaling values
        xWa = xW[1:n]
        if issparse(a)
            nza = nnz(a)
            (row,col) = find(a~=0);
            for k = 1:nza
                a[row[k],col[k]] = -a[row[k],col[k]]*xW[col[k]]
            end
        else
# ------------------------------------------------------------------------------
# 2.3.2.2 Prepare Jacobian for use by band-solver
            if mStor == 1
                for l1 = 1:n
                    for l2 = m2:-1:1
                        a[l2+ml,l1] = a[l2,l1]
                    end
                end
            end

# ------------------------------------------------------------------------------
# 2.4 Prepare solution of the linear system
# ------------------------------------------------------------------------------
# 2.4.1 Internal column scaling of matrix A
            if mStor == 0
                for k = 1:n
                    a[1:n,k] = -a[1:n,k]*xW[k]
                end
            else if mStro == 1
                for k = 1:n
                    l2 = max(1+m2-k,ml+1)
                    l3 = max(n+m2-k,m1)
                    a[l2:l3,k] = -a[l2:l3,k]*xW[k]
                end
            end
        end
# ------------------------------------------------------------------------------
# 2.4.2 Row scaling of matrix A
        if qScale
            if mStor == 0
                (a,fW) = n1scrf(n,n,a)
            else if mStro == 1
                (a,fW) = n1scrb(n,m1,ml,mu,a)
            end
        else
            fW = ones(n);
        end
    end
# ------------------------------------------------------------------------------
# 2.4.3 Save and scale values of F(n)
    fA = f;
    t1 = f.*fW;


end
