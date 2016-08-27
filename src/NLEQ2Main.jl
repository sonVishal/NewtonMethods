function n2int(n, fcn, x, xScal, rTol, nItmax, nonLin, iRank, cond, opt, retCode,
    m1, m2, nBroy, xIter, sumXall, dLevFall, sumXQall,tolAll, fcAll, fc, fcMin,
    sigma, sigma2, mStor, printWarn, printMon, printSol, printIOwarn, printIOmon,
    printIOsol, qBDamp)
    # --------------------------------------------------------------------------
    # Since wkNLEQ2 is module global
    # Create the local variables here rather than taking them as arguments
    a       = wkNLEQ2.options[WK_A]
    qa      = wkNLEQ2.options[WK_DXSAVE]
    dxSave  = wkNLEQ2.options[WK_DXSAVE]
    dx      = wkNLEQ2.options[WK_DX]
    dxQ     = wkNLEQ2.options[WK_DXQ]
    xa      = wkNLEQ2.options[WK_XA]
    xwa     = wkNLEQ2.options[WK_XWA]
    f       = wkNLEQ2.options[WK_F]
    fa      = wkNLEQ2.options[WK_FA]
    eta     = wkNLEQ2.options[WK_ETA]
    xw      = wkNLEQ2.options[WK_XW]
    fw      = wkNLEQ2.options[WK_FW]
    dxQa    = wkNLEQ2.options[WK_DXQA]
    qu      = wkNLEQ2.options[WK_QU]
    t1      = wkNLEQ2.options[WK_T1]
    t2      = wkNLEQ2.options[WK_T2]
    fcA     = wkNLEQ2.options[WK_FCA]
    fcKeep  = wkNLEQ2.options[WK_FCKEEP]
    fcPri   = wkNLEQ2.options[WK_FCPRI]
    dMyCor  = wkNLEQ2.options[WK_DMYCOR]
    conv    = wkNLEQ2.options[STATS_CONV]
    dLevF   = wkNLEQ2.options[STATS_DLEVF]
    nIter   = wkNLEQ2.options[STATS_NITER]
    nCorr   = wkNLEQ2.options[STATS_NCORR]
    nFcn    = wkNLEQ2.options[STATS_NFCN]
    nFcnJ   = wkNLEQ2.options[STATS_NFCNJ]
    nJac    = wkNLEQ2.options[STATS_NJAC]
    nRejR1  = wkNLEQ2.options[STATS_NREJR1]
    newt    = wkNLEQ2.options[STATS_NEW]

    # --------------------------------------------------------------------------
    # 0.1 Variables that need to be defined before since they appear in different
    # scopes. The declaration and usage are in different scopes.

    # --------------------------------------------------------------------------
    # 0.2 Persistent variables
    cLin0   = getOption!(wkNLEQ1,"P_CLIN0",0.0)
    alphaE  = getOption!(wkNLEQ1,"P_ALPHAE",0.0)
    # --------------------------------------------------------------------------

    epMach  = getMachineConstants(3)
    small   = getMachineConstants(6)
    n       = length(x)
    # Begin
    # --------------------------------------------------------------------------
    # 1 Initialization
    # --------------------------------------------------------------------------
    qSucc       = Bool(opt.options[OPT_QSUCC])
    qScale      = opt.options[OPT_NOROWSCAL] != 1
    qRank1      = Bool(opt.options[OPT_QRANK1])
    iOrMon      = getOption!(opt, OPT_IORMON, 0)
    if iOrMon == 0
        iOrMon = 2
    end
    iScal       = getOption!(opt, OPT_ISCAL, 0)
    mode        = getOption!(opt, OPT_MODE,  0)
    jacGen      = opt.options[OPT_JACGEN]
    qMixIO      = typeof(printIOmon) == typeof(printIOsol)
    if qMixIO && typeof(printIOmon) == IOStream && printIOmon.name != printIOsol.name
        qMixIO = false
    end
    qMixIO      &= mPrMon != 0 && mPrSol != 0
    # --------------------------------------------------------------------------
    # 1.2 Derived dimensional parameters
    if mStor == 0
        ml = 0
        mu = 0
    elseif mStor == 1
        ml = m1 - m2
        mu = m2 - 1 - ml
    end
    minRnk = max(1,n-max(round(Int,n/10.0),10))
    # --------------------------------------------------------------------------
    # 1.3 Derived internal parameters
    fcMin2  = fcMin*fcMin
    fcMinH  = sqrt(fcMin)
    tolMin  = sqrt(10.0*epMach)
    rSmall  = sqrt(10.0*rTol)
    # --------------------------------------------------------------------------
    # 1.4 Adaptation of input parameters, if necessary
    if fc < fcMin
        fc = fcMin
    end
    if fc > 1.0
        fc = 1.0
    end
    # --------------------------------------------------------------------------
    # 1.5 Initial preparations
    qJcRfr              = false
    qLInit              = false
    qRepeat             = false
    iFail               = 0
    fcBand              = 0.0
    if qBDamp
        fcBand = opt.options[OPT_FCBAND]
    end
    # --------------------------------------------------------------------------
    # 1.5.1 Numerical differentiation related initializations
    if jacGen == 2
        aJdel = getOption!(opt, OPT_AJDEL, 0.0)
        if aJdel <= small
            aJdel = sqrt(epMach*10.0)
        end
        aJmin = getOption!(opt, OPT_AJMIN, 0.0)
    elseif jacGen == 3
        etaDif = getOption!(opt, OPT_ETADIF, 0.0)
        if etaDif <= small
            etaDif = 1.0e-6
        end
        etaIni = getOption!(opt, OPT_ETAINI, 0.0)
        if etaIni <= small
            etaIni = 1.0e-6
        end
        epDiff = sqrt(epMach*10.0)
        etaMax = sqrt(epDiff)
        etaMin = epDiff*etaMax
    end
    # --------------------------------------------------------------------------
    # 1.5.2 Miscellaneous preparations of the first iteration step
    if !qSucc
        nIter  = 0
        nCorr  = 0
        nRejR1 = 0
        nFcn   = 0
        nFcnJ  = 0
        nJac   = 0
        qGenJ   = true
        qIniSc  = true
        fcKeep  = fc
        fcA     = fc
        fcPri   = fc
        fcMon   = fc
        fcK2    = fc
        conv   = 0.0
        if jacGen == 3
            eta = etaIni*ones(n)
        end

        xa[:] = x[:]

        iConv  = 0

        sumxa0  = 0.0
        sumxa1  = 0.0

        push!(xIter,x)
        wkNLEQ1.options[STATS_NITER] = nIter

        qMStop = false
        # ----------------------------------------------------------------------
        # 1.6 Print monitor header
        if mPrMon >= 2 && !qMixIO
            write(printIOmon,
            "\n",
            "  ******************************************************************",
            "\n",
            "        It       Normf           Normx         Damp.Fct.   New\n")
        end
        # ----------------------------------------------------------------------
        # 1.7 Startup step
        # ----------------------------------------------------------------------
        # 1.7.1 Computation of residual vector
        try
            fcn(f,x)
            nFcn += 1
            iFail = 0
        catch
            iFail = -1
            retCode = 82
        end
        # TODO: Perform check whether f and x are of same length in CheckOptionsNLEQ
    else
        qIniSc = false
    end
    # --------------------------------------------------------------------------
    # Main iteration loop

    # Repeat
    while qIter
        # ----------------------------------------------------------------------
        # 2 Startup of iteration step
        if !qJcRfr
            # ------------------------------------------------------------------
            # 2.1 Scaling of variables x(n)
            xw = n1scal(n , x, xa, xScal, iScal, qIniSc, opt)
            qIniSc = false
            if nIter != 0
                # Preliminary psuedo-rank
                iRankA = iRank
                # --------------------------------------------------------------
                # 2.2 Aposteriori estimate of damping factor
                dxQa = dxQ
                fcNumP = sum((dx./xw).^2)
                th = fc - 1.0
                fcDnm = sum(((dxQa+th*dx)./xw).^2)
                # ----------------------------------------------------------
                # 2.2.2 Decision criterion for Jacobian update technique
                # qGenJ == true   numerical differentation,
                # qGenJ == false  rank1 updating
                qGenJ = true
                if fc == fcPri
                    qGenJ = fc < 1.0 || fcA < 1.0 || dMyCor <= fc*sigma ||
                    !qRank1 || newt + 2 > nBroy

                    fcA = fc
                else
                    dMyCor = fcA*fcA*0.5*sqrt(fcNumP/fcDnm)
                    if nonLin <= 3
                        fcCor = min(1.0,dMyCor)
                    else
                        fcCor = min(1.0,0.5*dMyCor)
                    end
                    fcA = max(min(fc,fcCor),fcMin)

                    if mPrMon >= 5
                        write(printIOmon,"\n",
                        @sprintf("+++  aposteriori estimate  +++\n"),
                        @sprintf(" fcCor  = %18.10e\n",fcCor),
                        @sprintf(" fc     = %18.10e\n",fc),
                        @sprintf(" dMyCor = %18.10e\n",dMyCor),
                        @sprintf(" fcNumP = %18.10e\n",fcNumP),
                        @sprintf(" fcDnm  = %18.10e\n",fcDnm),
                        @sprintf("++++++++++++++++++++++++++++++\n"))
                    end
                end
                fck2 = fcA
                # ----------------------------------------------------------
                # 2.2.1 Computation of the numerator of damping
                # factor predictor
                fcNmp2 = sum((dxQa./xw).^2)
                fcNumP = fcNumP*fcNmp2
            end
        end
        qJcRfr = false
        # ----------------------------------------------------------------------
        # 2.3 Jacobian matrix
        # TODO: Multiple dispatch for the Jacobian evaluation
        # ----------------------------------------------------------------------
        # 2.3.1 Jacobian generation by routine jac or difference approximation
        # if qGenJ == true
        # - or -
        # Rank-1 update of Jacobian if qGenJ == false
        if qGenJ && (!qSimpl || nIter == 0)
            newt = 0
            if jacGen == 1
                jac = getOption(opt, OPT_JACFCN, 0)
                try
                    jac(a,x)
                    # iFail = 0
                catch err
                    retCode = 82
                    iFail   = -1
                    # throw(err)
                end
            else
                if mStor == 0
                    if jacGen == 3
                        (a,eta,nFcnJ,iFail) = n1jcf(fcn,n,n,x,f,xw,eta,etaMin,etaMax,etaDif,conv,nFcnJ)
                    end
                    if jacGen == 2
                        (a,nFcnJ,iFail) = n1jac(fcn,n,n,x,f,xw,aJdel,aJmin,nFcnJ)
                    end
                    if jacGen == 4
                        (a,iFail) = n1jacFAD(fcn,x)
                    end
                elseif mStor == 1
                    if jacGen == 3
                        (a,eta,nFcnJ,iFail) = n1jcfb(fcn,n,m1,ml,x,f,xw,eta,etaMin,etaMax,etaDif,conv,nFcnJ)
                    end
                    if jacGen == 2
                        (a,nFcnJ,iFail) = n1jacb(fcn,n,m1,ml,x,f,xw,aJdel,aJmin,nFcnJ)
                    end
                    if jacGen == 4
                        (a,iFail) = n1jacFAD(fcn,x)
                    end
                end
            end
            nJac += 1
            if jacGen == 1 && iFail < 0
                retCode = 83
                break
            end

            if jacGen != 1 && iFail != 0
                retCode = 82
                break
            end
        else
            newt += 1
        end

        if newt == 0
            # ------------------------------------------------------------------
            # 2.3.2 Save scaling values
            xwa = xw[1:n]
            # ------------------------------------------------------------------
            # 2.4 Prepare solution of the linear system
            # --------------------------------------------------------------
            # 2.4.1 Internal column scaling of matrix A
            if mStor == 0
                for k = 1:n
                    a[1:n,k] = -a[1:n,k]*xw[k]
                end
            elseif mStor == 1
                for k = 1:n
                    l2 = max(1+m2-k,ml+1)
                    l3 = min(n+m2-k,m1)
                    a[l2:l3,k] = -a[l2:l3,k]*xw[k]
                end
            end
            # ------------------------------------------------------------------
            # 2.4.2 Row scaling of matrix A
            if qScale
                if mStor == 0
                    (a,fw) = n1scrf(n,n,a)
                elseif mStor == 1
                    (a,fw) = n1scrb(n,m1,ml,mu,a)
                end
            else
                fw = ones(n)
            end
        end
        # ----------------------------------------------------------------------
        # 2.4.3 Save and scale values of F(n)
        fa[:] = f
        t1 = f.*fw
        iRankA = iRank
        if nIter != 0
            iRank = n
        end
        qNext = false
        # ----------------------------------------------------------------------
        # 3 Central part of iteration step
        # Pseudo-rank reduction loop
        # ==========================
        while
            # ----------------------------------------------------------------------
            # 3.1 Solution of the linear system
            # ----------------------------------------------------------------------
            # 3.1.1 Decomposition of (n,n) matrix A
        end
    end













    # TODO: Remove this. For testing only
    retCode = 1
    return (x, xScal, retCode)
end
