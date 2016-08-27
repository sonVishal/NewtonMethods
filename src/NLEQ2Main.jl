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
    end














    # TODO: Remove this. For testing only
    retCode = 1
    return (x, xScal, retCode)
end
