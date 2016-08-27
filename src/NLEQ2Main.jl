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

    # TODO: Remove this. For testing only
    retCode = 1
    return (x, xScal, retCode)
end
