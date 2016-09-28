using Debug
function nleq2(fcn, x::Vector{Float64}, xScal::Vector{Float64}, opt::OptionsNLEQ)

    # Initialize a common message string variable
    message = ""

    # Initialize error code 0
    retCode = 0

#-------------------------------------------------------------------------------
# Printing related stuff
#-------------------------------------------------------------------------------
    # Print warning messages?
    printWarn   = getOption(opt,OPT_PRINTWARNING,0)
    # Print iteration summary?
    printMon    = getOption!(opt,OPT_PRINTITERATION,0)
    # Print solution summary?
    printSol    = getOption!(opt,OPT_PRINTSOLUTION,0)
    # Where to print?
    # Defaults to STDOUT
    printIOwarn = getOption(opt,OPT_PRINTIOWARN,STDOUT)
    printIOmon  = getOption!(opt,OPT_PRINTIOMON,STDOUT)
    printIOsol  = getOption!(opt,OPT_PRINTIOSOL,STDOUT)
#-------------------------------------------------------------------------------

    # First call or successive call
    qSucc   = Bool(getOption!(opt,OPT_QSUCC,0))
    qIniMon = (printMon >= 1 && !qSucc)

    # TODO: Improve checkOptions and handle the errors properly!!

    # Check input parameters and options
    n = length(x)
    retCode = checkOptions(n,x,xScal,opt)

    # Exit if any parameter error was detected
    if retCode != 0
        error("Exit with return code $retCode")
    end

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    m1 = n
    m2 = n

    jacGen = getOption!(opt,OPT_JACGEN,0)
    if jacGen == 0
        jacGen = 2
        setOption!(opt, OPT_JACGEN, jacGen)
    end

    qRank1 = Bool(getOption!(opt, OPT_QRANK1, 0))

    if qRank1
        nBroy = getOption!(opt,OPT_NBROY,0)
        if nBroy == 0
            nBroy = max(m2,10)
            setOption!(opt,OPT_NBROY, nBroy)
        end
    else
        nBroy = 0
    end

    xIter    = getOption!(wkNLEQ2, "P_XITER", [])
    sumXall  = getOption!(wkNLEQ2, "P_SUMXALL", [])
    dLevFall = getOption!(wkNLEQ2, "P_DLEVFALL", [])
    sumXQall = getOption!(wkNLEQ2, "P_SUMXQALL", [])
    tolAll   = getOption!(wkNLEQ2, "P_TOLALL", [])
    fcAll    = getOption!(wkNLEQ2, "P_FCALL", [])
    # Check if this is a first call or successive call to nleq1
    # If first call then reset the workspace and persistent variables
    if !qSucc
        empty!(wkNLEQ2.options)
        initializeOptions(opt, wkNLEQ2, n, m1, nBroy, qRank1)
    end

    # Check for non linear option
    nonLin = getOption!(opt, OPT_NONLIN, 3)
    initOption!(opt, OPT_BOUNDEDDAMP, 0)

    if opt.options[OPT_BOUNDEDDAMP] == 0
        qBDamp = nonLin == 4
    elseif opt.options[OPT_BOUNDEDDAMP] == 1
        qBDamp = true
    elseif opt.options[OPT_BOUNDEDDAMP] == 2
        qBDamp = false
    end

    # Initialize bounded damping strategy restriction factor
    initOption!(opt, OPT_FCBAND, 0.0)
    if qBDamp
        if opt.options[OPT_FCBAND] < 1.0
            setOption!(opt, OPT_FCBAND, 10.0)
        end
    end

    # Maximum permitted number of iteration steps
    nItmax = getOption!(opt, OPT_NITMAX, 50)
    if nItmax <= 0
        nItmax = 50
        setOption!(opt, OPT_NITMAX, nItmax)
    end

    if qIniMon
        printInitialization(n, printIOmon, opt.options[OPT_RTOL], jacGen, 0,
        0, 0, opt.options[OPT_NOROWSCAL], qRank1, nonLin, qBDamp,
        opt.options[OPT_FCBAND], false, false, nItmax)
    end

    # Initial damping factor for highly nonlinear problems
    initOption!(opt, OPT_FCSTART, 0.0)
    qFcStart = opt.options[OPT_FCSTART] > 0.0
    if !qFcStart
        setOption!(opt, OPT_FCSTART, 1.0e-2)
        if nonLin == 4
            setOption!(opt, OPT_FCSTART, 1.0e-4)
        end
    end

    # Minimal permitted damping factor
    initOption!(opt,OPT_FCMIN,0.0)
    if opt.options[OPT_FCMIN] <= 0.0
        setOption!(opt, OPT_FCMIN, 1.0e-4)
        if nonLin == 4
            setOption!(opt, OPT_FCMIN, 1.0e-8)
        end
    end
    fcMin = getOption(opt,OPT_FCMIN,0.0)

    # Rank1 decision parameter SIGMA
    initOption!(opt,OPT_SIGMA,0.0)
    if opt.options[OPT_SIGMA] < 1.0
        setOption!(opt, OPT_SIGMA, 3.0)
    end
    if !qRank1
        setOption!(opt, OPT_SIGMA, 10.0/fcMin)
    end

    # Decision parameter about increasing too small predictor
    # to greater corrector value
    initOption!(opt,OPT_SIGMA2,0.0)
    if opt.options[OPT_SIGMA2] < 1.0
        setOption!(opt, OPT_SIGMA2, 10.0/fcMin)
    end

    # Starting value of damping factor (fcMin <= fc <= 1.0)
    if nonLin <= 2 && !qFcStart
        # for linear or mildly nonlinear problems
        fc = 1.0
    else
        # for highly or extremely nonlinear problems
        fc = getOption(opt, OPT_FCSTART, 0.0)
    end

    iRank = getOption!(opt, OPT_IRANK, 0)
    if iRank <= 0 || iRank > n
        iRank = n
        setOption!(opt, OPT_IRANK, iRank)
    end

    epMach = 1e-17
    cond = getOption!(opt, OPT_COND, 1.0/epMach)
    if cond < 1.0
        cond = 1.0/epMach
        setOption!(opt, OPT_COND, cond)
    end

    # Set starting damping factor
    setOption!(opt, OPT_FCSTART, fc)

    if printMon >= 2 && !qSucc
        write(printIOmon,"\nINFO: ","Internal parameters:",
        "\n\tStarting value for damping factor ",
        @sprintf("OPT_FCSTART\t= %1.2e",opt.options[OPT_FCSTART]),
        @sprintf("\n\tMinimum allowed damping factor OPT_FCMIN\t= %1.2e",fcMin),
        "\n\tRank-1 updates decision parameter ",
        @sprintf("OPT_SIGMA\t= %1.2e\n",opt.options[OPT_SIGMA]),
        @sprintf("\n\tInitial Jacobian pseudo-rank iRank\t\t= %1.2e", iRank),
        @sprintf("\n\tMaximum permitted subcondition cond\t\t= %1.2e\n", cond))
    end

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    retCode = -1

    if nBroy == 0
        nBroy = 1
    end

    # Call to n2int
    (x, xScal, retCode) = n2int(n, fcn, x, xScal, opt.options[OPT_RTOL], nItmax,
    nonLin, iRank, cond, opt, retCode, m1, m2, nBroy, xIter, sumXall, dLevFall,
    sumXQall, tolAll, fcAll, opt.options[OPT_FCSTART], opt.options[OPT_FCMIN],
    opt.options[OPT_SIGMA], opt.options[OPT_SIGMA2], printWarn, printMon,
    printSol, printIOwarn, printIOmon, printIOsol, qBDamp)

    # set stats variable
    stats = Dict{AbstractString,Any}()
    stats[STATS_XSCAL] = xScal
    if retCode == -1
        stats[STATS_RTOL] = tolAll[wkNLEQ2.options[STATS_NITER]]
    else
        stats[STATS_RTOL] = opt.options[OPT_RTOL]
    end
    stats[STATS_XITER]      = xIter
    stats[STATS_NATLEVEL]   = sumXall
    stats[STATS_SIMLEVEL]   = sumXQall
    stats[STATS_STDLEVEL]   = dLevFall
    stats[STATS_PRECISION]  = tolAll
    stats[STATS_DAMPINGFC]  = fcAll
    stats[STATS_NITER]      = wkNLEQ2.options[STATS_NITER]
    stats[STATS_NCORR]      = wkNLEQ2.options[STATS_NCORR]
    stats[STATS_NREJR1]     = wkNLEQ2.options[STATS_NREJR1]
    stats[STATS_NJAC]       = wkNLEQ2.options[STATS_NJAC]
    stats[STATS_NFCN]       = wkNLEQ2.options[STATS_NFCN]
    stats[STATS_NFCNJ]      = wkNLEQ2.options[STATS_NFCNJ]

    # Print statistics
    if printMon >= 2 && retCode != -1 && retCode != 10
        printStats(stats, printIOmon)
    end

    # Assign the persistent variables back
    setOption!(wkNLEQ2, "P_XITER", xIter)
    setOption!(wkNLEQ2, "P_SUMXALL", sumXall)
    setOption!(wkNLEQ2, "P_DLEVFALL", dLevFall)
    setOption!(wkNLEQ2, "P_SUMXQALL", sumXQall)
    setOption!(wkNLEQ2, "P_TOLALL", tolAll)
    setOption!(wkNLEQ2, "P_FCALL", fcAll)

    return (x, stats, retCode);
end

@debug function n2int(n, fcn, x, xScal, rTol, nItmax, nonLin, iRank, cond, opt, retCode,
    m1, m2, nBroy, xIter, sumXall, dLevFall, sumXQall, tolAll, fcAll, fc, fcMin,
    sigma, sigma2, mPrWarn, mPrMon, mPrSol, printIOwarn, printIOmon,
    printIOsol, qBDamp)
    # --------------------------------------------------------------------------
    # Since wkNLEQ2 is module global
    # Create the local variables here rather than taking them as arguments
    a       = wkNLEQ2.options[WK_A]
    qa      = wkNLEQ2.options[WK_QA_DXSAVE]
    dxSave  = wkNLEQ2.options[WK_QA_DXSAVE]
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
    iCnv = 0
    fck2 = fc
    # --------------------------------------------------------------------------
    # 0.2 Persistent variables
    cLin0   = getOption!(wkNLEQ2,"P_CLIN0",0.0)
    cLin1   = getOption!(wkNLEQ2,"P_CLIN1",0.0)
    cAlpha  = getOption!(wkNLEQ2,"P_CALPHA",0.0)
    alphaE  = getOption!(wkNLEQ2,"P_ALPHAE",0.0)
    alphaK  = getOption!(wkNLEQ2,"P_ALPHAK",0.0)
    alphaA  = getOption!(wkNLEQ2,"P_ALPHAA",0.0)
    qMStop  = getOption!(wkNLEQ2,"P_QMSTOP",false)
    sumxa2  = getOption!(wkNLEQ2,"P_SUMXA2",0.0)
    d       = getOption!(wkNLEQ2,"P_D",zeros(n))
    p       = getOption!(wkNLEQ2,"P_P",zeros(Int64,n))
    # --------------------------------------------------------------------------
    # 0.3 Machine related constants
    epMach  = 1e-17
    small   = 1e-150
    # --------------------------------------------------------------------------
    # Begin
    # --------------------------------------------------------------------------
    # 1 Initialization
    # --------------------------------------------------------------------------
    n           = length(x)
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
    qMixIO &= mPrMon != 0 && mPrSol != 0
    # --------------------------------------------------------------------------
    # 1.2 Derived dimensional parameters
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
    qIter               = true
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
        wkNLEQ2.options[STATS_NITER] = nIter

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
            # iFail = 0
        catch
            iFail = -1
            retCode = 82
            qIter = false
        end
        nFcn += 1
        if length(f) != n
            retCode = 22
            qIter   = false
            #=throw(DimensionMismatch("Dimension of the function output does not
            match the input dimension. Please check the function \"$fcn\" again."))=#
        end
    else
        qIniSc = false
    end

    # --------------------------------------------------------------------------
    # Main iteration loop
    # Repeat
    @bp
    while qIter
        # ----------------------------------------------------------------------
        # 2 Startup of iteration step
        if !qJcRfr
            # ------------------------------------------------------------------
            # 2.1 Scaling of variables x(n)
            xw = nScal(n , x, xa, xScal, iScal, qIniSc, opt)
            qIniSc = false
            if nIter != 0
                # Preliminary psuedo-rank
                iRankA = iRank
                # --------------------------------------------------------------
                # 2.2 Aposteriori estimate of damping factor
                dxQa[:] = dxQ
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
        if qGenJ
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
                if jacGen == 3
                    (a,eta,nFcnJ,iFail) = nJcf(fcn,n,n,x,f,xw,eta,etaMin,etaMax,etaDif,conv,nFcnJ)
                end
                if jacGen == 2
                    (a,nFcnJ,iFail) = nJac(fcn,n,n,x,f,xw,aJdel,aJmin,nFcnJ)
                end
                if jacGen == 4
                    (a,iFail) = n2jacFAD(fcn,x)
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
            for k = 1:n
                a[1:n,k] = -a[1:n,k]*xw[k]
            end
            # ------------------------------------------------------------------
            # 2.4.2 Row scaling of matrix A
            if qScale
                (a,fw) = nScrf(n,n,a)
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
        qPseudoRed = true
        # ----------------------------------------------------------------------
        # 3 Central part of iteration step
        # Pseudo-rank reduction loop
        # ==========================
        while qPseudoRed
            # ------------------------------------------------------------------
            # 3.1 Solution of the linear system
            # ------------------------------------------------------------------
            # 3.1.1 Decomposition of (n,n) matrix A
            if newt == 0
                cond1 = cond
                # TODO: Check if this is required or not
                if qRepeat
                    iRepeat = 1
                else
                    iRepeat = 0
                end
                # TODO: remember to set WK_SENS1 inside n2fact
                @bp
                (cond1,iFail) = n2fact(n,m1,n,1,1,a,qa,cond1,iRank,opt,p,d,iRepeat)
                if iFail != 0
                    retCode = 80
                    break
                end
                # Get the sensitivity of the Jacobian as estimated by n2fact
                sens1 = wkNLEQ2.options[WK_SENS1]
            end
            qLInit = true
            # ------------------------------------------------------------------
            # 3.1.2 Solution of linear (n,n) system
            if newt == 0
                iFail = n2solv(n,m1,n,1,1,a,qa,t1,t2,iRank,opt,iRepeat,d,p)
                if iFail != 0
                    retCode = 81
                    break
                end
                if !qRepeat && iRank != 0
                    qu[:] = t1
                end
            else
                alfa1 = sum((dx.*dxQ)./(xw.^2))
                alfa2 = sum((dx.^2)./(xw.^2))
                alfa  = alfa1/alfa2
                beta  = 1.0-alfa
                t2 = (dxQ+(fcA-1.0)*alfa*dx)/beta
                if newt == 1
                    dxSave[1:n,1] = dx
                end
                dxSave[1:n,newt+1] = t1
                dx = t1
                t1 = t1./xw
            end
            # ------------------------------------------------------------------
            # 3.2 Evaluation of scaled natural level function sumX
            # scaled maximum error norm conv
            # evaluation of (scaled) standard level function dlevf
            # and computation of ordinary Newton corrections dx[n]
            (dx,conv,sumX,dLevF) = nLvls(n,dx,t2,xw,f,mPrMon,newt == 0)
            wkNLEQ2.options[STATS_SUMX]   = sumX
            wkNLEQ2.options[STATS_DLEVF]  = dLevF
            xa[:]    = x
            sumXa    = sumX
            dLevXa   = sqrt(sumXa/n)
            conva    = conv
            dxANrm   = wnorm(n,dx,xw)
            push!(sumXall,dLevXa)
            push!(dLevFall,dLevF)
            # ------------------------------------------------------------------
            # 3.3 A-priori estimate of damping factor FC
            qRedu = false
            if nIter != 0 && nonLin != 1
                if newt == 0 || qRepeat
                    # ----------------------------------------------------------
                    # 3.3.1 Computation of the denominator of A-priori estimate
                    fcDnm = sum(((dx-dxQ)./xw).^2)
                    if iRank != n
                        # ------------------------------------------------------
                        # 3.3.2 Rank-deficient case (if previous rank was full)
                        # computation of the projected denominator of a-priori
                        # estimate
                        t1 = dxQ./xw
                        # Norm of projection of reduced component t1[n]
                        # TODO: Get "p" and "t2" from n2fact and n2solv
                        # "p" is the pivot vector
                        # "t2" is the diagonal returned by deccon
                        del = n2prjn(n, iRank, t1, d, qa, p)
                        fcDnm -= del
                    end
                    fcDnm *= sumX
                    # ----------------------------------------------------------
                    # 3.3.3 New damping factor
                    if fcDnm > fcNumP*fcMin2 || (nonLin == 4 && fcA*fcA*fcNumP < 4.0*fcDnm)
                        dMyPri = fcA*sqrt(fcNumP/fcDnm)
                        fcPri  = min(dMyPri,1.0)
                        if nonLin == 4
                            fcPri = min(0.5*dMyPri,1.0)
                        end
                    else
                        fcPri = 1.0
                        dMyPri = -1.0
                    end

                    if mPrMon >= 5
                        write(printIOmon,"\n",
                        @sprintf("+++++  apriori estimate  +++++\n"),
                        @sprintf(" fcPri  = %18.10e\n",fcPri),
                        @sprintf(" fc     = %18.10e\n",fc),
                        @sprintf(" fcA    = %18.10e\n",fcA),
                        @sprintf(" dMyPri = %18.10e\n",dMyPri),
                        @sprintf(" fcNumP = %18.10e\n",fcNumP),
                        @sprintf(" fcDnm  = %18.10e\n",fcDnm),
                        @sprintf("++++++++++++++++++++++++++++++\n"))
                    end

                    fc = fcPri
                    if iRank == n || iRank <= minRnk
                        fc = max(fc,fcMin)
                    end
                    if qBDamp
                        fcbh = fcA*fcBand
                        if fc > fcbh
                            fc = fcbh
                            if mPrMon >= 4
                                write(printIOmon, "*** Increase rest. act. (a priori)\n")
                            end
                        end
                        fcbh = fcA/fcBand
                        if fc < fcbh
                            fc = fcbh
                            if mPrMon >= 4
                                write(printIOmon, "*** Decrease rest. act. (a priori)\n")
                            end
                        end
                    end
                end
                qRedu = fc < fcMin
            end
            qRepeat = false
            if iOrMon >= 2
                sumxa2 = sumxa1
                sumxa1 = sumxa0
                sumxa0 = dLevXa
                if sumxa0 == 0.0
                    sumxa0 = small
                end
                # Check convergence rate (linear and superlinear)
                # iconv : Convergence idicator
                #           = 0: No convergence indicated yet
                #           = 1: Damping factor is 1.0e0
                #           = 2: Superlinear convergence detected (alpha >=1.2)
                #           = 3: Quadratic convergence detected (alpha > 1.8)
                fcMon = min(fc,fcMon)
                if fcMon < 1.0
                    iConv  = 0
                    alphaE = 0.0
                end
                if fcMon == 1.0 && iConv == 0
                    iConv = 1
                end
                if nIter >= 1
                    cLin1 = cLin0
                    cLin0 = sumxa0/sumxa1
                end
                if iConv >= 1 && nIter >= 2
                    alphaK = alphaE
                    alphaE = 0.0
                    if cLin1 <= 0.95
                        alphaE = log(cLin0)/log(cLin1)
                    end
                    if alphaK != 0.0
                        alphaK = 0.5*(alphaE+alphaK)
                    end
                    alphaA = min(alphaK,alphaE)
                    cAlphaK = cAlpha
                    cAlpha = 0.0
                    if alphaE != 0.0
                        cAlpha = sumxa1/sumxa2^alphaE
                    end
                    sumXte = sqrt(cAlpha*cAlphaK)*sumxa1^alphaK-sumxa0
                    if alphaA >= 1.2 && iConv == 1
                        iConv = 2
                    end
                    if alphaA > 1.8
                        iConv = 3
                    end
                    if mPrMon >= 4
                        write(printIOmon,"\n",
                        @sprintf(" ** iConv: %1i",iConv),
                        @sprintf("  alpha:       %9.2e",alphaE),
                        @sprintf("  const-alpha: %9.2e",cAlpha),
                        @sprintf("  const-lin:   %9.2e **\n",cLin0),
                        @sprintf(" ** alpha-post: %9.2e",alphaK),
                        @sprintf("  check:       %9.2e",sumXte),
                        @sprintf("                    **\n"))
                    end
                    if iConv >= 2 && alphaA < 0.9
                        if iOrMon == 3
                            retCode = 4
                            break
                        else
                            qMStop = true
                        end
                    end
                end
            end
            fcMon = fc

            if mPrMon >= 2
                n2prv1(dLevF, dLevXa, fcKeep, nIter, newt, iRank, mPrMon,
                    printIOmon, qMixIO, cond1, iRank)
            end

            if !qRedu
                # --------------------------------------------------------------
                # 3.4 Save natural level for later computations of corrector
                # and print iterate
                fcNumK = sumX
                nRed = 0
                qRep = false
                qDampRed = true
                # Damping-factor reduction loop
                # =============================
                while qDampRed
                    # ----------------------------------------------------------
                    # 3.5 Preliminary new iterate
                    x = xa + dx*fc
                    push!(fcAll,fc)
                    # ----------------------------------------------------------
                    # 3.5.2 Exit, if problem is specified as being linear
                    if nonLin == 1
                        retCode = 0
                        break
                    end
                    #-----------------------------------------------------------
                    # 3.6.1 Computation of the residual vector
                    try
                        fcn(f,x)
                        iFail = 0
                    catch
                        iFail = -1
                        retCode = 82
                    end
                    nFcn += 1
                    if iFail < 0
                        retCode = 82
                        break
                    end
                    if iFail == 1 || iFail == 2
                        if iFail == 1
                            fcRedu = 0.5
                        else
                            fcRedu = f[1]
                            if fcRedu <= 0.0 || fcRedu >= 1.0
                                retCode = 83
                                break
                            end
                        end
                        if mPrMon >= 2
                            write(printIOmon,
                            @sprintf("        %2i",nIter),
                            @sprintf(" %s could not be evaluated     %7.5f    %2i\n",fcn,fc,newt))
                        end
                        fch = fc
                        fc  = fcRedu*fc
                        if fch > fcMin
                            fc = max(fc,fcMin)
                        end
                        if qBDamp
                            fcbh = fch/fcBand
                            if fc < fcbh
                                fc = fcbh
                                if mPrMon >= 4
                                    write(printIOmon," *** Decrease rest. act. (fcn redu.) ***\n")
                                end
                            end
                        end
                        if fc < fcMin
                            retCode = 3
                            break
                        end
                    else
                        t1 = f.*fw
                        # ------------------------------------------------------
                        # 3.6.2 Solution of linear (n,n) system
                        if qRepeat
                            iRepeat = 1
                        else
                            iRepeat = 0
                        end
                        iFail = n2solv(n,m1,n,1,1,a,qa,t1,t2,iRank,opt,iRepeat,d,p)
                        if iFail != 0
                            retCode = 81
                            break
                        end
                        if newt > 0
                            dxQ = t2.*xwa
                            for iLoop = 1:newt
                                sum1 = sum((dxQ.*dxSave[1:n,iLoop])./xw.^2)
                                sum2 = sum((dxSave[1:n,iLoop]./xw).^2)
                                beta = sum1/sum2
                                dxQ = dxQ + beta*dxSave[1:n,iLoop+1]
                                t2 = dxQ./xw
                            end
                        end
                        # ------------------------------------------------------
                        # 3.6.3 Evaluation of scaled natural level function
                        #       sumX
                        #       scaled maximum error norm conv and evaluation
                        #       of (scaled) standard level function dLevFn
                        (dxQ,conv,sumX,dLevFn) =
                            nLvls(n,dxQ,t1,xw,f,mPrMon,newt==0)

                        push!(sumXQall,sqrt(sumX/n))
                        dxNrm = wnorm(n,dxQ,xw)
                        # ------------------------------------------------------
                        # 3.6.4 Convergence test
                        push!(tolAll,dxNrm)
                        if dxNrm <= rTol && dxANrm <= rSmall && fc == 1.0
                            retCode = 0
                            iCnv = 1
                            break
                        end

                        fcA = fc
                        # ------------------------------------------------------
                        # 3.6.5 Evaluation of reduced damping factor
                        th = fcA - 1.0
                        fcDnm = sum(((dxQ+th*dx)./xw).^2)
                        if fcDnm != 0.0
                            dMyCor = fcA*fcA*0.5*sqrt(fcNumK/fcDnm)
                        else
                            dMyCor = 1.0e+35
                        end
                        if nonLin <= 3
                            fcCor = min(1.0,dMyCor)
                        else
                            fcCor = min(1.0,0.5*dMyCor)
                        end

                        if mPrMon >= 5
                            write(printIOmon,
                            @sprintf(" +++ corrector computation +++\n"),
                            @sprintf("  fcCor    = %18.10e\n",fcCor),
                            @sprintf("  fc       = %18.10e\n",fc),
                            @sprintf("  dMyCor   = %18.10e\n",dMyCor),
                            @sprintf("  fcNumK   = %18.10e\n",fcNumK),
                            @sprintf("  fcDnm    = %18.10e\n",fcDnm),
                            @sprintf("  fcA      = %18.10e\n",fcA),
                            @sprintf("++++++++++++++++++++++++++++++\n"))
                        end
                        # ------------------------------------------------------
                        # 3.7 Natural monotonicity test
                        if sumX > sumXa
                            # --------------------------------------------------
                            # 3.8 Output of iterate
                            if mPrMon >= 3
                                nPrv2(dLevFn,sqrt(sumX/n),fc,nIter,mPrMon,
                                    printIOmon,qMixIO,"*")
                            end
                            if qMStop
                                retCode = 4
                                break
                            end
                            fc = min(fcCor,0.5*fc)
                            if (iRank == n || iRank == minRnk) && fcA > fcMin
                                fc = max(fc,fcMin)
                            end
                            if qBDamp
                                fcbh = fcA/fcBand
                                if fc < fcbh
                                    fc = fcbh
                                    if mPrMon >= 4
                                        write(printIOmon,
                                        " *** Decrease rest. act. (a posteriori) ***\n")
                                    end
                                end
                            end
                            fcMon = fc

                            if mPrMon >= 5
                                write(printIOmon,
                                " +++ corrector setting 1 +++\n",
                                @sprintf("fc    = %18.10e\n",fc),
                                " +++++++++++++++++++++++++++\n")
                            end

                            qRep = true
                            nCorr += 1
                            nRed += 1
                            # --------------------------------------------------
                            # 3.10 If damping factor is too small:
                            #      Refreash Jacobian, if current Jacobian
                            #      was computed by a Rank-1 update, else
                            #      reduce psuedo-rank
                            qRedu = fc < fcMin || newt > 0 && nRed > 1
                        else
                            if !qRep && fcCor > sigma2*fc
                                if mPrMon >= 3
                                    nPrv2(dLevFn,sqrt(sumX/n),fc,nIter,
                                    mPrMon,printIOmon,qMixIO,"+")
                                end
                                fc = fcCor

                                if mPrMon >= 5
                                    write(printIOmon,
                                    " +++ corrector setting 2 +++\n",
                                    @sprintf("fc    = %18.10e\n",fc),
                                    " +++++++++++++++++++++++++++\n")
                                end

                                qRep = true
                            else
                                qNext = true
                            end
                        end
                    end
                    qDampRed = !(qNext||qRedu)
                end
                # end of damping-factor reduction loop
                # ====================================
            end

            if nonLin == 1 || iCnv == 1 || (retCode != 0 && retCode != -1)
                break
            end

            if qRedu
                # --------------------------------------------------------------
                # 3.11 Restore former values for repeating step
                nRejR1 += 1
                x[:]   = xa
                f[:]   = fa
                dxQ[:] = dxQa
                if mPrMon >= 2
                    write(printIOmon,
                    @sprintf("        %2i not accepted damping factor %7.5f",nIter,fc),
                    @sprintf("    %2i %4i\n",newt,iRank))
                end
                fc  = fcKeep
                fcA = fck2
                if nIter == 0
                    fc = fcMin
                end
                if newt > 0
                    qGenJ  = true
                    qJcRfr = true
                    qRedu  = false
                else
                    # ----------------------------------------------------------
                    # 3.12 Psuedo-rank reduction
                    qRepeat = true
                    t1[:] = qu
                    iRank -= 1
                    if iRank < minRnk
                        retCode = 3
                        break
                    end
                end
            end
            qPseudoRed = !(!qRedu)
        end
        # end of psuedo-rank reduction loop
        # =================================
        if qNext
            # ------------------------------------------------------------------
            # 4 Preparations to start the following iteration step
            # ------------------------------------------------------------------
            # Print values
            if mPrMon >= 3
                nPrv2(dLevFn,sqrt(sumX/n),fc,nIter+1,mPrMon,printIOmon,qMixIO,"*")
            end
            # Print the natural level of the current iterate and
            # return it in one-step mode
            sumX = sumXa
            if mPrSol >= 2 && nIter != 0
                nSout(n,xa,2,opt,wkNLEQ2,mPrSol,printIOsol)
            elseif mPrSol >= 1 && nIter == 0
                nSout(n,xa,1,opt,wkNLEQ2,mPrSol,printIOsol)
            end
            nIter += 1
            wkNLEQ2.options[STATS_NITER] = nIter
            push!(xIter,x)
            dLevF = dLevFn
            if nIter >= nItmax
                retCode = 2
                break
            end
            fcKeep = fc
            # ------------------------------------------------------------------
            # 4.2 Return if in one-step mode
            if mode == 1
                qSucc = true
                setOption!(opt, OPT_QSUCC, Int(qSucc))
                setOption!(opt, OPT_FCSTART, fc)

                setOption!(wkNLEQ2, STATS_NITER,  nIter)
                setOption!(wkNLEQ2, STATS_NCORR,  nCorr)
                setOption!(wkNLEQ2, STATS_NFCN,   nFcn)
                setOption!(wkNLEQ2, STATS_NFCNJ,  nFcnJ)
                setOption!(wkNLEQ2, STATS_NJAC,   nJac)
                setOption!(wkNLEQ2, STATS_NREJR1, nRejR1)
                setOption!(wkNLEQ2, STATS_NEW,    newt)
                setOption!(wkNLEQ2, STATS_ICONV,  iConv)

                setOption!(wkNLEQ2, WK_A, a)
                setOption!(wkNLEQ2, WK_DXSAVE, dxSave)
                setOption!(wkNLEQ2, WK_DX, dx)
                setOption!(wkNLEQ2, WK_DXQ, dxQ)
                setOption!(wkNLEQ2, WK_DXQA, dxQa)
                setOption!(wkNLEQ2, WK_XA, xa)
                setOption!(wkNLEQ2, WK_XW, xw)
                setOption!(wkNLEQ2, WK_XWA, xwa)
                setOption!(wkNLEQ2, WK_F, f)
                setOption!(wkNLEQ2, WK_FA, fa)
                setOption!(wkNLEQ2, WK_FW, fw)
                setOption!(wkNLEQ2, WK_ETA, eta)
                setOption!(wkNLEQ2, WK_SUMXA0, sumxa0)
                setOption!(wkNLEQ2, WK_SUMXA1, sumxa1)
                setOption!(wkNLEQ2, WK_FCMON, fcMon)
                setOption!(wkNLEQ2, WK_FCA, fcA)
                setOption!(wkNLEQ2, WK_FCKEEP, fcKeep)
                setOption!(wkNLEQ2, WK_FCPRI, fcPri)
                setOption!(wkNLEQ2, WK_DMYCOR, dMyCor)
                setOption!(wkNLEQ2, STATS_CONV, conv)
                setOption!(wkNLEQ2, STATS_SUMX, sumX)
                setOption!(wkNLEQ2, WK_SUMXS, sumXs)
                setOption!(wkNLEQ2, STATS_DLEVF, dLevF)

                setOption!(wkNLEQ1, "P_CLIN0", cLin0)
                setOption!(wkNLEQ1, "P_CLIN1", cLin1)
                setOption!(wkNLEQ1, "P_CALPHA", cAlpha)
                setOption!(wkNLEQ1, "P_ALPHAE", alphaE)
                setOption!(wkNLEQ1, "P_ALPHAK", alphaK)
                setOption!(wkNLEQ1, "P_ALPHAA", alphaA)
                setOption!(wkNLEQ1, "P_QMSTOP", qMStop)
                setOption!(wkNLEQ1, "P_SUMXA2", sumxa2)
                setOption!(wkNLEQ1, "P_D", d)
                setOption!(wkNLEQ1, "P_P", p)

                return (x, xScal, retCode, wkNLEQ2)
            end
        end
    end
    # end of main iteration loop
    # ==========================
    # --------------------------------------------------------------------------
    # 9 Exits
    # --------------------------------------------------------------------------
    # 9.1 Solution exit
    aprec = -1.0

    if retCode == 0 || retCode == 4
        if nonLin != 1
            if retCode == 0
                aprec = sqrt(sumX/n)
                x += dxQ
                push!(xIter,x)
            else
                aprec = sqrt(sumXa/n)
                if alphaA > 0.0 && iOrMon == 3
                    x += dx
                    push!(xIter,x)
                end
            end
            if iRank < n
                retCode = 1
            end
            # Print final monitor output
            if mPrMon >= 2
                if retCode == 0
                    nPrv2(dLevFn,sqrt(sumX/n),fc,nIter+1,
                    mPrMon,printIOmon,qMixIO,"*")
                elseif iOrMon == 3
                    n2prv1(dLevFn,sqrt(sumXa/n),fc,nIter,newt,
                    mPrMon,printIOmon,qMixIO,cond1,iRank)
                end
            end
            if iOrMon >= 2
                if iConv <= 1 && alphaE != 0.0 && alphaK != 0.0
                    retCode = 5
                end
            end
            if mPrMon >= 1
                if retCode == 4
                    nOut = nIter
                else
                    nOut = nIter + 1
                end
                write(printIOmon,"\n\n\n ",
                @sprintf("Solution of nonlinear system of equations obtained within "),
                @sprintf("%3i iteration steps\n\n",nOut),
                @sprintf("Achieved relative accuracy %10.3e\n",aprec))
            end
        else
            if mPrMon >= 1
                write(printIOmon,"\n\n\n ",
                @sprintf("Solution of nonlinear system of equations obtained by NLEQ2\n"),
                @sprintf("No estimate available for the achieved relative accuracy\n"))
            end
        end
    end
    # --------------------------------------------------------------------------
    # 9.2 Fail exit messages
    # --------------------------------------------------------------------------
    # 9.2.1 Termination at stationary point
    if retCode == 1 && mPrWarn == 1
        write(printIOwarn,"\nIteration terminated at stationary point\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.2 Termination after more than nItmax iterations
    if retCode == 2 && mPrWarn == 1
        write(printIOwarn,"\n",
        @sprintf("Iteration terminates after nItmax %3i iteration steps\n",nItmax))
    end
    # --------------------------------------------------------------------------
    # 9.2.3 Newton method fails to converge
    if retCode == 3 && mPrWarn == 1
        write(printIOwarn,"\n",
        @sprintf("Newton method fails to converge\n"))
    end
    # --------------------------------------------------------------------------
    # 9.2.4.1 Superlinear convergence slowed down
    if retCode == 4 && mPrWarn == 1
        if iConv == 2
            ctyp = "superlinear"
        end
        if iConv == 3
            ctyp = "quadratic"
        end
        # TODO: write a wrapper for error(), warn(), info() such that it prints
        # either on STDOUT or in file
        if qMStop
            write(printIOwarn,"\nWARNING: Monotonicity test failed after ",ctyp,
            " convergence was already checked\nrTol requirement may be too",
            " stringent\n")
        else
            write(printIOwarn,"\nWARNING: ",ctyp, " convergence slowed down\n",
            "rTol requirement may be too stringent\n")
        end
    end
    # --------------------------------------------------------------------------
    # 9.2.4.2 Convergence criterion satisfied before superlinear
    #         convergence has been established
    if retCode == 5 && dLevFn == 0.0
        retCode = 0
    end
    if retCode == 5 && mPrWarn == 1
        write(printIOwarn,"\n",
        "WARNING: No quadratic or superlinear convergence established yet\n",
        "         your solution may perhaps be less accurate\n",
        "         as indicated by the standard error estimate\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.5 Error exit due to linear solver routine n2fact
    if retCode == 80 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by linear solver n2fact\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.6 Error exit due to linear solver routine n2solv
    if retCode == 81 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by linear solver n2solv\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.7 Error exit due to fail of user function fcn
    if retCode == 82 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by user function fcn\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.8 Error exit due to fail of user function Jacobian
    if retCode == 83 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by user function jac\n")
    end
    if (retCode == 82 || retCode == 83) && nIter <= 1 && mPrWarn == 1
        write(printIOwarn,"Try to find a better initial guess for the solution\n")
    end
    # --------------------------------------------------------------------------
    # 9.3 Common exit
    if mPrWarn == 1 && retCode != 0 && retCode != 4 && nonLin != 1
        write(printIOwarn,"\n",@sprintf("Achieved relative accuracy %10.3e\n",conva))
        aprec = conva
    end
    if mPrMon >= 1
        write(printIOmon,"\n",@sprintf("Subcondition ( 1, %4i ) %10.3e",iRank,cond1))
        write(printIOmon,"\n",@sprintf("Sensitivity  ( 1, %4i ) %10.3e",iRank,sens1))
    end
    rTol = aprec
    sumX = sumXa
    if mPrSol >= 2 && nIter != 0
        nSout(n,xa,2,opt,wkNLEQ2,mPrSol,printIOsol)
    elseif mPrSol >= 1 && nIter == 0
        nSout(n,xa,1,opt,wkNLEQ2,mPrSol,printIOsol)
    end
    if retCode != 4
        nIter += 1
    end
    dLevF = dLevFn
    # Print solution or final iteration vector
    if mPrSol >= 1
        if retCode == 0
            modefi = 3
        else
            modefi = 4
        end
        nSout(n,x,modefi,opt,wkNLEQ2,mPrSol,printIOsol)
    end
    # End of exits

    # --------------------------------------------------------------------------
    # 10 Prepare all the variables for returning
    xScal[:] = xw
    # TODO: Convert the setOptions to a function since it is called twice
    setOption!(opt, OPT_QSUCC, Int(qSucc))
    setOption!(opt, OPT_FCSTART, fc)

    setOption!(wkNLEQ2, STATS_NITER,  nIter)
    setOption!(wkNLEQ2, STATS_NCORR,  nCorr)
    setOption!(wkNLEQ2, STATS_NFCN,   nFcn)
    setOption!(wkNLEQ2, STATS_NFCNJ,  nFcnJ)
    setOption!(wkNLEQ2, STATS_NJAC,   nJac)
    setOption!(wkNLEQ2, STATS_NREJR1, nRejR1)
    setOption!(wkNLEQ2, STATS_NEW,    newt)
    setOption!(wkNLEQ2, STATS_ICONV,  iConv)

    setOption!(wkNLEQ2, WK_A, a)
    setOption!(wkNLEQ2, WK_DXSAVE, dxSave)
    setOption!(wkNLEQ2, WK_DX, dx)
    setOption!(wkNLEQ2, WK_DXQ, dxQ)
    setOption!(wkNLEQ2, WK_DXQA, dxQa)
    setOption!(wkNLEQ2, WK_XA, xa)
    setOption!(wkNLEQ2, WK_XW, xw)
    setOption!(wkNLEQ2, WK_XWA, xwa)
    setOption!(wkNLEQ2, WK_F, f)
    setOption!(wkNLEQ2, WK_FA, fa)
    setOption!(wkNLEQ2, WK_FW, fw)
    setOption!(wkNLEQ2, WK_ETA, eta)
    setOption!(wkNLEQ2, WK_SUMXA0, sumxa0)
    setOption!(wkNLEQ2, WK_SUMXA1, sumxa1)
    setOption!(wkNLEQ2, WK_FCMON, fcMon)
    setOption!(wkNLEQ2, WK_FCA, fcA)
    setOption!(wkNLEQ2, WK_FCKEEP, fcKeep)
    setOption!(wkNLEQ2, WK_FCPRI, fcPri)
    setOption!(wkNLEQ2, WK_DMYCOR, dMyCor)
    setOption!(wkNLEQ2, STATS_CONV, conv)
    setOption!(wkNLEQ2, STATS_SUMX, sumX)
    setOption!(wkNLEQ2, WK_SUMXS, sumXs)
    setOption!(wkNLEQ2, STATS_DLEVF, dLevF)

    setOption!(wkNLEQ1, "P_CLIN0", cLin0)
    setOption!(wkNLEQ1, "P_CLIN1", cLin1)
    setOption!(wkNLEQ1, "P_CALPHA", cAlpha)
    setOption!(wkNLEQ1, "P_ALPHAE", alphaE)
    setOption!(wkNLEQ1, "P_ALPHAK", alphaK)
    setOption!(wkNLEQ1, "P_ALPHAA", alphaA)
    setOption!(wkNLEQ1, "P_QMSTOP", qMStop)
    setOption!(wkNLEQ1, "P_SUMXA2", sumxa2)
    setOption!(wkNLEQ1, "P_D", d)
    setOption!(wkNLEQ1, "P_P", p)

    return (x, xScal, retCode)
    # End of function n2int
end
