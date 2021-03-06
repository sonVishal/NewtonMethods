"""
function nleq2{T}(fcn, x::Vector{T}, xScal::Vector{T}, opt::OptionsNLEQ)

Damped Newton-algorithm with rank strategy for systems of
highly nonlinear equations - damping strategy due to Ref.(1).

(The iteration is done by function N2INT currently. NLEQ2
itself does some house keeping and builds up workspace.)

Jacobian approximation by numerical differences, user
supplied function JAC or forward mode automatic differentation.

The numerical solution of the arising linear equations is
done by means of the subroutines DECCON and SOLCON (QR de-
composition with subcondition estimation, rank decision and
computation of the rank-deficient pseudoinverse) .
For special purposes these routines may be substituted.

This is a driver routine for the core solver N2INT.

Does not support Band mode for Jacobian storage.

T = Float64 or BigFloat

## Input parameters
| Variable   | Description                                                                              |
|:-----------|:-----------------------------------------------------------------------------------------|
| fcn        | Function for which zero is to be found. Should be in the form of fcn(y,x) with y = f(x). |
| x[1:n]     | Initial estimate of the solution.                                                        |
| xScal[1:n] | User scaling (lower threshold) of the iteration vector x                                 |
| opt        | Options for solving the nonlinear system. Valid options are listed below.                |

## Output parameters
| Variable | Description                                                                                   |
|:---------|:----------------------------------------------------------------------------------------------|
| x0[1:n]  | Solution values (or final values if exit before solution is reached).                         |
| stats    | A dictionary variable of additional output values. The fields are discussed below.            |
| retCode  | An integer value signifying the exit code. The meaning of the exit codes are discussed below. |
"""
function nleq2{T}(fcn, x::Vector{T}, xScal::Vector{T}, opt::OptionsNLEQ)
    # Begin
    # Check input parameters and options
    n = length(x)
    retCode = checkOptions(n, x, xScal, opt)

    # Create an empty statistics variable to be returned
    stats = Dict{AbstractString,Any}()

    # Exit if any parameter error was detected
    if retCode != 0
        println("Exit with return code $retCode")
        return (x, stats, retCode)
    end

#-------------------------------------------------------------------------------
# Printing related stuff
#-------------------------------------------------------------------------------
    # Print warning messages?
    printWarn   = opt.options[OPT_PRINTWARNING]
    # Print iteration summary?
    printMon    = getOption!(opt,OPT_PRINTITERATION,0)
    # Print solution summary?
    printSol    = getOption!(opt,OPT_PRINTSOLUTION,0)
#-------------------------------------------------------------------------------

    (m1, m2, nBroy, qBDamp) = initializeOptions(n, opt, 2, T)

    # Create a copy inside so that the original variable is untouched
    x0 = x[:]

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    retCode = -1

    # Call to n2int
    retCode = n2int(n, fcn, x0, xScal, opt.options[OPT_RTOL], opt.options[OPT_NITMAX],
        opt.options[OPT_NONLIN], opt.options[OPT_IRANK], opt.options[OPT_COND],
        opt, m1, m2, nBroy, opt.options[OPT_FCSTART], opt.options[OPT_FCMIN],
        opt.options[OPT_SIGMA], opt.options[OPT_SIGMA2], printWarn, printMon,
        printSol, opt.options[OPT_PRINTIOWARN], opt.options[OPT_PRINTIOMON],
        opt.options[OPT_PRINTIOSOL], qBDamp)

    # set stats variable
    stats[STATS_XSCAL] = xScal
    if retCode == -1
        stats[STATS_RTOL] = wkNLEQ2.options["P_TOLALL"][wkNLEQ2.options[STATS_NITER]]
    else
        stats[STATS_RTOL] = opt.options[OPT_RTOL]
    end
    if opt.options[OPT_STORE] == 1
        stats[STATS_XITER]      = wkNLEQ2.options["P_XITER"]
        stats[STATS_NATLEVEL]   = wkNLEQ2.options["P_SUMXALL"]
        stats[STATS_SIMLEVEL]   = wkNLEQ2.options["P_DLEVFALL"]
        stats[STATS_STDLEVEL]   = wkNLEQ2.options["P_SUMXQALL"]
        stats[STATS_DAMPINGFC]  = wkNLEQ2.options["P_FCALL"]
    end
    stats[STATS_PRECISION]  = wkNLEQ2.options["P_TOLALL"]
    stats[STATS_NITER]      = wkNLEQ2.options[STATS_NITER]
    stats[STATS_NCORR]      = wkNLEQ2.options[STATS_NCORR]
    stats[STATS_NREJR1]     = wkNLEQ2.options[STATS_NREJR1]
    stats[STATS_NJAC]       = wkNLEQ2.options[STATS_NJAC]
    stats[STATS_NFCN]       = wkNLEQ2.options[STATS_NFCN]
    stats[STATS_NFCNJ]      = wkNLEQ2.options[STATS_NFCNJ]
    stats[STATS_SUBCOND]    = wkNLEQ2.options[STATS_SUBCOND]
    stats[STATS_SENS]       = wkNLEQ2.options[STATS_SENS]
    stats[STATS_IFAIL]      = wkNLEQ2.options[STATS_IFAIL]

    # Print statistics
    if printMon >= 2 && retCode != -1 && retCode != 10
        printStats(stats, opt.options[OPT_PRINTIOMON])
    end

    return (x0, stats, retCode);
end

"""
function n2int{T}(n::Int64, fcn, x::Vector{T}, xScal::Vector{T},
    rTol::T, nItmax::Int64, nonLin::Int64, iRank::Int64, cond::T,
    opt::OptionsNLEQ, m1::Int64, m2::Int64, nBroy::Int64,
    fc::T, fcMin::T, sigma::T, sigma2::T, mPrWarn::Int64,
    mPrMon::Int64, mPrSol::Int64, printIOwarn, printIOmon, printIOsol, qBDamp::Bool)

Core routine for NLEQ2.

Damped Newton-algorithm with rank-strategy for systems of
highly nonlinear equations especially designed for
numerically sensitive problems.

T = Float64 or BigFloat

## Input parameters
| Variable     | Description                                                                              |
|:-------------|:-----------------------------------------------------------------------------------------|
| n            | Number of vector components                                                              |
| fcn          | Function for which zero is to be found. Should be in the form of fcn(y,x) with y = f(x). |
| x[1:n]*      | Initial estimate of the solution.                                                        |
| xScal[1:n]   | User scaling (lower threshold) of the iteration vector x                                 |
| rTol         | Required relative precision of solution components                                       |
| nItmax       | Maximum number of allowed iterations                                                     |
| nonLin       | Problem type specification. See OPT_NONLIN field in solver options                       |
| iRank        | Initially proposed (in) and final (out) rank of Jacobian                                 |
| cond         | Maximum permitted subcondition for rank-decision by linear solver.                       |
| opt          | Options for solving the nonlinear system. Valid options are listed below.                |
| m1           | Leading dimension of Jacobian array a                                                    |
| m2 = n       | For full mode                                                                            |
| = ml+mu+1    | For band mode                                                                            |
| nBroy        | Maximum possible consecutive iterative Broyden steps.                                    |
| fc           | Current Newton iteration damping factor.                                                 |
| fcMin        | Minimum permitted damping factor. fc < fcMin results in either of the following          |
|              | a. Recomputation of the Jacobian using difference approximation                          |
|              | b. Fail exit                                                                             |
| sigma        | Decision parameter for rank1-updates                                                     |
| sigma2       | Decision parameter for damping factor increasing to corrector                            |
| mStor        | Decision parameter for matrix storage. See option OPT_MSTOR in solver options.           |
| mPrWarn      | Decision parameter for printing warning messages                                         |
| mPrMon       | Decision parameter for printing iteration monitor                                        |
| mPrSol       | Decision parameter for printing solution                                                 |
| printIOwarn  | IO handle for printing warning                                                           |
| printIOmon   | IO handle for printing iteration monitor                                                 |
| printIOsol   | IO handle for printing solution                                                          |
| qBDamp       | Decision parameter for matrix storage. See option OPT_MSTOR in solver options.           |

(* marks inout parameters)

## Output parameters
| Variable | Description                                                                                   |
|:---------|:----------------------------------------------------------------------------------------------|
| x[1:n]*  | Solution values (or final values if exit before solution is reached).                         |
| retCode  | An integer value signifying the exit code. The meaning of the exit codes are discussed below. |
"""
function n2int{T}(n::Int64, fcn, x::Vector{T}, xScal::Vector{T},
    rTol::T, nItmax::Int64, nonLin::Int64, iRank::Int64, cond::T,
    opt::OptionsNLEQ, m1::Int64, m2::Int64, nBroy::Int64,
    fc::T, fcMin::T, sigma::T, sigma2::T, mPrWarn::Int64,
    mPrMon::Int64, mPrSol::Int64, printIOwarn, printIOmon, printIOsol, qBDamp::Bool)
    # Begin
    (epMach, small, _) = getMachineConstants(T)
    # --------------------------------------------------------------------------
    # Since wkNLEQ2 is module global
    # Create the local variables here rather than taking them as arguments
    if nBroy == 0
        qa       = wkNLEQ2.options[WK_A]
        dxSave   = wkNLEQ2.options[WK_A]
        nBroyNew = 1
    else
        qa      = wkNLEQ2.options[WK_QA_DXSAVE]
        dxSave  = wkNLEQ2.options[WK_QA_DXSAVE]
    end
    a        = wkNLEQ2.options[WK_A]
    dx       = wkNLEQ2.options[WK_DX]
    dxQ      = wkNLEQ2.options[WK_DXQ]
    xa       = wkNLEQ2.options[WK_XA]
    xwa      = wkNLEQ2.options[WK_XWA]
    f        = wkNLEQ2.options[WK_F]
    fa       = wkNLEQ2.options[WK_FA]
    eta      = wkNLEQ2.options[WK_ETA]
    xw       = wkNLEQ2.options[WK_XW]
    fw       = wkNLEQ2.options[WK_FW]
    dxQa     = wkNLEQ2.options[WK_DXQA]
    qu       = wkNLEQ2.options[WK_QU]
    sumxa0   = wkNLEQ2.options[WK_SUMXA0]
    sumxa1   = wkNLEQ2.options[WK_SUMXA1]
    fcMon    = wkNLEQ2.options[WK_FCMON]
    fcA      = wkNLEQ2.options[WK_FCA]
    fcKeep   = wkNLEQ2.options[WK_FCKEEP]
    fcPri    = wkNLEQ2.options[WK_FCPRI]
    dMyCor   = wkNLEQ2.options[WK_DMYCOR]
    conv     = wkNLEQ2.options[STATS_CONV]
    dLevF    = wkNLEQ2.options[STATS_DLEVF]
    nIter    = wkNLEQ2.options[STATS_NITER]
    nCorr    = wkNLEQ2.options[STATS_NCORR]
    nFcn     = wkNLEQ2.options[STATS_NFCN]
    nFcnJ    = wkNLEQ2.options[STATS_NFCNJ]
    nJac     = wkNLEQ2.options[STATS_NJAC]
    nRejR1   = wkNLEQ2.options[STATS_NREJR1]
    newt     = wkNLEQ2.options[STATS_NEW]
    iConv    = wkNLEQ2.options[STATS_ICONV]
    tolAll   = wkNLEQ2.options["P_TOLALL"]

    if opt.options[OPT_STORE] == 1
        xIter    = wkNLEQ2.options["P_XITER"]
        sumXall  = wkNLEQ2.options["P_SUMXALL"]
        dLevFall = wkNLEQ2.options["P_DLEVFALL"]
        sumXQall = wkNLEQ2.options["P_SUMXQALL"]
        fcAll    = wkNLEQ2.options["P_FCALL"]
    end
    # --------------------------------------------------------------------------
    # 0.1 Variables that need to be defined before since they appear in different
    # scopes. The declaration and usage are in different scopes.
    retCode = -1
    fck2   = fc
    dLevFn = zero(T)
    sumXa  = zero(T)
    conva  = zero(T)
    cond1  = zero(T)
    sens1  = zero(T)
    iRankC = 0
    sumX   = zero(T)
    t2     = zeros(T,n)
    # --------------------------------------------------------------------------
    # 0.2 Persistent variables
    cLin0   = getOption!(wkNLEQ2,"P_CLIN0",zero(T))
    cLin1   = getOption!(wkNLEQ2,"P_CLIN1",zero(T))
    cAlpha  = getOption!(wkNLEQ2,"P_CALPHA",zero(T))
    alphaE  = getOption!(wkNLEQ2,"P_ALPHAE",zero(T))
    alphaK  = getOption!(wkNLEQ2,"P_ALPHAK",zero(T))
    alphaA  = getOption!(wkNLEQ2,"P_ALPHAA",zero(T))
    qMStop  = getOption!(wkNLEQ2,"P_QMSTOP",false)
    sumxa2  = getOption!(wkNLEQ2,"P_SUMXA2",zero(T))
    d       = getOption!(wkNLEQ2,"P_D",zero(x))
    p       = getOption!(wkNLEQ2,"P_P",zeros(Int64,n))
    # --------------------------------------------------------------------------
    # Begin
    # --------------------------------------------------------------------------
    # 1 Initialization
    # --------------------------------------------------------------------------
    qBreak      = false
    qSucc       = Bool(opt.options[OPT_QSUCC])
    qScale      = opt.options[OPT_NOROWSCAL] != 1
    qRank1      = Bool(opt.options[OPT_QRANK1])
    iOrMon      = getOption!(opt, OPT_IORMON, 2)
    iScal       = getOption!(opt, OPT_ISCAL, 0)
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
        fc = one(T)
    end
    # --------------------------------------------------------------------------
    # 1.5 Initial preparations
    qJcRfr              = false
    qRepeat             = false
    qIter               = true
    iFail               = 0
    fcBand              = zero(T)
    if qBDamp
        fcBand = opt.options[OPT_FCBAND]
    end
    # --------------------------------------------------------------------------
    # 1.5.1 Numerical differentiation related initializations
    if jacGen == 2
        aJdel = getOption!(opt, OPT_AJDEL, zero(T))
        if aJdel <= small
            aJdel = sqrt(epMach*10.0)
        end
        aJmin = getOption!(opt, OPT_AJMIN, zero(T))
    elseif jacGen == 3
        etaDif = getOption!(opt, OPT_ETADIF, zero(T))
        if etaDif <= small
            etaDif = T(1.0e-6)
        end
        etaIni = getOption!(opt, OPT_ETAINI, zero(T))
        if etaIni <= small
            etaIni = T(1.0e-6)
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
        fcKeep  = fc
        fcA     = fc
        fcPri   = fc
        fcMon   = fc
        fcK2    = fc
        conv    = zero(T)
        if jacGen == 3
            eta[:] = etaIni*ones(T,n)
        end

        xa[:] = x

        iConv  = 0

        sumxa0 = zero(T)
        sumxa1 = zero(T)

        if opt.options[OPT_STORE] == 1
            push!(xIter,x)
        end
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
            flush(printIOmon)
        end
        # ----------------------------------------------------------------------
        # 1.7 Startup step
        # ----------------------------------------------------------------------
        # 1.7.1 Computation of residual vector
        try
            fcn(f,x)
        catch
            iFail = -1
        end
        nFcn += 1
        if length(f) != n || iFail != 0
            retCode = 82
            qIter = false
        end
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
            nScal(n, x, xa, xScal, iScal, mPrMon, printIOmon, xw)
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
                    !qRank1 || newt + 2 > nBroyNew

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
                        flush(printIOmon)
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
                catch
                    iFail = -1
                end
            else
                if jacGen == 2
                    (nFcnJ,iFail) = nJacFD(fcn,n,n,x,f,xw,aJdel,aJmin,nFcnJ,a)
                end
                if jacGen == 3
                    (nFcnJ,iFail) = nJcf(fcn,n,n,x,f,xw,eta,etaMin,etaMax,
                                            etaDif,conv,nFcnJ,a)
                end
                # Forward mode automatic differentiation using ForwardDiff
                if jacGen == 4
                    fd = zero(x)
                    try
                        fcn(fd,x)
                    catch
                        iFail = -1
                    end
                    nFcnJ += 1
                    if iFail == 0
                        try
                            a[:,:] = ForwardDiff.jacobian(fcn,fd,x)
                        catch
                            iFail = -1
                        end
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
            xwa[:] = xw
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
                nScrf(n,n,a,fw)
            else
                fw[:] = ones(n)
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
                if qRepeat
                    iRepeat = 1
                else
                    iRepeat = 0
                end
                (cond1,iRankC,iFail) = nFact(n,m1,n,1,1,a,qa,cond1,iRank,
                                        opt,p,d,iRepeat,iRankC)
                if iFail != 0
                    retCode = 80
                    qBreak = true
                    break
                end
                # Get the sensitivity of the Jacobian as estimated by nFact
                sens1 = wkNLEQ2.options[WK_SENS1]
            end
            # ------------------------------------------------------------------
            # 3.1.2 Solution of linear (n,n) system
            if newt == 0
                iFail = nSolv(n,m1,n,1,1,a,qa,t1,t2,iRank,iRepeat,d,p,iRankC)
                if iFail != 0
                    retCode = 81
                    qBreak = true
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
                dx[:] = t1
                t1[:] = t1./xw
            end
            # ------------------------------------------------------------------
            # 3.2 Evaluation of scaled natural level function sumX
            # scaled maximum error norm conv
            # evaluation of (scaled) standard level function dlevf
            # and computation of ordinary Newton corrections dx[n]
            (conv,sumX,dLevF) = nLvls(n,dx,t2,xw,f,newt == 0)
            wkNLEQ2.options[STATS_SUMX]   = sumX
            wkNLEQ2.options[STATS_DLEVF]  = dLevF
            xa[:]    = x
            sumXa    = sumX
            dLevXa   = sqrt(sumXa/n)
            conva    = conv
            dxANrm   = wnorm(n,dx,xw)
            if opt.options[OPT_STORE] == 1
                push!(sumXall,dLevXa)
                push!(dLevFall,dLevF)
            end
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
                        del = n2prjn(n, iRank, t1, d, qa, p, t2)
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
                        fcPri = one(T)
                        dMyPri = -one(T)
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
                        flush(printIOmon)
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
                                flush(printIOmon)
                            end
                        end
                        fcbh = fcA/fcBand
                        if fc < fcbh
                            fc = fcbh
                            if mPrMon >= 4
                                write(printIOmon, "*** Decrease rest. act. (a priori)\n")
                                flush(printIOmon)
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
                    alphaE = zero(T)
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
                    alphaE = zero(T)
                    if cLin1 <= 0.95
                        alphaE = log(cLin0)/log(cLin1)
                    end
                    if alphaK != 0.0
                        alphaK = 0.5*(alphaE+alphaK)
                    end
                    alphaA = min(alphaK,alphaE)
                    cAlphaK = cAlpha
                    cAlpha = zero(T)
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
                        flush(printIOmon)
                    end
                    if iConv >= 2 && alphaA < 0.9
                        if iOrMon == 3
                            retCode = 4
                            qBreak = true
                            break
                        else
                            qMStop = true
                        end
                    end
                end
            end
            fcMon = fc

            if mPrMon >= 2
                nPrv1(dLevF, dLevXa, fcKeep, nIter, newt, iRank, mPrMon,
                    printIOmon, qMixIO, cond1)
            end

            if !qRedu
                # --------------------------------------------------------------
                # 3.4 Save natural level for later computations of corrector
                # and print iterate
                fcNumK   = sumX
                nRed     = 0
                qRep     = false
                qDampRed = true
                # Damping-factor reduction loop
                # =============================
                while qDampRed
                    # ----------------------------------------------------------
                    # 3.5 Preliminary new iterate
                    x[:] = xa + dx*fc
                    if opt.options[OPT_STORE] == 1
                        push!(fcAll,fc)
                    end
                    # ----------------------------------------------------------
                    # 3.5.2 Exit, if problem is specified as being linear
                    if nonLin == 1
                        retCode = 0
                        qBreak = true
                        break
                    end
                    #-----------------------------------------------------------
                    # 3.6.1 Computation of the residual vector
                    try
                        fcn(f,x)
                    catch
                        iFail = -1
                    end
                    nFcn += 1
                    if iFail < 0
                        retCode = 82
                        qBreak = true
                        break
                    end
                    if iFail == 1 || iFail == 2
                        if iFail == 1
                            fcRedu = 0.5*one(T)
                        else
                            fcRedu = f[1]
                            if fcRedu <= 0.0 || fcRedu >= 1.0
                                retCode = 83
                                qBreak = true
                                break
                            end
                        end
                        if mPrMon >= 2
                            write(printIOmon,
                            @sprintf("        %2i",nIter),
                            @sprintf(" %s could not be evaluated     %7.5f    %2i\n",fcn,fc,newt))
                            flush(printIOmon)
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
                                    flush(printIOmon)
                                end
                            end
                        end
                        if fc < fcMin
                            retCode = 3
                            qBreak = true
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
                        iFail = nSolv(n,m1,n,1,1,a,qa,t1,t2,iRank,iRepeat,d,p,iRankC)
                        if iFail != 0
                            retCode = 81
                            qBreak = true
                            break
                        end
                        if newt > 0
                            dxQ[:] = t2.*xwa
                            for iLoop = 1:newt
                                sum1 = sum((dxQ.*dxSave[1:n,iLoop])./xw.^2)
                                sum2 = sum((dxSave[1:n,iLoop]./xw).^2)
                                beta = sum1/sum2
                                dxQ[:] = dxQ + beta*dxSave[1:n,iLoop+1]
                                t2 = dxQ./xw
                            end
                        end
                        # ------------------------------------------------------
                        # 3.6.3 Evaluation of scaled natural level function
                        #       sumX
                        #       scaled maximum error norm conv and evaluation
                        #       of (scaled) standard level function dLevFn
                        (conv,sumX,dLevFn) =
                            nLvls(n,dxQ,t2,xw,f,newt==0)

                        dxNrm = wnorm(n,dxQ,xw)
                        push!(tolAll,dxNrm)

                        if opt.options[OPT_STORE] == 1
                            push!(sumXQall,sqrt(sumX/n))
                        end
                        # ------------------------------------------------------
                        # 3.6.4 Convergence test
                        if dxNrm <= rTol && dxANrm <= rSmall && fc == 1.0
                            retCode = 0
                            qBreak = true
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
                            dMyCor = T(1.0e+35)
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
                            flush(printIOmon)
                        end
                        # ------------------------------------------------------
                        # 3.7 Natural monotonicity test
                        if sumX > sumXa
                            # --------------------------------------------------
                            # 3.8 Output of iterate
                            if mPrMon >= 3
                                nPrv2(dLevFn,sqrt(sumX/n),fc,nIter,
                                    printIOmon,qMixIO,"*")
                            end
                            if qMStop
                                retCode = 4
                                qBreak = true
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
                                        flush(printIOmon)
                                    end
                                end
                            end
                            fcMon = fc

                            if mPrMon >= 5
                                write(printIOmon,
                                " +++ corrector setting 1 +++\n",
                                @sprintf("fc    = %18.10e\n",fc),
                                " +++++++++++++++++++++++++++\n")
                                flush(printIOmon)
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
                                    printIOmon,qMixIO,"+")
                                end
                                fc = fcCor

                                if mPrMon >= 5
                                    write(printIOmon,
                                    " +++ corrector setting 2 +++\n",
                                    @sprintf("fc    = %18.10e\n",fc),
                                    " +++++++++++++++++++++++++++\n")
                                    flush(printIOmon)
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

            # Any break inside the damping-factor reduction loop comes here
            if qBreak
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
                    @sprintf("        %2i Not accepted damping factor         %7.5f",nIter,fc),
                    @sprintf("    %2i      %4i\n",newt,iRank))
                    flush(printIOmon)
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
                        qBreak = true
                        break
                    end
                end
            end
            qPseudoRed = !(!qRedu)
        end
        # end of psuedo-rank reduction loop
        # =================================

        # Any break inside the pseudo-rank reduction loop comes here
        if qBreak
            break
        end

        if qNext
            # ------------------------------------------------------------------
            # 4 Preparations to start the following iteration step
            # ------------------------------------------------------------------
            # Print values
            if mPrMon >= 3
                nPrv2(dLevFn,sqrt(sumX/n),fc,nIter+1,printIOmon,qMixIO,"*")
            end
            # Print the natural level of the current iterate and
            # return it in one-step mode
            sumX = sumXa
            if mPrSol >= 2 && nIter != 0
                nSout(n,xa,2,mPrSol,printIOsol,nIter,dLevF,sumX)
            elseif mPrSol >= 1 && nIter == 0
                nSout(n,xa,1,mPrSol,printIOsol,nIter,dLevF,sumX)
            end
            nIter += 1
            wkNLEQ2.options[STATS_NITER] = nIter
            if opt.options[OPT_STORE] == 1
                push!(xIter,x)
            end
            dLevF = dLevFn
            if nIter >= nItmax
                retCode = 2
                break
            end
            fcKeep = fc
            # ------------------------------------------------------------------
            # 4.2 Return if in one-step mode
            if opt.options[OPT_MODE] == 1
                qSucc = true
                setOption!(opt, OPT_QSUCC, Int(qSucc))
                setOption!(opt, OPT_FCSTART, fc)
                setOption!(opt, OPT_IRANK, iRank)

                setOption!(wkNLEQ2, STATS_NITER,  nIter)
                setOption!(wkNLEQ2, STATS_NCORR,  nCorr)
                setOption!(wkNLEQ2, STATS_NFCN,   nFcn)
                setOption!(wkNLEQ2, STATS_NFCNJ,  nFcnJ)
                setOption!(wkNLEQ2, STATS_NJAC,   nJac)
                setOption!(wkNLEQ2, STATS_NREJR1, nRejR1)
                setOption!(wkNLEQ2, STATS_NEW,    newt)
                setOption!(wkNLEQ2, STATS_ICONV,  iConv)
                setOption!(wkNLEQ2, STATS_SUBCOND, cond1)
                setOption!(wkNLEQ2, STATS_SENS, sens1)
                setOption!(wkNLEQ2, STATS_IFAIL,  iFail)

                setOption!(wkNLEQ2, WK_SUMXA0, sumxa0)
                setOption!(wkNLEQ2, WK_SUMXA1, sumxa1)
                setOption!(wkNLEQ2, WK_FCMON, fcMon)
                setOption!(wkNLEQ2, WK_FCA, fcA)
                setOption!(wkNLEQ2, WK_FCKEEP, fcKeep)
                setOption!(wkNLEQ2, WK_FCPRI, fcPri)
                setOption!(wkNLEQ2, WK_DMYCOR, dMyCor)
                setOption!(wkNLEQ2, STATS_CONV, conv)
                setOption!(wkNLEQ2, STATS_SUMX, sumX)
                setOption!(wkNLEQ2, STATS_DLEVF, dLevF)

                setOption!(wkNLEQ2, "P_CLIN0", cLin0)
                setOption!(wkNLEQ2, "P_CLIN1", cLin1)
                setOption!(wkNLEQ2, "P_CALPHA", cAlpha)
                setOption!(wkNLEQ2, "P_ALPHAE", alphaE)
                setOption!(wkNLEQ2, "P_ALPHAK", alphaK)
                setOption!(wkNLEQ2, "P_ALPHAA", alphaA)
                setOption!(wkNLEQ2, "P_QMSTOP", qMStop)
                setOption!(wkNLEQ2, "P_SUMXA2", sumxa2)
                return retCode
            end
        end
    end
    # end of main iteration loop
    # ==========================
    # --------------------------------------------------------------------------
    # 9 Exits
    # --------------------------------------------------------------------------
    # 9.1 Solution exit
    aprec = -one(T)

    if retCode == 0 || retCode == 4
        if nonLin != 1
            if retCode == 0
                aprec = sqrt(sumX/n)
                x[:] += dxQ
                if opt.options[OPT_STORE] == 1
                    push!(xIter,x)
                end
            else
                aprec = sqrt(sumXa/n)
                if alphaA > 0.0 && iOrMon == 3
                    x[:] += dx
                    if opt.options[OPT_STORE] == 1
                        push!(xIter,x)
                    end
                end
            end
            if iRank < n
                retCode = 1
            end
            # Print final monitor output
            if mPrMon >= 2
                if retCode == 0
                    nPrv2(dLevFn,sqrt(sumX/n),fc,nIter+1,
                    printIOmon,qMixIO,"*")
                elseif iOrMon == 3
                    nPrv1(dLevFn,sqrt(sumXa/n),fc,nIter,newt,
                    mPrMon,printIOmon,qMixIO,cond1,iRank)
                end
            end
            if iOrMon >= 2
                if iConv <= 1 && alphaE != 0.0 && alphaK != 0.0
                    retCode = 5
                end
            end
            if mPrMon >= 1 && retCode != 1
                if retCode == 4
                    nOut = nIter
                else
                    nOut = nIter + 1
                end
                write(printIOmon,"\n\n\n ",
                @sprintf("Solution of nonlinear system of equations obtained within "),
                @sprintf("%3i iteration steps\n\n",nOut),
                @sprintf("Achieved relative accuracy %10.3e\n",aprec))
                flush(printIOmon)
            end
        else
            if mPrMon >= 1
                write(printIOmon,"\n\n\n ",
                @sprintf("Solution of nonlinear system of equations obtained by NLEQ2\n"),
                @sprintf("No estimate available for the achieved relative accuracy\n"))
                flush(printIOmon)
            end
        end
    end
    # --------------------------------------------------------------------------
    # 9.2 Fail exit messages
    # --------------------------------------------------------------------------
    # 9.2.1 Termination at stationary point
    if retCode == 1 && mPrWarn == 1
        write(printIOwarn,"\nIteration terminated at stationary point\n")
        flush(printIOwarn)
    end
    # --------------------------------------------------------------------------
    # 9.2.2 Termination after more than nItmax iterations
    if retCode == 2 && mPrWarn == 1
        write(printIOwarn,"\n",
        @sprintf("Iteration terminates after nItmax %3i iteration steps\n",nItmax))
        flush(printIOwarn)
    end
    # --------------------------------------------------------------------------
    # 9.2.3 Newton method fails to converge
    if retCode == 3 && mPrWarn == 1
        write(printIOwarn,"\n",
        @sprintf("Newton method fails to converge\n"))
        flush(printIOwarn)
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
        if qMStop
            write(printIOwarn,"\nWARNING: Monotonicity test failed after ",ctyp,
            " convergence was already checked\nrTol requirement may be too",
            " stringent\n")
            flush(printIOwarn)
        else
            write(printIOwarn,"\nWARNING: ",ctyp, " convergence slowed down\n",
            "rTol requirement may be too stringent\n")
            flush(printIOwarn)
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
        flush(printIOwarn)
    end
    # --------------------------------------------------------------------------
    # 9.2.5 Error exit due to linear solver routine nFact
    if retCode == 80 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by linear solver nFact\n")
        flush(printIOwarn)
    end
    # --------------------------------------------------------------------------
    # 9.2.6 Error exit due to linear solver routine nSolv
    if retCode == 81 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by linear solver nSolv\n")
        flush(printIOwarn)
    end
    # --------------------------------------------------------------------------
    # 9.2.7 Error exit due to fail of user function fcn
    if retCode == 82 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by user function fcn\n")
        flush(printIOwarn)
    end
    # --------------------------------------------------------------------------
    # 9.2.8 Error exit due to fail of user function Jacobian
    if retCode == 83 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by user function jac\n")
        flush(printIOwarn)
    end
    if (retCode == 82 || retCode == 83) && nIter <= 1 && mPrWarn == 1
        write(printIOwarn,"Try to find a better initial guess for the solution\n")
        flush(printIOwarn)
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
        flush(printIOmon)
    end
    rTol = aprec
    sumX = sumXa
    if mPrSol >= 2 && nIter != 0
        nSout(n,xa,2,mPrSol,printIOsol,nIter,dLevF,sumX)
    elseif mPrSol >= 1 && nIter == 0
        nSout(n,xa,1,mPrSol,printIOsol,nIter,dLevF,sumX)
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
        nSout(n,x,modefi,mPrSol,printIOsol,nIter,dLevF,sumX)
    end
    # End of exits

    # --------------------------------------------------------------------------
    # 10 Prepare all the variables for returning
    xScal[:] = xw
    setOption!(opt, OPT_QSUCC, Int(qSucc))
    setOption!(opt, OPT_FCSTART, fc)
    setOption!(opt, OPT_IRANK, iRank)

    setOption!(wkNLEQ2, STATS_NITER,  nIter)
    setOption!(wkNLEQ2, STATS_NCORR,  nCorr)
    setOption!(wkNLEQ2, STATS_NFCN,   nFcn)
    setOption!(wkNLEQ2, STATS_NFCNJ,  nFcnJ)
    setOption!(wkNLEQ2, STATS_NJAC,   nJac)
    setOption!(wkNLEQ2, STATS_NREJR1, nRejR1)
    setOption!(wkNLEQ2, STATS_NEW,    newt)
    setOption!(wkNLEQ2, STATS_ICONV,  iConv)
    setOption!(wkNLEQ2, STATS_SUBCOND, cond1)
    setOption!(wkNLEQ2, STATS_SENS, sens1)
    setOption!(wkNLEQ2, STATS_IFAIL,  iFail)

    setOption!(wkNLEQ2, WK_SUMXA0, sumxa0)
    setOption!(wkNLEQ2, WK_SUMXA1, sumxa1)
    setOption!(wkNLEQ2, WK_FCMON, fcMon)
    setOption!(wkNLEQ2, WK_FCA, fcA)
    setOption!(wkNLEQ2, WK_FCKEEP, fcKeep)
    setOption!(wkNLEQ2, WK_FCPRI, fcPri)
    setOption!(wkNLEQ2, WK_DMYCOR, dMyCor)
    setOption!(wkNLEQ2, STATS_CONV, conv)
    setOption!(wkNLEQ2, STATS_SUMX, sumX)
    setOption!(wkNLEQ2, STATS_DLEVF, dLevF)

    setOption!(wkNLEQ2, "P_CLIN0", cLin0)
    setOption!(wkNLEQ2, "P_CLIN1", cLin1)
    setOption!(wkNLEQ2, "P_CALPHA", cAlpha)
    setOption!(wkNLEQ2, "P_ALPHAE", alphaE)
    setOption!(wkNLEQ2, "P_ALPHAK", alphaK)
    setOption!(wkNLEQ2, "P_ALPHAA", alphaA)
    setOption!(wkNLEQ2, "P_QMSTOP", qMStop)
    setOption!(wkNLEQ2, "P_SUMXA2", sumxa2)
    return retCode
    # End of function n2int
end
