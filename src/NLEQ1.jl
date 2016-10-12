"""
function nleq1(fcn, x::Vector{Float64}, xScal::Vector{Float64}, opt::OptionsNLEQ)

Damped Newton-algorithm for systems of highly nonlinear
equations - damping strategy due to Ref. (1).

(The iteration is done by function N1INT currently. NLEQ1
itself does some house keeping and builds up workspace.)

Jacobian approximation by numerical differences, user
supplied function JAC or forward mode automatic differentation.

The numerical solution of the arising linear equations is
done by means of Julia builtin lu-factorization routine in the
dense or sparse matrix case; or by the functions DGBFA and DGBSL,
which have been converted from the equal named Fortran LINPACK
subroutines into Julia code, in the band matrix case.
For special purposes these routines may be substituted.

This is a driver routine for the core solver N1INT.

## Input parameters
| Variable   | Description                                                                              |
|------------|------------------------------------------------------------------------------------------|
| fcn        | Function for which zero is to be found. Should be in the form of fcn(y,x) with y = f(x). |
| x[1:n]     | Initial estimate of the solution.                                                        |
| xScal[1:n] | User scaling (lower threshold) of the iteration vector x                                 |
| opt        | Options for solving the nonlinear system. Valid options are listed below.                |

## Output parameters
| Variable | Description                                                                                   |
|----------|-----------------------------------------------------------------------------------------------|
| x0[1:n]  | Solution values (or final values if exit before solution is reached).                         |
| stats    | A dictionary variable of additional output values. The fields are discussed below.            |
| retCode  | An integer value signifying the exit code. The meaning of the exit codes are discussed below. |

"""
function nleq1(fcn, x::Vector{Float64}, xScal::Vector{Float64}, opt::OptionsNLEQ)

    # Initialize error code 0
    retCode = 0

#-------------------------------------------------------------------------------
# Printing related stuff
#-------------------------------------------------------------------------------
    # Print warning messages?
    printWarn   = getOption!(opt,OPT_PRINTWARNING,0)
    # Print iteration summary?
    printMon    = getOption!(opt,OPT_PRINTITERATION,0)
    # Print solution summary?
    printSol    = getOption!(opt,OPT_PRINTSOLUTION,0)
    # Where to print?
    # Defaults to STDOUT
    printIOwarn = getOption!(opt,OPT_PRINTIOWARN,STDOUT)
    printIOmon  = getOption!(opt,OPT_PRINTIOMON,STDOUT)
    printIOsol  = getOption!(opt,OPT_PRINTIOSOL,STDOUT)

#-------------------------------------------------------------------------------

    # First call or successive call
    qSucc   = Bool(getOption!(opt,OPT_QSUCC,0))
    qIniMon = (printMon >= 1 && !qSucc)

    # Create an empty statistics variable to be returned
    stats = Dict{AbstractString,Any}()

    # Check input parameters and options
    n = length(x)
    retCode = checkOptions(n, x, xScal, opt)

    # Exit if any parameter error was detected
    if retCode != 0
        println("Exit with return code $retCode")
        return (x, stats, retCode)
    end

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    mStor = getOption!(opt,OPT_MSTOR,0)
    ml = getOption!(opt,"OPT_ML",0)
    mu = getOption!(opt,"OPT_MU",0)
    if mStor == 0
        m1 = n
        m2 = n
    elseif mStor == 1
        m1 = 2*ml + mu + 1
        m2 = ml + mu + 1
    end

    jacGen = opt.options[OPT_JACGEN]

    qRank1 = Bool(getOption!(opt, OPT_QRANK1, 0))
    qOrdi  = Bool(getOption!(opt, OPT_QORDI,  0))
    qSimpl = Bool(getOption!(opt, OPT_QSIMPL, 0))

    if qRank1
        nBroy = getOption!(opt,OPT_NBROY,0)
        if nBroy == 0
            nBroy = max(m2,10)
            setOption!(opt,OPT_NBROY, nBroy)
        end
    else
        nBroy = 0
    end

    # Check if this is a first call or successive call to nleq1
    # If first call then reset the workspace and persistent variables
    if !qSucc
        empty!(wkNLEQ1.options)
        initializeOptions(opt, wkNLEQ1, n, m1, nBroy, qRank1, 1)
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
        printInitialization(n, printIOmon, opt.options[OPT_RTOL], jacGen, mStor,
        ml, mu, opt.options[OPT_NOROWSCAL], qRank1, nonLin, qBDamp,
        opt.options[OPT_FCBAND], qOrdi, qSimpl, nItmax)
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

    # Simplified Newton iteration implies ordinary Newton iteration mode
    if qSimpl
        setOption!(opt, OPT_QORDI, 1)
    end

    # If ordinary Newton iteration, damping factor is always 1
    if opt.options[OPT_QORDI] == 1
        fc = 1.0
    end

    # Set starting damping factor
    setOption!(opt, OPT_FCSTART, fc)

    if printMon >= 2 && !qSucc
        write(printIOmon,"\nINFO: ","Internal parameters:",
        "\n\tStarting value for damping factor ",
        @sprintf("OPT_FCSTART\t= %1.2e",opt.options[OPT_FCSTART]),
        @sprintf("\n\tMinimum allowed damping factor OPT_FCMIN\t= %1.2e",fcMin),
        "\n\tRank-1 updates decision parameter ",
        @sprintf("OPT_SIGMA\t= %1.2e\n",opt.options[OPT_SIGMA]))
    end

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    retCode = -1

    if nBroy == 0
        nBroy = 1
    end

    # Create a copy inside so that the original variable is untouched
    x0 = x[:]

    # Call to n1int
    retCode = n1int(n, fcn, x0, xScal, opt.options[OPT_RTOL], nItmax,
        nonLin, opt, m1, m2, nBroy, opt.options[OPT_FCSTART], opt.options[OPT_FCMIN],
        opt.options[OPT_SIGMA], opt.options[OPT_SIGMA2], mStor, printWarn,
        printMon, printSol, printIOwarn, printIOmon, printIOsol, qBDamp)

    # set stats variable
    stats[STATS_XSCAL] = xScal
    if retCode == -1
        stats[STATS_RTOL] = wkNLEQ1.options["P_TOLALL"][wkNLEQ1.options[STATS_NITER]]
    else
        stats[STATS_RTOL] = opt.options[OPT_RTOL]
    end
    if opt.options[OPT_STORE] == 1
        stats[STATS_XITER]      = wkNLEQ1.options["P_XITER"]
        stats[STATS_NATLEVEL]   = wkNLEQ1.options["P_SUMXALL"]
        stats[STATS_SIMLEVEL]   = wkNLEQ1.options["P_DLEVFALL"]
        stats[STATS_STDLEVEL]   = wkNLEQ1.options["P_SUMXQALL"]
        stats[STATS_DAMPINGFC]  = wkNLEQ1.options["P_FCALL"]
    end
    stats[STATS_PRECISION]  = wkNLEQ1.options["P_TOLALL"]
    stats[STATS_NITER]      = wkNLEQ1.options[STATS_NITER]
    stats[STATS_NCORR]      = wkNLEQ1.options[STATS_NCORR]
    stats[STATS_NREJR1]     = wkNLEQ1.options[STATS_NREJR1]
    stats[STATS_NJAC]       = wkNLEQ1.options[STATS_NJAC]
    stats[STATS_NFCN]       = wkNLEQ1.options[STATS_NFCN]
    stats[STATS_NFCNJ]      = wkNLEQ1.options[STATS_NFCNJ]
    stats[STATS_IFAIL]      = wkNLEQ1.options[STATS_IFAIL]

    # Print statistics
    if printMon >= 2 && retCode != -1 && retCode != 10
        printStats(stats, printIOmon)
    end

    return (x0, stats, retCode);
end

"""
function n1int(n::Int64, fcn, x::Vector{Float64}, xScal::Vector{Float64},
    rTol::Float64, nItmax::Int64, nonLin::Int64, opt::OptionsNLEQ,
    m1::Int64, m2::Int64, nBroy::Int64, fc::Float64, fcMin::Float64,
    sigma::Float64, sigma2::Float64, mStor::Int64, mPrWarn::Int64, mPrMon::Int64,
    mPrSol::Int64, printIOwarn, printIOmon, printIOsol, qBDamp::Bool)

Core routine for NLEQ1.

Damped Newton-algorithm for systems of highly nonlinear
equations especially designed for numerically sensitive
problems.

## Input parameters
| Variable     | Description                                                                              |
|--------------|------------------------------------------------------------------------------------------|
| n            | Number of vector components                                                              |
| fcn          | Function for which zero is to be found. Should be in the form of fcn(y,x) with y = f(x). |
| x[1:n]*      | Initial estimate of the solution.                                                        |
| xScal[1:n]   | User scaling (lower threshold) of the iteration vector x                                 |
| rTol         | Required relative precision of solution components                                       |
| nItmax       | Maximum number of allowed iterations                                                     |
| nonLin       | Problem type specification. See OPT_NONLIN field in solver options                       |
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
|----------|-----------------------------------------------------------------------------------------------|
| x[1:n]*  | Solution values (or final values if exit before solution is reached).                         |
| retCode  | An integer value signifying the exit code. The meaning of the exit codes are discussed below. |
"""
function n1int(n::Int64, fcn, x::Vector{Float64}, xScal::Vector{Float64},
    rTol::Float64, nItmax::Int64, nonLin::Int64, opt::OptionsNLEQ,
    m1::Int64, m2::Int64, nBroy::Int64, fc::Float64, fcMin::Float64,
    sigma::Float64, sigma2::Float64, mStor::Int64, mPrWarn::Int64, mPrMon::Int64,
    mPrSol::Int64, printIOwarn, printIOmon, printIOsol, qBDamp::Bool)

    # --------------------------------------------------------------------------
    # Since wkNLEQ1 is module global
    # Create the local variables here rather than taking them as arguments
    a        = wkNLEQ1.options[WK_A]
    dxSave   = wkNLEQ1.options[WK_DXSAVE]
    dx       = wkNLEQ1.options[WK_DX]
    dxQ      = wkNLEQ1.options[WK_DXQ]
    xa       = wkNLEQ1.options[WK_XA]
    xwa      = wkNLEQ1.options[WK_XWA]
    f        = wkNLEQ1.options[WK_F]
    fa       = wkNLEQ1.options[WK_FA]
    eta      = wkNLEQ1.options[WK_ETA]
    xw       = wkNLEQ1.options[WK_XW]
    fw       = wkNLEQ1.options[WK_FW]
    dxQa     = wkNLEQ1.options[WK_DXQA]
    sumxa0   = wkNLEQ1.options[WK_SUMXA0]
    sumxa1   = wkNLEQ1.options[WK_SUMXA1]
    fcMon    = wkNLEQ1.options[WK_FCMON]
    fcA      = wkNLEQ1.options[WK_FCA]
    fcKeep   = wkNLEQ1.options[WK_FCKEEP]
    fcPri    = wkNLEQ1.options[WK_FCPRI]
    dMyCor   = wkNLEQ1.options[WK_DMYCOR]
    conv     = wkNLEQ1.options[STATS_CONV]
    sumXs    = wkNLEQ1.options[WK_SUMXS]
    dLevF    = wkNLEQ1.options[STATS_DLEVF]
    nIter    = wkNLEQ1.options[STATS_NITER]
    nCorr    = wkNLEQ1.options[STATS_NCORR]
    nFcn     = wkNLEQ1.options[STATS_NFCN]
    nFcnJ    = wkNLEQ1.options[STATS_NFCNJ]
    nJac     = wkNLEQ1.options[STATS_NJAC]
    nRejR1   = wkNLEQ1.options[STATS_NREJR1]
    newt     = wkNLEQ1.options[STATS_NEW]
    iConv    = wkNLEQ1.options[STATS_ICONV]
    tolAll   = wkNLEQ1.options["P_TOLALL"]

    if opt.options[OPT_STORE] == 1
        xIter    = wkNLEQ1.options["P_XITER"]
        sumXall  = wkNLEQ1.options["P_SUMXALL"]
        dLevFall = wkNLEQ1.options["P_DLEVFALL"]
        sumXQall = wkNLEQ1.options["P_SUMXQALL"]
        fcAll    = wkNLEQ1.options["P_FCALL"]
    end
    # --------------------------------------------------------------------------
    # 0.1 Variables that need to be defined before since they appear in different
    # scopes. The declaration and usage are in different scopes.
    retCode = -1
    dLevFn  = 0.0
    sumXa   = 0.0
    conva   = 0.0
    # --------------------------------------------------------------------------
    # 0.2 Persistent variables
    cLin0   = getOption!(wkNLEQ1,"P_CLIN0",0.0)
    cLin1   = getOption!(wkNLEQ1,"P_CLIN1",0.0)
    cAlpha  = getOption!(wkNLEQ1,"P_CALPHA",0.0)
    alphaE  = getOption!(wkNLEQ1,"P_ALPHAE",0.0)
    alphaK  = getOption!(wkNLEQ1,"P_ALPHAK",0.0)
    alphaA  = getOption!(wkNLEQ1,"P_ALPHAA",0.0)
    qMStop  = getOption!(wkNLEQ1,"P_QMSTOP",false)
    sumxa2  = getOption!(wkNLEQ1,"P_SUMXA2",0.0)
    l       = getOption!(wkNLEQ1,"P_L",zero(a))
    u       = getOption!(wkNLEQ1,"P_U",zero(a))
    p       = getOption!(wkNLEQ1,"P_P",zeros(Int64,n))
    # --------------------------------------------------------------------------
    # Begin
    # --------------------------------------------------------------------------
    # 1 Initialization
    # --------------------------------------------------------------------------
    # 1.1 Control variables
    qSucc       = Bool(opt.options[OPT_QSUCC])
    qScale      = opt.options[OPT_NOROWSCAL] != 1
    qOrdi       = Bool(opt.options[OPT_QORDI])
    qSimpl      = Bool(opt.options[OPT_QSIMPL])
    qRank1      = Bool(opt.options[OPT_QRANK1])
    iOrMon      = getOption!(opt, OPT_IORMON, 2)
    iScal       = getOption!(opt, OPT_ISCAL, 0)
    mode        = getOption!(opt, OPT_MODE,  0)
    jacGen      = opt.options[OPT_JACGEN]
    qMixIO      = typeof(printIOmon) == typeof(printIOsol)
    if qMixIO && typeof(printIOmon) == IOStream && printIOmon.name != printIOsol.name
        qMixIO = false
    end
    qMixIO      &= mPrMon != 0 && mPrSol != 0
    qLU         = !qSimpl
    # --------------------------------------------------------------------------
    # 1.2 Derived dimensional parameters
    if mStor == 0
        ml = 0
        mu = 0
    elseif mStor == 1
        ml = m1 - m2
        mu = m2 - 1 - ml
    end
    # --------------------------------------------------------------------------
    # 1.3 Derived internal parameters
    fcMin2  = fcMin*fcMin
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
        if opt.options[OPT_STORE] == 1
            push!(xIter,x)
        end

        nIter  = 0
        wkNLEQ1.options[STATS_NITER] = nIter
        nCorr  = 0
        nRejR1 = 0
        nFcn   = 0
        nFcnJ  = 0
        nJac   = 0
        iConv  = 0
        conv   = 0.0

        fcKeep  = fc
        fcA     = fc
        fcPri   = fc
        fcMon   = fc
        sumxa0  = 0.0
        sumxa1  = 0.0

        qGenJ   = true
        fcK2    = fc

        if jacGen == 3
            eta[:] = etaIni*ones(n)
        end

        xa[:] = x

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
        catch
            iFail = -1
            retCode = 82
            qIter   = false
        end
        nFcn += 1
        if length(f) != n
            retCode = 22
            qIter   = false
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
                # --------------------------------------------------------------
                # 2.2 Aposteriori estimate of damping factor
                dxQa[:] = dxQ
                if !qOrdi
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
        end
        qJcRfr = false
        # ----------------------------------------------------------------------
        # 2.3 Jacobian matrix
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
                catch
                    iFail   = -1
                end
            else
                if mStor == 0
                    if jacGen == 2
                        (nFcnJ,iFail) = nJacFD(fcn,n,n,x,f,xw,aJdel,aJmin,nFcnJ,a)
                    end
                    if jacGen == 3
                        (nFcnJ,iFail) = nJcf(fcn,n,n,x,f,xw,eta,etaMin,
                                            etaMax,etaDif,conv,nFcnJ,a)
                    end
                    # Forward mode automatic differentiation
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
                                iFail = 0
                            catch
                                iFail = -1
                            end
                        end
                    end
                elseif mStor == 1
                    if jacGen == 2
                        (nFcnJ,iFail) = nJacFDb(fcn,n,m1,ml,x,f,xw,aJdel,aJmin,
                                            nFcnJ,a)
                    end
                    if jacGen == 3
                        (nFcnJ,iFail) =
                            nJcfb(fcn,n,m1,ml,x,f,xw,eta,etaMin,etaMax,etaDif,conv,nFcnJ,a)
                    end
                    # Forward mode automatic differentiation
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
                                iFail = 0
                            catch
                                iFail = -1
                            end
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

        elseif !qSimpl
            newt += 1
        end

        if newt == 0 && (qLU || nIter == 0)
            # ------------------------------------------------------------------
            # 2.3.2.1 Save scaling values
            xwa[:] = xw
            # ------------------------------------------------------------------
            if issparse(a)
                nza = nnz(a)
                (row,col) = findn(a)
                for k = 1:nza
                    a[row[k],col[k]] = -a[row[k],col[k]]*xw[col[k]]
                end
            else
                # --------------------------------------------------------------
                # 2.3.2.2 Prepare Jacobian for use by band-solver
                if mStor == 1
                    for l1 = 1:n
                        for l2 = m2:-1:1
                            a[l2+ml,l1] = a[l2,l1]
                        end
                    end
                end
                # --------------------------------------------------------------
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
            end
            # ------------------------------------------------------------------
            # 2.4.2 Row scaling of matrix A
            if qScale
                if mStor == 0
                    nScrf(n,n,a,fw)
                elseif mStor == 1
                    nScrb(n,m1,ml,mu,a,fw)
                end
            else
                fw[:] = 1.0
            end
        end
        # ----------------------------------------------------------------------
        # 2.4.3 Save and scale values of F(n)
        fa[:] = f
        t1 = f.*fw
        # ----------------------------------------------------------------------
        # 3 Central part of iteration step
        # ----------------------------------------------------------------------
        # 3.1 Solution of the linear system
        # ----------------------------------------------------------------------
        # 3.1.1 Decomposition of (n,n) matrix A
        if newt == 0 && (qLU || nIter == 0)
            iFail = nFact(n,m1,ml,mu,a,mStor,l,u,p)
            if iFail != 0
                if iFail == 1
                    retCode = 1
                else
                    retCode = 80
                end
                break
            end
        end
        # ----------------------------------------------------------------------
        # 3.1.2 Solution of (n,n) system
        if newt == 0
            iFail = nSolv(n,m1,ml,mu,l,u,p,t1,mStor)
            if iFail != 0
                retCode = 81
                break
            end
        else
            alfa1 = sum(dx.*dxQ./xw.^2)
            alfa2 = sum(dx.^2./xw.^2)
            alfa = alfa1/alfa2
            beta = 1.0 - alfa
            t1 = (dxQ+(fcA-one)*alfa*dx)/beta
            if newt == 1
                dxSave[1:n,1] = dx
            end
            dxSave[1:n,newt+1:newt+1] = t1
            dx[:] = t1
            t1 = t1./xw
        end
        # ----------------------------------------------------------------------
        # 3.2 Evaluation of scaled natural level function sumX
        # scaled maximum error norm conv
        # evaluation of (scaled) standard level function dlevf
        # and computation of ordinary Newton corrections dx[n]
        if !qSimpl
            (conv,sumX,dLevF) = nLvls(n,dx,t1,xw,f,newt == 0)
        else
            (conv,sumX,dLevF) = nLvls(n,dx,t1,xwa,f,newt == 0)
        end
        wkNLEQ1.options[STATS_SUMX]   = sumX
        wkNLEQ1.options[STATS_DLEVF]  = dLevF
        xa[:]    = x
        sumXa    = sumX
        dLevXa   = sqrt(sumXa/n)
        conva    = conv
        dxANrm   = wnorm(n,dx,xw)
        if opt.options[OPT_STORE] == 1
            push!(sumXall,dLevXa)
            push!(dLevFall,dLevF)
        end

        # ----------------------------------------------------------------------
        # 3.3 A - priori estimate of damping factor fc
        if nIter != 0 && nonLin != 1 && newt == 0 && !qOrdi
            # 3.3.1 Computation of the denominatior of a-priori estimate
            fcDnm = sum(((dx-dxQa)./xw).^2)*sumX
            # ------------------------------------------------------------------
            # 3.3.2 New damping factor
            if fcDnm > fcNumP*fcMin2 || (nonLin == 4 && fcA^2*fcNumP < 4.0*fcDnm)
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

            fc = max(fcPri,fcMin)
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

        if iOrMon >= 2
            sumxa2 = sumxa1
            sumxa1 = sumxa0
            sumxa0 = dLevXa
            if sumxa0 == 0.0
                sumxa0 = small
            end
            # Check convergence rate (linear or superlinear)
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
        # ----------------------------------------------------------------------
        # 3.4 Save natural level for later computations of corrector
        # and print iterate
        fcNumK = sumX
        if mPrMon >= 2
            nPrv1(dLevF,dLevXa,fcKeep,nIter,newt,mPrMon,printIOmon,qMixIO)
        end
        nRed    = 0
        qNext   = false
        qRep    = false
        qRed    = true
        iCnv    = 0

        # Damping-factor reduction loop
        # ======================================================================
        while qRed
            # ------------------------------------------------------------------
            # 3.5 Preliminary new iterate
            x[:] = xa + dx*fc
            if opt.options[OPT_STORE] == 1
                push!(fcAll,fc)
            end
            # ------------------------------------------------------------------
            # 3.5.2 Exit, if problem is specified as being linear
            if nonLin == 1
                retCode = 0
                break
            end
            #-------------------------------------------------------------------
            # 3.6.1 Computation of the residual vector
            try
                fcn(f,x)
            catch
                iFail = -1
            end
            nFcn += 1
            # TODO: Understand what is happening here
            # and handle the failure properly
            # What does iFail = 1 and iFail = 2 mean??
            if iFail < 0
                retCode = 82
                break
            end
            if iFail == 1 || iFail == 2
                if qOrdi
                    retCode = 82
                    break
                end
                if iFail == 1
                    fcRedu = 0.5
                else
                    fcRedu = f[1]

                    if fcRedu <= 0 || fcRedu >= 1
                        retCode = 82
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
                if qOrdi
                    # ----------------------------------------------------------
                    # 3.6.2 Convergence test for ordinary Newton iteration
                    push!(tolAll,dxANrm)
                    if dxANrm <= rTol
                        retCode = 0
                        iCnv    = 1
                        break
                    end
                else
                    t1 = f.*fw
                    # ------------------------------------------------------
                    # 3.6.3 Solution of linear (n,n) system
                    iFail = nSolv(n,m1,ml,mu,l,u,p,t1,mStor)
                    if iFail != 0
                        retCode = 81
                        break
                    end
                    if newt > 0
                        dxQ[:] = t1.*xwa
                        for iLoop = 1:newt
                            sum1 = sum((dxQ.*dxSave[1:n,iLoop:iLoop])./xw.^2)
                            sum2 = sum((dxSave[1:n,iLoop:iLoop]./xw).^2)
                            beta = sum1/sum2
                            dxQ[:] = dxQ + beta*dxSave[1:n,iLoop+1:iLoop+1]
                            t1 = dxQ./xw
                        end
                    end
                    # ------------------------------------------------------
                    # 3.6.4 Evaluation of scaled natural level function
                    #       sumX
                    #       scaled maximum error norm conv and evaluation
                    #       of (scaled) standard level function dLevFn
                    if !qSimpl
                        (conv,sumX,dLevFn) =
                        nLvls(n,dxQ,t1,xw,f,newt==0)
                    else
                        (conv,sumX,dLevFn) =
                        nLvls(n,dxQ,t1,xwa,f,newt==0)
                    end
                    dxNrm = wnorm(n,dxQ,xw)
                    push!(tolAll,dxNrm)
                    if opt.options[OPT_STORE] == 1
                        push!(sumXQall,sqrt(sumX/n))
                    end
                    # ------------------------------------------------------
                    # 3.6.5 Convergence test
                    if dxNrm <= rTol && dxANrm <= rSmall && fc == 1.0
                        retCode = 0
                        iCnv = 1
                        break
                    end

                    fcA = fc
                    # ------------------------------------------------------
                    # 3.6.6 Evaluation of reduced damping factor
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
                end
                # ----------------------------------------------------------
                # 3.7 Natural monotonicity test
                if sumX > sumXa && !qOrdi
                    # ------------------------------------------------------
                    # 3.8 Output of iterate
                    if mPrMon >= 3
                        nPrv2(dLevFn,sqrt(sumX/n),fc,nIter,printIOmon,qMixIO,"*")
                    end
                    if qMStop
                        retCode = 4
                        break
                    end
                    fch = min(fcCor,0.5*fc)
                    if fc > fcMin
                        fc = max(fch,fcMin)
                    else
                        fc = fch
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
                    # ------------------------------------------------------
                    # 3.10 If damping factor is too small:
                    #      Refreash Jacobian, if current Jacobian was computed
                    #      by a Rank-1 update, else fail and exit
                    qJcRfr = fc < fcMin || (newt > 0 && nRed > 1)

                    if qJcRfr && newt == 0
                        retCode = 3
                        break
                    end
                else
                    if !qOrdi && !qRep && fcCor > sigma2*fc
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
                        end

                        qRep = true
                    else
                        qNext = true
                    end
                end
            end
            qRed = !(qNext||qJcRfr)
        end
        # End of damping factor reduction loop
        if nonLin == 1 || iCnv == 1 || (retCode != 0 && retCode != -1)
            break
        end
        # ======================================================================
        if qJcRfr
            # ------------------------------------------------------------------
            # 3.11 Restore former values for repeating iteration step
            nRejR1 += 1
            x[:] = xa
            f[:] = fa
            if mPrMon >= 2
                write(printIOmon,
                @sprintf("        %2i not accepted damping factor %7.5f",nIter,fc),
                @sprintf("    %2i\n",newt))
            end
            fc  = fcKeep
            fcA = fck2
            if nIter == 0
                fc = fcMin
            end
            qGenJ = true
        else
            # ------------------------------------------------------------------
            # 4 Preparations to start the following iteration step
            # ------------------------------------------------------------------
            # Print values
            if mPrMon >= 3 && !qOrdi
                nPrv2(dLevFn,sqrt(sumX/n),fc,nIter+1,printIOmon,qMixIO,"*")
            end
            # Print the natural level of the current iterate and
            # return it in one-step mode
            sumXs = sumX
            sumX = sumXa
            if mPrSol >= 2 && nIter != 0
                nSout(n,xa,2,mPrSol,printIOsol,nIter,dLevF,sumX)
            elseif mPrSol >= 1 && nIter == 0
                nSout(n,xa,1,mPrSol,printIOsol,nIter,dLevF,sumX)
            end
            nIter += 1
            wkNLEQ1.options[STATS_NITER] = nIter
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
            if mode == 1
                qSucc = true
                setOption!(opt, OPT_QSUCC, Int(qSucc))
                setOption!(opt, OPT_FCSTART, fc)

                setOption!(wkNLEQ1, STATS_NITER,  nIter)
                setOption!(wkNLEQ1, STATS_NCORR,  nCorr)
                setOption!(wkNLEQ1, STATS_NFCN,   nFcn)
                setOption!(wkNLEQ1, STATS_NFCNJ,  nFcnJ)
                setOption!(wkNLEQ1, STATS_NJAC,   nJac)
                setOption!(wkNLEQ1, STATS_NREJR1, nRejR1)
                setOption!(wkNLEQ1, STATS_NEW,    newt)
                setOption!(wkNLEQ1, STATS_ICONV,  iConv)
                setOption!(wkNLEQ1, STATS_IFAIL,  iFail)

                setOption!(wkNLEQ1, WK_SUMXA0, sumxa0)
                setOption!(wkNLEQ1, WK_SUMXA1, sumxa1)
                setOption!(wkNLEQ1, WK_FCMON, fcMon)
                setOption!(wkNLEQ1, WK_FCA, fcA)
                setOption!(wkNLEQ1, WK_FCKEEP, fcKeep)
                setOption!(wkNLEQ1, WK_FCPRI, fcPri)
                setOption!(wkNLEQ1, WK_DMYCOR, dMyCor)
                setOption!(wkNLEQ1, STATS_CONV, conv)
                setOption!(wkNLEQ1, STATS_SUMX, sumX)
                setOption!(wkNLEQ1, WK_SUMXS, sumXs)
                setOption!(wkNLEQ1, STATS_DLEVF, dLevF)

                setOption!(wkNLEQ1, "P_CLIN0", cLin0)
                setOption!(wkNLEQ1, "P_CLIN1", cLin1)
                setOption!(wkNLEQ1, "P_CALPHA", cAlpha)
                setOption!(wkNLEQ1, "P_ALPHAE", alphaE)
                setOption!(wkNLEQ1, "P_ALPHAK", alphaK)
                setOption!(wkNLEQ1, "P_ALPHAA", alphaA)
                setOption!(wkNLEQ1, "P_QMSTOP", qMStop)
                setOption!(wkNLEQ1, "P_SUMXA2", sumxa2)

                return retCode
            end
        end
    end
    # End repeat
    # End of main iteration loop
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # 9 Exits
    # --------------------------------------------------------------------------
    # 9.1 Solution exit
    aprec = -1.0

    if retCode == 0 || retCode == 4
        if nonLin != 1
            if !qOrdi
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
                # Print final monitor output
                if mPrMon >= 2
                    if retCode == 0
                        nPrv2(dLevFn,sqrt(sumX/n),fc,nIter+1,
                        printIOmon,qMixIO,"*")
                    elseif iOrMon == 3
                        nPrv1(dLevFn,sqrt(sumXa/n),fc,nIter,newt,
                        mPrMon,printIOmon,qMixIO)
                    end
                end
                if iOrMon >= 2
                    if iConv <= 1 && alphaE != 0.0 && alphaK != 0.0
                        retCode = 5
                    end
                end
            else
                # if qOrdi is true
                aprec = sqrt(sumXa/n)
            end
            if mPrMon >= 1
                if qOrdi || retCode == 4
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
                @sprintf("Solution of nonlinear system of equations obtained by NLEQ1\n"),
                @sprintf("No estimate available for the achieved relative accuracy\n"))
            end
        end
    end
    # --------------------------------------------------------------------------
    # 9.2 Fail exit messages
    # --------------------------------------------------------------------------
    # 9.2.1 Termination, since Jacobian matrix became singular
    if retCode == 1 && mPrWarn == 1
        write(printIOwarn,"\nIteration terminated due to singular Jacobian matrix\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.2 Termination after more than nItmax iterations
    if retCode == 2 && mPrWarn == 1
        write(printIOwarn,"\n",
        @sprintf("Iteration terminates after nItmax %3i iteration steps\n",nItmax))
    end
    # --------------------------------------------------------------------------
    # 9.2.3 Damping factor fc became too small
    if retCode == 3 && mPrWarn == 1
        write(printIOwarn,"\n",
        @sprintf("Damping factor has become too small: lambda = %10.3e\n",fc))
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
    # 9.2.5
    if retCode == 22 && mPrWarn == 1
        write(printIOwarn,"ERROR: dimensions of startvector and problem function ",
        "output differ:\n",
        @sprintf("      length(x0) = %5i,     length(fcn(x0)) = %5i\n",n,length(f)))
    end
    # --------------------------------------------------------------------------
    # 9.2.6 Error exit due to linear solver routine nFact
    if retCode == 80 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by linear solver nFact\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.7 Error exit due to linear solver routine nSolv
    if retCode == 81 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by linear solver nSolv\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.8 Error exit due to fail of user function fcn
    if retCode == 82 && mPrWarn == 1
        write(printIOwarn,@sprintf("ERROR: %5i",iFail)," signalled by user function fcn\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.9 Error exit due to fail of user function Jacobian
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
    rTol = aprec
    sumX = sumXa
    if mPrSol >= 2 && nIter != 0
        mode = 2
        if qOrdi
            mode = 3
        end
        nSout(n,xa,mode,mPrSol,printIOsol,nIter,dLevF,sumX)
    elseif mPrSol >= 1 && nIter == 0
        nSout(n,xa,1,mPrSol,printIOsol,nIter,dLevF,sumX)
    end
    if !qOrdi
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
    end
    # End of exits

    # --------------------------------------------------------------------------
    # 10 Prepare all the variables for returning
    xScal[:] = xw

    setOption!(opt, OPT_QSUCC, Int(qSucc))
    setOption!(opt, OPT_FCSTART, fc)

    setOption!(wkNLEQ1, STATS_NITER,  nIter)
    setOption!(wkNLEQ1, STATS_NCORR,  nCorr)
    setOption!(wkNLEQ1, STATS_NFCN,   nFcn)
    setOption!(wkNLEQ1, STATS_NFCNJ,  nFcnJ)
    setOption!(wkNLEQ1, STATS_NJAC,   nJac)
    setOption!(wkNLEQ1, STATS_NREJR1, nRejR1)
    setOption!(wkNLEQ1, STATS_NEW,    newt)
    setOption!(wkNLEQ1, STATS_ICONV,  iConv)
    setOption!(wkNLEQ1, STATS_IFAIL,  iFail)

    setOption!(wkNLEQ1, WK_SUMXA0, sumxa0)
    setOption!(wkNLEQ1, WK_SUMXA1, sumxa1)
    setOption!(wkNLEQ1, WK_FCMON, fcMon)
    setOption!(wkNLEQ1, WK_FCA, fcA)
    setOption!(wkNLEQ1, WK_FCKEEP, fcKeep)
    setOption!(wkNLEQ1, WK_FCPRI, fcPri)
    setOption!(wkNLEQ1, WK_DMYCOR, dMyCor)
    setOption!(wkNLEQ1, STATS_CONV, conv)
    setOption!(wkNLEQ1, STATS_SUMX, sumX)
    setOption!(wkNLEQ1, WK_SUMXS, sumXs)
    setOption!(wkNLEQ1, STATS_DLEVF, dLevF)

    setOption!(wkNLEQ1, "P_CLIN0", cLin0)
    setOption!(wkNLEQ1, "P_CLIN1", cLin1)
    setOption!(wkNLEQ1, "P_CALPHA", cAlpha)
    setOption!(wkNLEQ1, "P_ALPHAE", alphaE)
    setOption!(wkNLEQ1, "P_ALPHAK", alphaK)
    setOption!(wkNLEQ1, "P_ALPHAA", alphaA)
    setOption!(wkNLEQ1, "P_QMSTOP", qMStop)
    setOption!(wkNLEQ1, "P_SUMXA2", sumxa2)

    return retCode
    # End of function n1int
end
