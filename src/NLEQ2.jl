include("NLEQ2Main.jl")

function nleq2(fcn::Function, x, xScal, opt::OptionsNLEQ)

    # TODO: Get rid of this assertion.
    assert(typeof(x[1]) == Float64 && typeof(xScal[1]) == Float64)

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
    mStor = getOption!(opt,OPT_MSTOR,0)
    if mStor == 0
        m1 = n
        m2 = n
    elseif mStor == 1
        ml = getOption!(opt,"OPT_ML",0)
        mu = getOption!(opt,"OPT_MU",0)
        m1 = 2*ml + mu + 1
        m2 = ml + mu + 1
    end

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

    # Check if this is a first call or successive call to nleq1
    # If first call then reset the workspace and persistent variables
    if !qSucc
        empty!(wkNLEQ2.options)
        empty!(xIter)
        empty!(sumXall)
        empty!(dLevFall)
        empty!(sumXQall)
        empty!(tolAll)
        empty!(fcAll)
        initializeOptions(opt, wkNLEQ2, n, m1, qRank1)
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

    epMach = getMachineConstants(3)
    cond = getOption!(opt, OPT_COND, 1/epMach)
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
    nonLin, iRank, cond, opt, retCode, m1, m2, nBroy,
    xIter, sumXall, dLevFall, sumXQall, tolAll, fcAll,
    opt.options[OPT_FCSTART], opt.options[OPT_FCMIN], opt.options[OPT_SIGMA],
    opt.options[OPT_SIGMA2], mStor, printWarn,
    printMon, printSol, printIOwarn, printIOmon, printIOsol, qBDamp)

    # set stats variable
    stats = Dict{ASCIIString,Any}()
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

    return (x, stats, retCode);
end
