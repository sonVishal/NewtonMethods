# TODO: Make everything type independent
# Currently everything is assumed to be Float64

# function nleq1(fcn::Function,x::Vector,xScal::Vector,opt::OptionsNLEQ)
#     Input parameters:
#         fcn:        Function of the type fcn(x,f,retCode) -> nothing
#         x:          Initial guess
#         xScal:      Scaling vector
#         opt:        Options for the solver
#         wk:         Structured workspace TODO: get rid of this
#     Output parameters:
#         x:          Solution if found
#         xScal:      Lat internal scaling vector if successful solve
#         stats:      Statistics of the solution
#         retCode:    Return code signifying success or failure
# This function does the bookkeeping. The actual computation are done in n1int

include("CheckOptionsNLEQ1.jl")
include("Jacobian.jl")
include("Error.jl")
function nleq1(fcn::Function,x::Vector,xScal::Vector,opt::OptionsNLEQ,
    wk::OptionsNLEQ)

    # TODO: persistent variables
    # Might be a good idea to store them inside opt
    # so that they can be used again and again.

    # Initialize output statistics
    stats = Dict{ASCIIString,Any}()

    # Initialize a common message string variable
    message = ""

    # Initialize error code 0
    retCode = 0

#-------------------------------------------------------------------------------
# Printing related stuff
#-------------------------------------------------------------------------------
    # Print warning messages?
    printFlag   = getOption!(opt,OPT_PRINTWARNING,0)
    # Print iteration summary?
    printIt     = getOption!(opt,OPT_PRINTITERATION,0)
    # Print solution summary?
    printSol    = getOption!(opt,OPT_PRINTSOLUTION,0)
    # Where to print?
    # Defaults to STDOUT
    printIO     = getOption!(opt,OPT_PRINTIO,STDOUT)
    if printIO == "FILE"
        # If not STDOUT then print to file
        # Default file name is log.txt and the file is opened for writing
        printFileName   = getOption!(opt,OPT_PRINTFILENAME,"log.txt")
        printFileMode   = getOption!(opt,OPT_PRINTFILEMODE,"w")
        if printFileMode != "w" || printFileMode != "a"
            throw(InvalidOption("OPT_PRINTFILEMODE",printFileMode))
        end
        f = open(printFileName,printFileMode)
    end
#-------------------------------------------------------------------------------

    # First call or successive call
    qSucc   = getOption!(opt,OPT_QSUCC,0)
    qIniMon = (printIt >= 1 && qSucc == 0)

    # Check input parameters and options
    n = length(x)
    retCode = checkOptions(n,x,xScal,opt)

    # return (stats,retCode)

    # Exit if any parameter error was detected
    if retCode != 0
        error("Exit with return code $retCode")
    end

    # Check if this is a first call or successive call
    # to nleq1
    if qSucc == 0
        # If this is the first call then assign memory to the variables
        xIter       = []
        sumXall     = []
        dLevFall    = []
        sumXQall    = []
        tolAll      = []
        fcAll       = []
    end

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    mStor = getOption!(opt,OPT_MSTOR,0)
    if mStor == 0
        m1 = n
        m2 = n
    elseif mStor == 1
        ml = getOption!(opt,"ML",0)
        mu = getOption!(opt,"MU",0)
        m1 = 2*ml + mu + 1
        m2 = ml + mu + 1
    end

    # Assign the Jacobian depending on user input
    # Multiple dispatch calls the required function based on
    # the storage requirement of the user
    jacGen = getOption!(opt,OPT_JACGEN,0)
    if jacGen == 0
        jacGen = 2
        opt.options[OPT_JACGEN] = jacGen;
    end

    qRank1 = getOption!(opt, OPT_QRANK1,    0)
    qOrdi  = getOption!(opt, OPT_QORDI,     0)
    qSimpl = getOption!(opt, OPT_QSIMPL,    0)

    if qRank1 == 1
        nBroy = getOption!(opt,OPT_NBROY,0)
        if nBroy == 0
            nBroy = max(m2,10)
            opt.options[OPT_NBROY] = nBroy
        end
    else
        nBroy = 0
    end

    # Initialize options
    initOption!(opt, OPT_FCMIN,     0.0)
    initOption!(opt, OPT_SIGMA,     0.0)
    initOption!(opt, OPT_SIGMA2,    0.0)
    initOption!(opt, OPT_NOROWSCAL, 0)

    # Workspace: WK
    initOption!(wk,WK_A,0.0)

    if qRank1 == 1
        initOption!(wk,WK_DXSAVE,zeros(n,nBroy))
    else
        initOption!(wk,WK_DXSAVE,0.0)
    end

    # Initialize temporary workspace
    initOption!(wk, WK_DX  , zeros(n))
    initOption!(wk, WK_DXQ , zeros(n))
    initOption!(wk, WK_XA  , zeros(n))
    initOption!(wk, WK_XWA , zeros(n))
    initOption!(wk, WK_F   , zeros(n))
    initOption!(wk, WK_FA  , zeros(n))
    initOption!(wk, WK_ETA , zeros(n))
    initOption!(wk, WK_XW  , zeros(n))
    initOption!(wk, WK_FW  , zeros(n))
    initOption!(wk, WK_DXQA, zeros(n))

    initOption!(wk, WK_SUMXA0, 0.0)
    initOption!(wk, WK_SUMXA1, 0.0)
    initOption!(wk, WK_FCMON,  0.0)
    initOption!(wk, WK_FCA,    0.0)
    initOption!(wk, WK_FCKEEP, 0.0)
    initOption!(wk, WK_FCPRI,  0.0)
    initOption!(wk, WK_DMYCOR, 0.0)
    initOption!(wk, WK_SUMXS,  0.0)

    # Initialize statistics
    stats[STATS_NITER]  = 0
    stats[STATS_NCORR]  = 0
    stats[STATS_NFCN]   = 0
    stats[STATS_NJAC]   = 0
    stats[STATS_NREJR1] = 0
    stats[STATS_NEW]    = 0
    stats[STATS_ICONV]  = 0
    stats[STATS_CONV]   = 0.0
    stats[STATS_SUMX]   = 0.0
    stats[STATS_DLEVF]  = 0.0

    if qIniMon
        write(printIO,"INFO: ","N = $n\n")
        write(printIO,"INFO: ","Prescribed relative precision ",
        "$(opt.options[OPT_RTOL])\n")
        if jacGen == 1
            message = "a user function"
        elseif jacGen == 2
            message = "numerical differentation (without feedback strategy)"
        elseif jacGen == 3
            message = "numerical differentation (feedback strategy included)"
        end
        write(printIO,"INFO: ","The Jacobian is supplied by $message\n")
        if mStor == 0
            message = "full"
        elseif mStor == 1
            message = "banded"
        end
        write(printIO,"INFO: ","The Jacobian will be stored in $message mode\n")
        if mStor == 1
            write(printIO,"INFO: ","Lower bandwidth : $ml \t",
            "Upper bandwidth : $mu\n")
        end
        if opt.options[OPT_NOROWSCAL] == 1
            message = "inhibited"
        else
            message = "allowed"
        end
        write(printIO,"INFO: ",
        "Automatic row scaling of the jacobian is $message\n")
    end

    # Check for non linear option
    nonLin = getOption!(opt,OPT_NONLIN,3)
    initOption!(opt,OPT_BOUNDEDDAMP,0)

    if opt.options[OPT_BOUNDEDDAMP] == 0
        qBDamp = Int(nonLin == 4)
    elseif opt.options[OPT_BOUNDEDDAMP] == 1
        qBDamp = 1
    elseif opt.options[OPT_BOUNDEDDAMP] == 2
        qBDamp = 0
    end

    # Initialize bounded damping strategy restriction factor
    initOption!(opt,OPT_FCBAND,0.0)
    if qBDamp == 1
        if opt.options[OPT_FCBAND] < 1.0
            opt.options[OPT_FCBAND] = 10.0
        end
    end

    if qIniMon
        if qRank1 == 1
            message = "allowed"
        else
            message = "inhibited"
        end
        write(printIO,"INFO: ","Rank-1 updates are $message\n")
        if nonLin == 1
            message = "linear"
        elseif nonLin == 2
            message = "mildly nonlinear"
        elseif nonLin == 3
            message = "highly nonlinear"
        elseif nonLin == 4
            message = "extremely nonlinear"
        end
        write(printIO,"INFO: ","Problem is specified as being $message\n")
        if qBDamp == 1
            write(printIO,"INFO: ","Bounded damping strategy is active\n",
            "bounding factor is $(opt.options[OPT_FCBAND])\n")
        else
            write(printIO,"INFO: ","Bounded damping strategy is off\n")
        end
        if qOrdi == 1
            write(printIO,"INFO: ","Special mode: ",
            "Ordinary Newton iteration will be done\n")
        end
        if qSimpl == 1
            write(printIO,"INFO: ","Special mode: ",
            "Simplified Newton iteration will be done\n")
        end
    end

    # Maximum permitted number of iteration steps
    nItmax = getOption!(opt,OPT_NITMAX,50)
    if nItmax <= 0
        nItmax = 50
        opt.options[OPT_NITMAX] = nItmax
    end

    if qIniMon
        write(printIO,"INFO: ","Maximum permitted number of ",
        "iteration steps : $nItmax\n")
    end

    # Initial damping factor for highly nonlinear problems
    initOption!(opt,OPT_FCSTART,0.0)
    qFcStart = opt.options[OPT_FCSTART] > 0.0
    if !qFcStart
        opt.options[OPT_FCSTART] = 1.0e-2
        if nonLin == 4
            opt.options[OPT_FCSTART] = 1.0e-4
        end
    end

    # Minimal permitted damping factor
    initOption!(opt,OPT_FCMIN,0.0)
    if opt.options[OPT_FCMIN] <= 0.0
        opt.options[OPT_FCMIN] = 1.0e-4
        if nonLin == 4
            opt.options[OPT_FCMIN] = 1.0e-8
        end
    end
    fcMin = getOption(opt,OPT_FCMIN,0.0)

    # Rank1 decision parameter SIGMA
    initOption!(opt,OPT_SIGMA,0.0)
    if opt.options[OPT_SIGMA] < 1.0
        opt.options[OPT_SIGMA] = 3.0
    end
    if qRank1 == 0
        opt.options[OPT_SIGMA] = 10.0/fcMin
    end

    # Decision parameter about increasing too small predictor
    # to greater corrector value
    initOption!(opt,OPT_SIGMA2,0.0)
    if opt.options[OPT_SIGMA2] < 1.0
        opt.options[OPT_SIGMA2] = 10.0/fcMin
    end

    # Starting value of damping factor (fcMin <= fc <= 1.0)
    if nonLin <= 2 && !qFcStart
        # for linear or mildly nonlinear problems
        fc = 1.0
    else
        # for highly or extremely nonlinear problems
        fc = getOption!(opt,OPT_FCSTART)
    end

    # Simplified Newton iteration implies ordinary Newton iteration mode
    if getOption!(opt,OPT_QSIMPL,0) == 1
        fc = 1.0
    end

    # If ordinary Newton iteration, damping factor is always 1
    if getOption!(opt,OPT_QORDI,0) == 1
        fc = 1.0
    end

    setOption!(opt,OPT_FCSTART,fc)

    if printIt >= 2 && qSucc == 0
        write(printIO,"INFO: ","Internal parameters:",
        "\n Starting value for damping factor ",
        "OPT_FCSTART = $(opt.options[OPT_FCSTART])",
        "\n Minimum allowed damping factor OPT_FCMIN = $fcMin",
        "\n Rank-1 updates decision parameter ",
        "OPT_SIGMA = $(opt.options[OPT_SIGMA])\n")
    end

    retCode = -1

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    if nBroy == 0
        nBroy = 1
    end

    # TODO: Initialize statistics before passing to n1int

    # Call to n1int
    (x,xScal,retCode, stats[STATS_CONV], stats[STATS_DLEVF], stats[STATS_NCORR],
    stats[STATS_NFCN], stats[STATS_NJAC], stats[STATS_NREJR1],
    stats[STATS_NEW],stats[STATS_ICONV]) = n1int(n, fcn, jac, x, xScal,
    opt.options[OPT_RTOL], nItmax, nonLin, opt, retCode, wk, m1, m2, nBroy,
    xIter, sumXall, dLevFall, sumXQall, tolAll, fcAll, wk.options[WK_A],
    wk.options[WK_DXSAVE], wk.options[WK_DX], wk.options[WK_DXQ],
    wk.options[WK_XA],wk.options[WK_XWA],wk.options[WK_F],wk.options[WK_FA],
    wk.options[WK_ETA], wk.options[WK_XW], wk.options[WK_FW], wk.options[WK_DXQA],
    wk.options[WK_SUMXA0],wk.options[WK_SUMXA1],wk.options[WK_FCMON],
    opt.options[OPT_FCSTART], opt.options[OPT_FCMIN], opt.options[OPT_SIGMA],
    opt.options[OPT_SIGMA2], wk.options[WK_FCA], wk.options[WK_FCKEEP],
    wk.options[WK_FCPRI], wk.options[WK_DMYCOR], stats[STATS_CONV],
    wk.options[WK_SUMXS], stats[STATS_DLEVF], mStor, printWarning,
    printIterationMonitor, printSolution, printIO, stats[STATS_NITER],
    stats[STATS_NCORR], stats[STATS_NFCN], stats[STATS_NJAC], stats[STATS_NREJR1],
    stats[STATS_NEW],stats[STATS_ICONV], qBDamp, stats)

    # TODO: This is supposed to happen inside n1int
    push!(xIter,1)
    push!(sumXall,1)
    push!(sumXQall,1)
    push!(dLevFall,1)
    push!(tolAll,1)
    push!(fcAll,1)
    stats[STATS_NITER]  = 1
    stats[STATS_NCORR]  = 1
    stats[STATS_NREJR1] = 1
    stats[STATS_NFCN]   = 1
    stats[STATS_NJAC]   = 1

    # set stats variable
    stats[STATS_XSCAL] = xScal;
    if retCode == -1
        stats[STATS_RTOL] = tolAll[stats[STATS_NITER]]
    else
        stats[STATS_RTOL] = opt.options[OPT_RTOL]
    end
    stats[STATS_XITER]          = xIter
    stats[STATS_NATLEVEL]       = sumXall
    stats[STATS_SIMLEVEL]       = sumXQall
    stats[STATS_STDLEVEL]       = dLevFall
    stats[STATS_PRECISION]      = tolAll
    stats[STATS_DAMPINGFC]      = fcAll

    # Print statistics
    if printIt >= 2 && retCode != -1 && retCode != 10
        write(printIO,"\n",
        @sprintf("*************   Statistics   ************\n"),
        @sprintf("***  Newton-iterations     : %7i  ***\n", (stats[STATS_NITER])),
        @sprintf("***  Corrector steps       : %7i  ***\n", (stats[STATS_NCORR])),
        @sprintf("***  Rejected Rank-1 steps : %7i  ***\n", (stats[STATS_NREJR1])),
        @sprintf("***  Jacobian evaluations  : %7i  ***\n", (stats[STATS_NJAC])),
        @sprintf("***  Function evaluations  : %7i  ***\n", (stats[STATS_NFCN])),
        @sprintf("*****************************************\n"))
    end

    return (stats, retCode)
end
