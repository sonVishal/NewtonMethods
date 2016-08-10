__precompile__(false)
module NewtonMethods
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

include("Error.jl")
include("Helper.jl")
include("Jacobian.jl")
include("Options.jl")
include("Constants.jl")
include("CheckOptionsNLEQ1.jl")
include("Error.jl")
include("NLEQ1Main.jl")
function nleq1(fcn::Function,x::Vector,xScal::Vector,opt::OptionsNLEQ,wk::OptionsNLEQ)

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
    printMon     = getOption!(opt,OPT_PRINTITERATION,0)
    # Print solution summary?
    printSol    = getOption!(opt,OPT_PRINTSOLUTION,0)
    # Where to print?
    # Defaults to STDOUT
    printIOwarn = getOption(opt,OPT_PRINTIOWARN,STDOUT)
    printIOmon  = getOption!(opt,OPT_PRINTIOMON,STDOUT)
    printIOsol  = getOption!(opt,OPT_PRINTIOSOL,STDOUT)
    # if printIO == "FILE"
    #     # If not STDOUT then print to file
    #     # Default file name is log.txt and the file is opened for writing
    #     printFileName   = getOption!(opt,OPT_PRINTFILENAME,"log.txt")
    #     printFileMode   = getOption!(opt,OPT_PRINTFILEMODE,"w")
    #     if printFileMode != "w" || printFileMode != "a"
    #         throw(InvalidOption("OPT_PRINTFILEMODE",printFileMode))
    #     end
    #     f = open(printFileName,printFileMode)
    # end
#-------------------------------------------------------------------------------

    # First call or successive call
    qSucc   = Bool(getOption!(opt,OPT_QSUCC,0))
    qIniMon = (printMon >= 1 && !qSucc)

    # Check input parameters and options
    n = length(x)
    retCode = checkOptions(n,x,xScal,opt)

    # Exit if any parameter error was detected
    if retCode != 0
        error("Exit with return code $retCode")
    end

    # Check if this is a first call or successive call
    # to nleq1
    # if !qSucc
    # If this is the first call then assign memory to the variables
    xIter    = getOption!(wk,"persistent_xIter",[])
    sumXall  = getOption!(wk,"persistent_sumXall",[])
    dLevFall = getOption!(wk,"persistent_dLevFall",[])
    sumXQall = getOption!(wk,"persistent_sumXQall",[])
    tolAll   = getOption!(wk,"persistent_tolAll",[])
    fcAll    = getOption!(wk,"persistent_fcAll",[])
    # end

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
    end
    setOption!(opt, OPT_JACGEN, jacGen)

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

    # Initialize options
    initOption!(opt, OPT_FCMIN,     0.0)
    initOption!(opt, OPT_SIGMA,     0.0)
    initOption!(opt, OPT_SIGMA2,    0.0)
    initOption!(opt, OPT_NOROWSCAL, 0)

    # Workspace: WK
    initOption!(wk, WK_A, zeros(m1,n))

    if qRank1
        initOption!(wk, WK_DXSAVE, zeros(n,nBroy))
    else
        initOption!(wk, WK_DXSAVE, 0.0)
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

    initOption!(wk, STATS_NITER,  0)
    initOption!(wk, STATS_NCORR,  0)
    initOption!(wk, STATS_NFCN,   0)
    initOption!(wk, STATS_NFCNJ,  0)
    initOption!(wk, STATS_NJAC,   0)
    initOption!(wk, STATS_NREJR1, 0)
    initOption!(wk, STATS_NEW,    0)
    initOption!(wk, STATS_ICONV,  0)
    initOption!(wk, STATS_CONV,   0.0)
    initOption!(wk, STATS_SUMX,   0.0)
    initOption!(wk, STATS_DLEVF,  0.0)
    initOption!(wk, STATS_RTOL,   0.0)

    if qIniMon
        write(printIOmon,"\nINFO: ","N = $n\n")
        write(printIOmon,"\nINFO: ","Prescribed relative precision ",
        "$(opt.options[OPT_RTOL])\n")
        if jacGen == 1
            message = "a user function"
        elseif jacGen == 2
            message = "numerical differentation (without feedback strategy)"
        elseif jacGen == 3
            message = "numerical differentation (feedback strategy included)"
        end
        write(printIOmon,"\nINFO: ","The Jacobian is supplied by $message\n")
        if mStor == 0
            message = "full"
        elseif mStor == 1
            message = "banded"
        end
        write(printIOmon,"INFO: ","The Jacobian will be stored in $message mode\n")
        if mStor == 1
            write(printIOmon,"INFO: ","Lower bandwidth : $ml \t",
            "Upper bandwidth : $mu\n")
        end
        if opt.options[OPT_NOROWSCAL] == 1
            message = "inhibited"
        else
            message = "allowed"
        end
        write(printIOmon,"INFO: ",
        "Automatic row scaling of the jacobian is $message\n")
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

    if qIniMon
        if qRank1
            message = "allowed"
        else
            message = "inhibited"
        end
        write(printIOmon,"\nINFO: ","Rank-1 updates are $message\n")
        if nonLin == 1
            message = "linear"
        elseif nonLin == 2
            message = "mildly nonlinear"
        elseif nonLin == 3
            message = "highly nonlinear"
        elseif nonLin == 4
            message = "extremely nonlinear"
        end
        write(printIOmon,"INFO: ","Problem is specified as being $message\n")
        if qBDamp
            write(printIOmon,"INFO: ","Bounded damping strategy is active\n",
            "bounding factor is $(opt.options[OPT_FCBAND])\n")
        else
            write(printIOmon,"INFO: ","Bounded damping strategy is off\n")
        end
        if qOrdi
            write(printIOmon,"INFO: ","Special mode: ",
            "Ordinary Newton iteration will be done\n")
        end
        if qSimpl
            write(printIOmon,"INFO: ","Special mode: ",
            "Simplified Newton iteration will be done\n")
        end
    end

    # Maximum permitted number of iteration steps
    nItmax = getOption!(opt, OPT_NITMAX, 50)
    if nItmax <= 0
        nItmax = 50
        setOption!(opt, OPT_NITMAX, nItmax)
    end

    if qIniMon
        write(printIOmon,"INFO: ","Maximum permitted number of ",
        "iteration steps : $nItmax\n")
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

    setOption!(opt, OPT_FCSTART, fc)

    if printMon >= 2 && !qSucc
        write(printIOmon,"\nINFO: ","Internal parameters:",
        "\n\tStarting value for damping factor ",
        @sprintf("OPT_FCSTART\t= %1.2e",opt.options[OPT_FCSTART]),
        @sprintf("\n\tMinimum allowed damping factor OPT_FCMIN\t= %1.2e",fcMin),
        "\n\tRank-1 updates decision parameter ",
        @sprintf("OPT_SIGMA\t= %1.2e\n",opt.options[OPT_SIGMA]))
    end

    retCode = -1

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    if nBroy == 0
        nBroy = 1
    end

    # Call to n1int
    (x, xScal, retCode) = n1int(n, fcn, x, xScal,
    opt.options[OPT_RTOL], nItmax, nonLin, opt, retCode, wk, m1, m2, nBroy,
    xIter, sumXall, dLevFall, sumXQall, tolAll, fcAll, wk.options[WK_A],
    wk.options[WK_DXSAVE], wk.options[WK_DX], wk.options[WK_DXQ],
    wk.options[WK_XA],wk.options[WK_XWA],wk.options[WK_F],wk.options[WK_FA],
    wk.options[WK_ETA], wk.options[WK_XW], wk.options[WK_FW], wk.options[WK_DXQA],
    wk.options[WK_SUMXA0],wk.options[WK_SUMXA1],wk.options[WK_FCMON],
    opt.options[OPT_FCSTART], opt.options[OPT_FCMIN], opt.options[OPT_SIGMA],
    opt.options[OPT_SIGMA2], wk.options[WK_FCA], wk.options[WK_FCKEEP],
    wk.options[WK_FCPRI], wk.options[WK_DMYCOR], wk.options[STATS_CONV],
    wk.options[WK_SUMXS], wk.options[STATS_DLEVF], mStor, printWarn,
    printMon, printSol, printIOwarn, printIOmon, printIOsol, wk.options[STATS_NITER],
    wk.options[STATS_NCORR], wk.options[STATS_NFCN], wk.options[STATS_NFCNJ], wk.options[STATS_NJAC], wk.options[STATS_NREJR1],
    wk.options[STATS_NEW],wk.options[STATS_ICONV], qBDamp)

    # set stats variable
    stats = Dict{ASCIIString,Any}()
    stats[STATS_XSCAL] = xScal
    if retCode == -1
        stats[STATS_RTOL] = tolAll[wk.options[STATS_NITER]]
    else
        stats[STATS_RTOL] = opt.options[OPT_RTOL]
    end
    stats[STATS_XITER]      = xIter
    stats[STATS_NATLEVEL]   = sumXall
    stats[STATS_SIMLEVEL]   = sumXQall
    stats[STATS_STDLEVEL]   = dLevFall
    stats[STATS_PRECISION]  = tolAll
    stats[STATS_DAMPINGFC]  = fcAll
    stats[STATS_NITER]      = wk.options[STATS_NITER]
    stats[STATS_NCORR]      = wk.options[STATS_NCORR]
    stats[STATS_NREJR1]     = wk.options[STATS_NREJR1]
    stats[STATS_NJAC]       = wk.options[STATS_NJAC]
    stats[STATS_NFCN]       = wk.options[STATS_NFCN]
    stats[STATS_NFCNJ]      = wk.options[STATS_NFCNJ]
    # Print statistics
    if printMon >= 2 && retCode != -1 && retCode != 10
        write(printIOmon,"\n",
        @sprintf("*************   Statistics   ************\n"),
        @sprintf("***  Newton-iterations     : %7i  ***\n", (stats[STATS_NITER])),
        @sprintf("***  Corrector steps       : %7i  ***\n", (stats[STATS_NCORR])),
        @sprintf("***  Rejected Rank-1 steps : %7i  ***\n", (stats[STATS_NREJR1])),
        @sprintf("***  Jacobian evaluations  : %7i  ***\n", (stats[STATS_NJAC])),
        @sprintf("***  Function evaluations  : %7i  ***\n", (stats[STATS_NFCN])),
        @sprintf("***  ... for Jacobain eval : %7i  ***\n", (stats[STATS_NFCNJ])),
        @sprintf("*****************************************\n"))
    end

    # Assign the persistent variables back
    setOption!(wk, "persistent_xIter", xIter)
    setOption!(wk, "persistent_sumXall", sumXall)
    setOption!(wk, "persistent_dLevFall", dLevFall)
    setOption!(wk, "persistent_sumXQall", sumXQall)
    setOption!(wk, "persistent_tolAll", tolAll)
    setOption!(wk, "persistent_fcAll", fcAll)

    return (x, stats, retCode);
end
export nleq1,OptionsNLEQ;
export OPT_RTOL,
OPT_QSUCC,
OPT_MODE,
OPT_JACGEN,
OPT_JACFCN,
OPT_MSTOR,
OPT_ML,
OPT_MU,
OPT_ISCAL,
OPT_PRINTWARNING,
OPT_PRINTITERATION,
OPT_PRINTIOWARN,
OPT_PRINTIOMON,
OPT_PRINTIOSOL,
OPT_PRINTSOLUTION,
OPT_NONLIN,
OPT_QRANK1,
OPT_QORDI,
OPT_QSIMPL,
OPT_NOROWSCAL,
OPT_BOUNDEDDAMP,
OPT_IORMON,
OPT_NITMAX,
OPT_FCBAND,
OPT_SIGMA,
OPT_SIGMA2,
OPT_AJDEL,
OPT_AJMIN,
OPT_ETADIF,
OPT_ETAINI,
OPT_NBROY,
OPT_FCSTART,
OPT_FCMIN;
end
