# TODO: Make everything type independent
# Currently everything is assumed to be Float64
# TODO: Better comments.
# TODO: Logging and printing of Logs!!

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

    # Initialize error code 0
    retCode = 0

#-------------------------------------------------------------------------------
# Printing related stuff
#-------------------------------------------------------------------------------
    # Print warning messages?
    printFlag   = getOption!(opt,OPT_PRINTWARNING,0)
    # Print iteration summary?
    printItMon  = getOption!(opt,OPT_PRINTITERATIONMON,0)
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
    qIniMon = (printItMon >= 1 && qSucc == 0)

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
    if qSucc != 0
        # If this is the first call then assign memory to the variables
        xIter       = []
        sumXall     = []
        dLevFall    = []
        sumXQall    = []
        tolAll      = []
        fCall       = []
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
    jacMethod = getOption!(opt,OPT_JACMETHOD,0)
    if jacMethod == 0
        jacMethod = 2
        opt.options[OPT_JACMETHOD] = jacMethod;
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

    # Workspace: WK
    initOption!(wk,WK_A,0.0)

    if qRank1 == 1
        initOption!(wk,WK_DXSAVE,zeros(n,nBroy))
    else
        initOption!(wk,WK_DXSAVE,0.0)
    end

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
    initOption!(wk, WK_FCMIN,  0.0)
    initOption!(wk, WK_SIGMA,  0.0)
    initOption!(wk, WK_SIGMA2, 0.0)
    initOption!(wk, WK_FCA,    0.0)
    initOption!(wk, WK_FCKEEP, 0.0)
    initOption!(wk, WK_FCPRI,  0.0)
    initOption!(wk, WK_DMYCOR, 0.0)
    initOption!(wk, WK_CONV,   0.0)
    initOption!(wk, WK_SUMX,   0.0)
    initOption!(wk, WK_SUMXS,  0.0)
    initOption!(wk, WK_DLEVF,  0.0)
    initOption!(wk, WK_NITER,  0)
    initOption!(wk, WK_NCORR,  0)
    initOption!(wk, WK_NFCN,   0)
    initOption!(wk, WK_NJAC,   0)
    initOption!(wk, WK_NFCNJ,  0)
    initOption!(wk, WK_NREJR1, 0)
    initOption!(wk, WK_NEW,    0)
    initOption!(wk, WK_ICONV,  0)

    initOption!(opt,OPT_NOROWSCAL => 0)

    if qIniMon
        write(printIO,"\nN = $n\n")
        write(printIO,"\nPrescribed relative precision $(opt.options[OPT_RTOL])\n")
        if jacMethod == 1
            message = "a user function"
        elseif jacMethod == 2
            message = "numerical differentation (without feedback strategy)"
        elseif jacMethod == 3
            message = "numerical differentation (feedback strategy included)"
        end
        write(printIO,"\nThe Jacobian is supplied by $message\n")
        if mStor == 0
            message = "full"
        elseif mStor == 1
            message = "banded"
        end
        write(printIO,"The Jacobian will be stored in $message mode\n")
        if mStor == 1
            write(printIO,"Lower bandwidth : $ml \t Upper bandwidth : $mu\n")
        end
        if opt.options[OPT_NOROWSCAL] == 1
            message = "inhibited"
        else
            message = "allowed"
        end
        write(printIO,"Automatic row scaling of the jacobian is $message\n")
    end

    # Line 742 starts here
    nonLin = getOption!(opt,OPT_NONLIN,3)
    initOption!(opt,OPT_IBDAMP,0)

    if opt.options[OPT_IBDAMP] == 0
        qBDamp = Int(nonLin == 4)
    elseif opt.options[OPT_IBDAMP] == 1
        qBDamp = 1
    elseif opt.options[OPT_IBDAMP] == 2
        qBDamp = 0
    end

    initOption!(opt,OPT_FCBAND,0.0)

    if qBDamp == 1
        if opt.options[OPT_FCBAND] < 1.0
            opt.options[OPT_FCBAND] = 10.0
        end
    end

    # TODO: print somethings

    # Maximum permitted number of iteration steps
    nItmax = getOption!(opt,OPT_NITMAX,50)
    if nItmax <= 0
        nItmax = 50
        opt.options[OPT_NITMAX] = nItmax
    end

    # TODO: Print somethings

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
    fcmin = getOption!(opt,OPT_FCMIN,0.0)

    # Rank1 decision parameter SIGMA
    initOption!(opt,OPT_SIGMA,0.0)
    if opt.options[OPT_SIGMA] < 1.0
        opt.options[OPT_SIGMA] = 3.0
    end
    if qRank1 == 0
        opt.options[OPT_SIGMA] = 10.0/fcmin
    end

    # Decision parameter about increasing too small predictor
    # to greater corrector value
    initOption!(opt,OPT_SIGMA2,0.0)
    if opt.options[OPT_SIGMA2] < 1.0
        opt.options[OPT_SIGMA2] = 10.0/fcmin
    end

    # Starting value of damping factor (fcmin <= fc <= 1.0)
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

    opt.options[OPT_FCSTART] = fc

    # Print something

    # Time measurement

    retCode = -1

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    if nBroy == 0
        nBroy = 1
    end

    # Call to n1int
    # n1int(n,fcn,jacFcn,x,xScal,)

    # Print statistics

    # set stats variable

    return (stats, retCode)
end
