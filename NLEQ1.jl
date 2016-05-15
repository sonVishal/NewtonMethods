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
#     Output parameters:
#         x:          Solution if found
#         xScal:      Lat internal scaling vector if successful solve
#         stats:      Statistics of the solution
#         retCode:    Return code signifying success or failure
# This function does the bookkeeping. The actual computation are done in n1int

include("CheckOptionsNLEQ1.jl")
include("Jacobian.jl")
function nleq1(fcn::Function,x::Vector,xScal::Vector,opt::OptionsNLEQ)

    # TODO: persistent variables
    # Might be a good idea to store them inside opt
    # so that they can be used again and again.

    # Initialize output statistics
    stats = Dict{ASCIIString,Any}()

    # Initialize error code 0
    retCode = 0

    # Print warning messages?
    printFlag = getOption!(opt,OPT_PRINTWARNING,0)

    # First call or successive call
    qsucc = getOption!(opt,OPT_QSUCC,0)

    # Check input parameters and options
    n = length(x)
    retCode = checkOptions(n,x,xScal,opt)

    # Exit if any parameter error was detected
    if retCode != 0
        stats["Error_Code"] = retCode
        error("Exit with return code $retCode")
    end

    # Check if this is a first call or successive call
    # to nleq1
    if qsucc != 0
        # TODO: check what this should be doing
    end

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    mstor = getOption!(opt,OPT_MSTOR,0)
    if mstor == 0
        m1 = n
        m2 = n
    elseif mstor == 1
        ml = getOption!(opt,"ML",0)
        mu = getOption!(opt,"MU",0)
        m1 = 2*ml + mu + 1
        m2 = ml + mu + 1
    end

    # Assign the Jacobian depending on user input
    # Multiple dispatch calls the required function based on
    # the storage requirement of the user
    jacFcn = getOption!(opt,OPT_JACOBIFCN,0)
    if jacFcn == 0 || jacFcn == 2
        jacFcn = numDiffJac
    elseif jacFcn == 3
        jacFcn = numDiffJacFB
    end
    opt.options[OPT_JACOBIFCN] = jacFcn;

    qrank1 = getOption!(opt,OPT_QRANK1,0)
    qordi  = getOption!(opt,OPT_QORDI,0)
    qsimpl = getOption!(opt,OPT_QSIMPL,0)

    if qrank1 == 1
        nbroy = getOption!(opt,OPT_NBROY,0)
        if nbroy == 0
            nbroy = max(m2,10)
            opt.options[OPT_NBROY] = nbroy
        end
    else
        nbroy = 0
    end

    # Workspace: WK
    # wk_A = zeros(m1,n)
    #
    # if qrank1 == 1
    #     wk_DXSAVE = zeros(n,nbroy)
    # else
    #     wk_DXSAVE = 0.0
    # end
    #
    # wk_DX       = zeros(n)
    # wk_DXQ      = zeros(n)
    # wk_XA       = zeros(n)
    # wk_XWA      = zeros(n)
    # wk_F        = zeros(n)
    # wk_FA       = zeros(n)
    # wk_ETA      = zeros(n)
    # wk_XW       = zeros(n)
    # wk_FW       = zeros(n)
    # wk_DXQA     = zeros(n)
    # wk_SUMXA0   = 0.0
    # wk_SUMXA1   = 0.0

    initOption!(opt,OPT_QNSCAL => 0)

    # TODO: Print log of things done till now

    # Line 742 starts here
    nonLin = getOption!(opt,OPT_NONLIN,3)
    initOption!(opt,OPT_IBDAMP,0)

    if opt.options[OPT_IBDAMP] == 0
        qbdamp = Int(nonLin == 4)
    elseif opt.options[OPT_IBDAMP] == 1
        qbdamp = 1
    elseif opt.options[OPT_IBDAMP] == 2
        qbdamp = 0
    end

    initOption!(opt,OPT_FCBAND,0.0)

    if qbdamp == 1
        if opt.options[OPT_FCBAND] < 1.0
            opt.options[OPT_FCBAND] = 10.0
        end
    end

    # TODO: print somethings

    # Maximum permitted number of iteration steps
    nItmax = getOption!(opt,OPT_NITMAX,50)
    if nItmax <= 0
        nItmax = 50
    end
    opt.options[OPT_NITMAX] = nItmax

    # TODO: Print somethings

    # Initial damping factor for highly nonlinear problems
    initOption!(opt,OPT_FCSTART,0.0)
    qfcstr = opt.options[OPT_FCSTART] > 0.0
    if !qfcstr
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
    if qrank1 == 0
        opt.options[OPT_SIGMA] = 10.0/fcmin
    end

    # Decision parameter about increasing too small predictor
    # to greater corrector value
    initOption!(opt,OPT_SIGMA2,0.0)
    if opt.options[OPT_SIGMA2] < 1.0
        opt.options[OPT_SIGMA2] = 10.0/fcmin
    end

    # Starting value of damping factor (fcmin <= fc <= 1.0)
    if nonLin <= 2 && !qfcstr
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
    if nbroy == 0
        nbroy = 1
    end

    # Call to n1int

    # Print statistics

    # set stats variable

    return (stats, retCode)
end
