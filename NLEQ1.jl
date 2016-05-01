include("CheckOptionsNLEQ1.jl");

# Bookkeeping
# TODO: Make everything type independent
# Currently everything is assumed to be Float64
function nleq1(f,x,xscal,opt::OptionsNLEQ,wk::OptionsNLEQ)

    # TODO: persistent variables

    # First call or successive call
    qsucc = getOrSet(opt,"QSUCC",0);

    # Initialize output statistics
    stats = Dict{ASCIIString,Any}();

    # Initialize error code 0
    retCode = 0;

    # To print or not?
    printFlag = getOrSet(opt,"PR_ERR",0);

    # Check input parameters and options
    n = length(x);
    retCode = checkOptions(n,x,xscal,opt);

    # Exit if any parameter error was detected
    if retCode != 0
        stats["Error_Code"] = retCode;
        return (stats, retCode)
    end

    # Check if this is a first call or successive call
    # to nleq1
    if qsucc != 0
        # TODO: check what this should be doing
    end

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    mstor = getOrSet(opt,"MSTOR",0);
    if mstor == 0
        m1 = n;
        m2 = n;
    elseif mstor == 1
        ml = getOrSet(opt,"ML",0);
        mu = getOrSet(opt,"MU",0);
        m1 = 2*ml + mu + 1;
        m2 = ml + mu + 1;
    end

    # User given Jacobian or not
    jacgen = getOrSet(opt,"JACGEN",0);
    if jacgen == 0
        jacgen = 2;
    end
    opt.options["JACGEN"] = jacgen;

    qrank1 = getOrSet(opt,"QRANK1",0);
    qordi  = getOrSet(opt,"QORDI",0);
    qsimpl = getOrSet(opt,"QSIMPL",0);

    if qrank1 == 1
        nbroy = getOrSet(wk,"NBROY",0);
        if nbroy == 0
            nbroy = max(m2,10);
        end
        wk.options["NBROY"] = nbroy;
    else
        nbroy = 0;
    end

    # Workspace: WK

    initOption!(wk,"A",0.0);

    if qrank1 == 1
        wk.options["DXSAVE"] = zeros(n,nbroy);
    else
        wk.options["DXSAVE"] = 0.0;
    end

    initOptions!(wk,"DX" => zeros(n),"DXQ" => zeros(n),"XA" => zeros(n));
    initOptions!(wk,"XWA" => zeros(n),"F" => zeros(n),"FA" => zeros(n));
    initOptions!(wk,"ETA" => zeros(n),"XW" => zeros(n),"FW" => zeros(n));
    initOption!(wk,"DXQA",zeros(n));

    initOptions!(wk,"SUMXA0" => 0.0, "SUMXA1" => 0.0, "FCMON" => 0.0);
    initOptions!(wk,"FCMIN" => 0.0, "SIGMA" => 0.0, "SIGMA2" => 0.0);
    initOptions!(wk,"FCA" => 0.0, "FCKEEP" => 0.0, "FCPRI" => 0.0);
    initOptions!(wk,"DMYCOR" => 0.0, "CONV" => 0.0, "SUMX" => 0.0);
    initOptions!(wk,"DLEVF" => 0.0, "NITER" => 0, "NCORR" => 0);
    initOptions!(wk,"NFCN" => 0, "NJAC" => 0, "NFCNJ" => 0);
    initOptions!(wk,"NREJR1" => 0, "NEW" => 0, "ICONV" => 0);

    initOption!(opt,"NOROWSCAL" => 0);

    # TODO: Print log of things done till now

    # Line 742 starts here
    nonlin = getOrSet(opt,"NONLIN",3);
    initOption!(opt,"BOUNDEDDAMP",0);

    if opt.options["BOUNDEDDAMP"] == 0
        qbdamp = Int(nonlin == 4);
    elseif opt.options["BOUNDEDDAMP"] == 1
        qbdamp = 1;
    elseif opt.options["BOUNDEDDAMP"] == 2
        qbdamp = 0;
    end

    initOption!(wk,"FCBND",0.0);

    if qbdamp == 1
        if wk.options["FCBND"] < 1.0
            wk.options["FCBND"] = 10.0;
        end
    end

    # TODO: print somethings

    # Maximum permitted number of iteration steps
    nitmax = getOrSet(wk,"NITMAX",50);
    if nitmax <= 0
        nitmax = 50;
    end
    wk.options["NITMAX"] = nitmax;

    # TODO: Print somethings

    # Initial damping factor for highly nonlinear problems
    initOption!(wk,"FCSTART",0.0);
    qfcstr = wk.options["FCSTART"] > 0.0;
    if !qfcstr
        wk.options["FCSTART"] = 1.0e-2;
        if nonlin == 4
            wk.options["FCSTART"] = 1.0e-4;
        end
    end

    # Minimal permitted damping factor
    initOption!(wk,"FCMIN",0.0);
    if wk.options["FCMIN"] <= 0.0
        wk.options["FCMIN"] = 1.0e-4;
        if nonlin == 4
            wk.options["FCMIN"] = 1.0e-8;
        end
    end
    fcmin = getOrSet(wk,"FCMIN",0.0);

    # Rank1 decision parameter SIGMA
    initOption!(wk,"SIGMA",0.0);
    if wk.options["SIGMA"] < 1.0
        wk.options["SIGMA"] = 3.0;
    end
    if qrank1 == 0
        wk.options["SIGMA"] = 10.0/fcmin;
    end

    # Decision parameter about increasing too small predictor
    # to greater corrector value
    initOption!(wk,"SIGMA2",0.0);
    if wk.options["SIGMA2"] < 1.0
        wk.options["SIGMA2"] = 10.0/fcmin;
    end

    # Starting value of damping factor (fcmin <= fc <= 1.0)
    if nonlin <= 2 && !qfcstr
        # for linear or mildly nonlinear problems
        fc = 1.0;
    else
        # for highly or extremely nonlinear problems
        fc = getOrSet(wk,"FCSTART");
    end

    # Simplified Newton iteration implies ordinary Newton iteration mode
    if getOrSet(opt,"QSIMPL",0) == 1
        fc = 1.0;
    end

    # If ordinary Newton iteration, damping factor is always 1
    if getOrSet(opt,"QORDI",0) == 1
        fc = 1.0;
    end

    wk.options["FCSTART"] = fc

    # Print something

    # Time measurement

    retCode = -1;

    # If retCode is unmodified on exit, successive steps are required
    # to complete the Newton iterations
    if nbroy == 0
        nbroy = 1;
    end

    # Call to n1int

    # Print statistics

    # set stats variable

    return (stats, retCode)
end
