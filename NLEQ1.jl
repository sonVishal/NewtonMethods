include("CheckOptionsNLEQ1.jl");

# Bookkeeping
# TODO: Make everything type independent
# Currently everything is assumed to be Float64
function nleq1(f,x,xscal,opt::OptionsNLEQ,wk::OptionsNLEQ)
    # Initialize output statistics
    stats = Dict{ASCIIString,Any}();

    # Initialize the return code for the solver
    retcode = 0;

    # Initialize error code 0
    retCode = 0;

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
    qsucc = getOption(opt,"QSUCC",0);
    if qsucc != 0
        # TODO: check what this should be doing
    end

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    mstor = getOption(opt,"MSTOR",0);
    if mstor == 0
        m1 = n;
        m2 = n;
    elseif mstor ==1
        ml = getOption(opt,"ML",0);
        mu = getOption(opt,"MU",0);
        m1 = 2*ml + mu + 1;
        m2 = ml + mu + 1;
    end

    # User given Jacobian or not
    jacgen = getOption(opt,"JACGEN",0);
    if jacgen == 0
        jacgen = 2;
    end
    setOption!(opt,"JACGEN",jacgen);

    qrank1 = getOption(opt,"QRANK1",0);
    qordi  = getOption(opt,"QORDI",0);
    qsimpl = getOption(opt,"QSIMPL",0);

    if qrank1 == 1
        nbroy = getOption(wk,"NBROY",0);
        if nbroy == 0
            nbroy = max(m2,10);
        end
        setOption!(wk,"NBROY",nbroy);
    else
        nbroy = 0;
    end

    # Workspace: WK

    initOption!(wk,"A",0.0);

    if qrank1 == 1
        setOption!(wk,"DXSAVE",zeros(n,nbroy));
    else
        setOption!(wk,"DXSAVE",0.0);
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

    initOption!(opt,"NOROWSCAL",0);

    # TODO: Print log of things done till now

    # Maximum permitted number of iteration steps
    nitmax = getOption(wk,"NITMAX",50);
    if nitmax <= 0
        nitmax = 50;
    end
    setOption!(wk,"NITMAX",nitmax);

    # TODO: Print something

    return (stats, retCode)
end
