"""
function checkOptions(n::Int64, x::Vector{Float64}, xScal::Vector{Float64},
    opt::OptionsNLEQ)

Checking of common input parameters and options.

## Input parameters
-------------------
| Variable | Description             |
|----------|-------------------------|
| n        | Size of the problem     |
| x        | Initial guess           |
| xScal*   | Initial scaling vector  |
| opt*     | Options set by the user |

(* marks inout parameters)

## Output parameters
--------------------
| Variable | Description                 |
|----------|-----------------------------|
| retCode  | Exit code in case of errors |
"""
function checkOptions(n::Int64, x::Vector{Float64}, xScal::Vector{Float64},
    opt::OptionsNLEQ)
    # Begin
    # Check the IO
    pIOwarn = getOption!(opt, OPT_PRINTIOWARN, STDOUT)
    if typeof(pIOwarn) != IOStream && pIOwarn != STDOUT
        retCode = 30
        write(STDOUT, "ERROR: Please provide a file stream for writing warnings or\n",
        "use the default option.\n")
        return retCode
    end

    pIOmon = getOption!(opt, OPT_PRINTIOMON, STDOUT)
    if typeof(pIOmon) != IOStream && pIOmon != STDOUT
        retCode = 30
        write(STDOUT, "ERROR: Please provide a file stream for writing iterations or\n",
        "use the default option.\n")
        return retCode
    end

    pIOsol = getOption!(opt, OPT_PRINTIOSOL, STDOUT)
    if typeof(pIOsol) != IOStream && pIOsol != STDOUT
        retCode = 30
        write(STDOUT, "ERROR: Please provide a file stream for writing solution or\n",
        "use the default option.\n")
        return retCode
    end
    # Check whether warnings need to be printed
    printWarn = getOption!(opt,OPT_PRINTWARNING,0)

    # Initialize return code to 0
    retCode = 0

    # Check dimensional parameter n
    if n <= 0
        retCode = 20
        write(pIOwarn,"ERROR: Bad input to dimensional parameter n supplied","\n",
            "Choose n positive, your input is: n = $n\n")
        return retCode
    end

    # Problem type specification by user
    nonLin = getOption!(opt,OPT_NONLIN,3)

    # Checking and conditional adaptation of user given RTOL
    # if RTOL is not set, set it to 1e-6
    rTol = getOption!(opt,OPT_RTOL,1e-6)
    if rTol <= 0.0
        retCode = 21
        write(pIOwarn,"ERROR: Nonpositive $OPT_RTOL supplied\n")
        return retCode
    else
        tolMin = epMach*10.0*n
        if rTol < tolMin
            rTol = tolMin
            setOption!(opt,OPT_RTOL,rTol)
            if printWarn == 1
                write(pIOwarn,"WARNING: User prescribed $OPT_RTOL increased to a reasonable smallest value RTOL = $rTol\n")
            end
        end

        tolMax = 1.0e-1
        if rTol > tolMax
            rTol = tolMax
            setOption!(opt,OPT_RTOL,rTol)
            if printWarn == 1
                write(pIOwarn,"WARNING: User prescribed $OPT_RTOL decreased to a reasonable largest value RTOL = $rTol\n")
            end
        end
    end

    # Test user prescribed accuracy and scaling on proper values
    if nonLin >= 3
        defScal = rTol
    else
        defScal = 1.0
    end

    for i = 1:n
        # Scaling Values cannot be negative
        # Positive scaling values give scale invariance
        if xScal[i] < 0.0
            retCode = 22
            write(pIOwarn,"ERROR: Negative value in xScal[$i] supplied\n")
            return retCode
        end

        if xScal[i] == 0.0
            xScal[i] = defScal
        end
        # Avoid overflow due to division by xScal[i]
        if xScal[i] > 0.0 && xScal[i] < small
            if printWarn == 1
                write(pIOwarn,"WARNING: xScal[$i] = $xScal[i] too small, increased to $small\n")
            end
            xScal[i] = small
        end
        # Avoid underflow due to division by xScal[i]
        if xScal[i] > great
            if printWarn == 1
                write(pIOwarn,"WARNING: xScal[$i] = $xScal[i] too big, increased to $great\n")
            end
            xScal[i] = great
        end
    end

    # Assign the Jacobian depending on user input
    # By default Forward mode automatic differentiation is used
    jacGen = getOption!(opt,OPT_JACGEN,4)
    if jacGen == 1
        jacFcn = getOption!(opt,OPT_JACFCN,0)
        if jacFcn == 0
            retCode = 30
            write(pIOwarn,"ERROR: The Jacobian function OPT_JACFCN is not supplied. ",
            "Please supply a Jacobian function or use OPT_JACGEN = 2 or 3 for numerical differentiation based jacobian evaluation.\n")
            return retCode
        end
    end

    return retCode
end

"""
function initializeOptions(opt::OptionsNLEQ, wk::OptionsNLEQ, n::Int64,
    m1::Int64, nBroy::Int64, qRank1::Bool, solver::Int64)

Initialization of options based on the solver input argument.

## Input parameters
-------------------
| Variables  | Description                                                    |
|------------|----------------------------------------------------------------|
| opt*       | Options set by the user                                        |
| wk*        | Internal workspace specific to the solver                      |
| n          | Size of the problem                                            |
| m1         | In full mode = n and in band mode = 2*ml+mu+1                  |
| nBroy      | Maximum number of possible consecutive iterative Broyden steps |
| qRank1     | Decision parameter for Rank-1 updates                          |
| solver = 1 | Specifies that the solver is NLEQ1                             |
|        = 2 | Specifies that the solver is NLEQ2                             |

(* marks inout parameters)
"""
function initializeOptions(opt::OptionsNLEQ, wk::OptionsNLEQ, n::Int64,
    m1::Int64, nBroy::Int64, qRank1::Bool, solver::Int64)
    # Begin
    # Initialize options: OPT
    initOption!(opt, OPT_FCMIN,     0.0)
    initOption!(opt, OPT_SIGMA,     0.0)
    initOption!(opt, OPT_SIGMA2,    0.0)
    initOption!(opt, OPT_NOROWSCAL, 0)

    # Initialize workspace: WK
    if solver == 1
        initOption!(wk, WK_A, zeros(m1,n))
        if qRank1
            initOption!(wk, WK_DXSAVE, zeros(n,nBroy))
        else
            initOption!(wk, WK_DXSAVE, 0.0)
        end
    elseif solver == 2
        initOption!(wk, WK_QU, zeros(n))
        initOption!(wk, WK_A, zeros(n,n))
        if nBroy != 0
            initOption!(wk, WK_QA_DXSAVE, zeros(n,n))
        end
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
    initOption!(wk, "P_TOLALL", Vector{Float64}())

    # Check whether these iteration variables should be stored or not
    initOption!(opt, OPT_STORE, 0)
    if opt.options[OPT_STORE] == 1
        initOption!(wk, "P_XITER", Vector{Vector{Float64}}())
        initOption!(wk, "P_SUMXALL", Vector{Float64}())
        initOption!(wk, "P_DLEVFALL", Vector{Float64}())
        initOption!(wk, "P_SUMXQALL", Vector{Float64}())
        initOption!(wk, "P_FCALL", Vector{Float64}())
    end

    return nothing
end

"""
function initializeOptions(opt::OptionsNLEQ, solver::Int64)

Initialization of options based on the solver input argument.

## Input parameters
-------------------
| Variables  | Description                                                    |
|------------|----------------------------------------------------------------|
| opt*       | Options set by the user                                        |
| solver = 1 | Specifies that the solver is NLEQ1                             |
|        = 2 | Specifies that the solver is NLEQ2                             |

(* marks inout parameters)
"""
function initializeOptions(n::Int64, opt::OptionsNLEQ, solver::Int64)
    # First call or successive call
    qSucc   = Bool(getOption!(opt,OPT_QSUCC,0))
    qIniMon = (opt.options[OPT_PRINTITERATION] >= 1 && !qSucc)

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

    # If in stepwise mode and this is the first call then clear the workspace
    initOption!(opt, OPT_MODE, 0)

    if opt.options[OPT_MODE] == 1 && !qSucc
        if solver == 1
            empty!(wkNLEQ1.options)
        elseif solver == 2
            empty!(wkNLEQ1.options)
        end
    end

    # Check if this is a first call or successive call to nleq1
    # If first call then reset the workspace and persistent variables
    if !qSucc
        if solver == 1
            empty!(wkNLEQ1.options)
            initializeOptions(opt, wkNLEQ1, n, m1, nBroy, qRank1, 1)
        elseif solver == 2
            empty!(wkNLEQ2.options)
            initializeOptions(opt, wkNLEQ2, n, m1, nBroy, qRank1, 2)
        end
    end

    # Check for non linear option
    nonLin = opt.options[OPT_NONLIN]
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
        if solver == 1
            printInitialization(n, opt.options[OPT_PRINTIOMON], opt.options[OPT_RTOL],
            opt.options[OPT_JACGEN], mStor, ml, mu, opt.options[OPT_NOROWSCAL],
            qRank1, nonLin, qBDamp, opt.options[OPT_FCBAND], qOrdi, qSimpl, nItmax)
        elseif solver == 2
            printInitialization(n, opt.options[OPT_PRINTIOMON], opt.options[OPT_RTOL],
            opt.options[OPT_JACGEN], 0, 0, 0, opt.options[OPT_NOROWSCAL], qRank1,
            nonLin, qBDamp, opt.options[OPT_FCBAND], false, false, nItmax)
        end
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

    if solver == 2
        iRank = getOption!(opt, OPT_IRANK, 0)
        if iRank <= 0 || iRank > n
            iRank = n
            setOption!(opt, OPT_IRANK, iRank)
        end

        cond = getOption!(opt, OPT_COND, 1.0/epMach)
        if cond < 1.0
            cond = 1.0/epMach
            setOption!(opt, OPT_COND, cond)
        end
    end

    if opt.options[OPT_PRINTITERATION] >= 2 && !qSucc
        write(opt.options[OPT_PRINTIOMON],"\nINFO: ","Internal parameters:",
        "\n\tStarting value for damping factor ",
        @sprintf("OPT_FCSTART\t= %1.2e",opt.options[OPT_FCSTART]),
        @sprintf("\n\tMinimum allowed damping factor OPT_FCMIN\t= %1.2e",fcMin),
        "\n\tRank-1 updates decision parameter ",
        @sprintf("OPT_SIGMA\t= %1.2e\n",opt.options[OPT_SIGMA]))
        if solver == 2
            write(opt.options[OPT_PRINTIOMON], @sprintf("\n\tInitial Jacobian pseudo-rank iRank\t\t= %6i", iRank),
            @sprintf("\n\tMaximum permitted subcondition cond\t\t= %1.2e\n", cond))
        end
    end

    return (m1, m2, nBroy, qBDamp)
end

"""
function printInitialization(n::Int64, printIOmon, rTol::Float64, jacGen::Int64,
    mStor::Int64, ml::Int64, mu::Int64, qNoRowScal::Int64, qRank1::Bool, nonLin::Int64,
    qBDamp::Bool, fcBand::Float64, qOrdi::Bool, qSimpl::Bool, nItmax::Int64)

Print a summary of the initialization.

## Input parameters
-------------------
| Variables  | Description                                        |
|------------|----------------------------------------------------|
| n          | Size of the problem                                |
| printIOmon | IO handle for printing                             |
| rTol       | Relative tolerance                                 |
| jacGen     | Method of Jacobian generation                      |
| mStor      | Dense or band mode storage of Jacobian             |
| ml         | Lower bandwidth in case of band storage            |
| mu         | Upper bandwidth in case of band storage            |
| qNoRowScal | Decision parameter for automatic row scaling       |
| qRank1     | Decision parameter for Rank-1 updates              |
| nonLin     | Problem type specification                         |
| qBDamp     | Decision parameter for bounded damping strategy    |
| fcBand     | Bounded damping strategy restriction factor        |
| qOrdi      | Decision parameter for ordinary Newton iteration   |
| qSimpl     | Decision parameter for simplified Newton iteration |
| nItmax     | Maximum permitted Newton iterations                |
"""
function printInitialization(n::Int64, printIOmon, rTol::Float64, jacGen::Int64,
    mStor::Int64, ml::Int64, mu::Int64, qNoRowScal::Int64, qRank1::Bool, nonLin::Int64,
    qBDamp::Bool, fcBand::Float64, qOrdi::Bool, qSimpl::Bool, nItmax::Int64)
    # Begin
    message = ""
    write(printIOmon,"\nINFO: ","N = $n\n")
    write(printIOmon,"\nINFO: ","Prescribed relative precision ",
    "$rTol\n")
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
    if qNoRowScal == 1
        message = "inhibited"
    else
        message = "allowed"
    end
    write(printIOmon,"INFO: ",
    "Automatic row scaling of the jacobian is $message\n")

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
        "bounding factor is $fcBand\n")
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

    write(printIOmon,"INFO: ","Maximum permitted number of ",
    "iteration steps : $nItmax\n")
end

"""
function printStats(stats::Dict{AbstractString,Any}, printIOmon)

Print a summary of the statistics.

## Input parameters
-------------------
| Variables  | Description                                      |
|------------|--------------------------------------------------|
| stats      | Dictionary variable containing solver statistics |
| printIOmon | IO handle for printing                           |
"""
function printStats(stats::Dict{AbstractString,Any}, printIOmon)
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

"""
function nScal(n::Int64, x::Vector{Float64}, xa::Vector{Float64}, xScal::Vector{Float64},
    iScal::Int64, mPr::Int64, printIO, xw::Vector{Float64})

To be used in connection with NLEQ1 and NLEQ2.
Computation of the internal scaling vector XW used for the
Jacobian matrix, the iterate vector and it's related
vectors - especially for the solution of the linear system
and the computations of norms to avoid numerical overflow.

## Input parameters
-------------------
| Variables | Description                                                                |
|-----------|----------------------------------------------------------------------------|
| n         | Size of the problem                                                        |
| x         | Current iterate                                                            |
| xa        | Previous iterate                                                           |
| xScal     | Scaling vector                                                             |
| iScal     | Decision parameter for scaling                                             |
| mPr       | Decision parameter for printing                                            |
| printIO   | IO handle for printing                                                     |

## Output parameters
-------------------
| Variables | Description                                                                |
|-----------|----------------------------------------------------------------------------|
| xw        | Scaling vector computed by this routine<br>All components must be positive.|

"""
function nScal(n::Int64, x::Vector{Float64}, xa::Vector{Float64}, xScal::Vector{Float64},
    iScal::Int64, mPr::Int64, printIO, xw::Vector{Float64})
    # Begin
    if iScal == 1
        xw[:] = xScal
    else
        for l1 = 1:n
            xw[l1] = max(xScal[l1],max((abs(x[l1])+abs(xa[l1]))*0.5,small))
        end
    end

    if mPr >= 6
        write(printIO,"\n\n",
        "+++++++++++++++++++++++++++++++++++++++++++++++++\n",
        "      x-components         scaling-components\n")
        for l1 = 1:n
            write(printIO,"  %18.10e   %18.10e\n",x[l1],xw[l1])
        end
        write(printIO,"+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
    end
    return nothing
end

"""
function nScrf(m::Int64, n::Int64, a::Array{Float64,2}, fw::Vector{Float64})

Row scaling of a (m,n)-matrix in full storage mode

## Input parameters
-------------------
| Variables | Description                    |
|-----------|--------------------------------|
| m         | Number of rows of the matrix   |
| n         | Numer of columns of the matrix |
| a[m,n]*   | Matrix to be scaled            |

(* marks inout parameters)

## Output parameters
-------------------
| Variables | Description                    |
|-----------|--------------------------------|
| fw        | Row scaling factors.           |
"""
function nScrf(m::Int64, n::Int64, a::Array{Float64,2}, fw::Vector{Float64})
    # Begin
    if issparse(a)
        nza = nnz(a)
        (row,col) = findn(a)
        aout = sparse(row,col,zeros(nza),m,n)
        for j = 1:nza
            rj = row[j]
            aaj = abs(a[rj,col[j]])
            fw[rj] = max(fw[rj],aaj)
        end
        for k = 1:m
            if fw[k] > 0.0
                fw[k] = 1.0/fw[k]
            else
                fw[k] = 1.0
            end
        end
        for j = 1:nza
            aout[row[j],col[j]] = a[row[j],col[j]]*fw[row[j]]
        end
        a[:,:] = aout
    else
        aout = zeros(a)
        for k = 1:m
            s1 = maximum(abs(a[k,1:n]))
            if s1 > 0.0
                s1 = 1.0/s1
                fw[k] = s1
                aout[k,1:n] = a[k,1:n]*s1
            else
                fw[k] = 1.0
                aout[k,1:n] = a[k,1:n]
            end
        end
        a[:,:] = aout
    end
end

"""
function nScrb(n::Int64, lda::Int64, ml::Int64, mu::Int64, a::Array{Float64,2},
    fw::Vector{Float64})

Row scaling of a (n,n)-matrix in band storage mode

## Input parameters
-------------------
| Variables | Description                              |
|-----------|------------------------------------------|
| n         | Number of rows and columns of the matrix |
| lda       | Leading dimension of the matrix array    |
| ml        | Lower bandwidth of the matrix            |
| mu        | Upper bandwidth of the matrix            |
| a[lda,n]* | Matrix to be scaled                      |

(* marks inout parameters)

## Output parameters
-------------------
| Variables | Description                              |
|-----------|------------------------------------------|
| fw        | Row scaling factors.                     |
"""
function nScrb(n::Int64, lda::Int64, ml::Int64, mu::Int64, a::Array{Float64,2},
    fw::Vector{Float64})
    # Begin
    aout = zeros(a)
    m2 = ml + mu + 1
    for k = 1:n
        s1 = 0.0
        l2 = max(1,k-ml)
        l3 = min(n,k+mu)
        k1 = m2 + k
        for l1 = l2:l3
            s1 = max(s1,abs(a[k1-l1,l1]))
        end
        if s1 > 0.0
            s1 = 1.0/s1
            fw[k] = s1
            for l1 = l2:l3
                aout[k1-l1,l1] = a[k1-l1,l1]*s1
            end
        else
            fw[k] = 1.0
        end
    end
    a[:,:] = aout
end

"""
function nLvls(n::Int64, dxq::Vector{Float64}, dx1::Vector{Float64},
    xw::Vector{Float64}, f::Vector{Float64}, qdscal::Bool)

Provides descaled solution, error norm and level functions
To be used in connection with NLEQ1 and NLEQ2.

## Input parameters
-------------------
| Variables | Description                                   |
|-----------|-----------------------------------------------|
| n         | Number of parameters to be estimated          |
| dx1       | Scaled Newton correction                      |
| xw        | Vector of scaling values                      |
| f         | Residual vector                               |
| qdscal    | true of descaling of dx1 required, else false |

## Output parameters
-------------------
| Variables | Description                              |
|-----------|------------------------------------------|
| dxq       | Leading dimension of the matrix array    |
| conv      | Scaled maximum norm of Newton correction |
| sumX      | Scaled natural level function value      |
| dLevF     | Standard level function value            |

"""
function nLvls(n::Int64, dxq::Vector{Float64}, dx1::Vector{Float64},
    xw::Vector{Float64}, f::Vector{Float64}, qdscal::Bool)
    # Begin
    if qdscal
        # ----------------------------------------------------------------------
        # 1.2 Descaling of solution dx1 (stored to dxq)
        dxq[:] = dx1.*xw
    end
    # --------------------------------------------------------------------------
    # 2 Evaluation of scaled natural level function sumx and scaled maximum
    # error norm conv
    conv = maximum(abs(dx1))
    sumx = sum(dx1.^2)
    # --------------------------------------------------------------------------
    # 3 Evaluation of (scaled) standard level function dlevf
    dlevf = sqrt(sum(f.^2)/n)
    return (conv,sumx,dlevf)
end

# The following function nPrv1 is a multiple dispatch function
# The first one corresponds to nleq1 and the second corresponds to nleq2
"""
function nPrv1(dlevf::Float64, dlevx::Float64, fc::Float64, niter::Int64,
    newt::Int64, mPr::Int64, printIO, qMixIO::Bool)

Printing of intermediate values (Type 1 routine)

## Parameters
-------------
For all the parameters check n1int
"""
function nPrv1(dlevf::Float64, dlevx::Float64, fc::Float64, niter::Int64,
    newt::Int64, mPr::Int64, printIO, qMixIO::Bool)
    # Begin
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
        if mPr >= 3
            write(printIO,"        It       Normf           Normx                     New\n")
        end
        if mPr == 2
            write(printIO,"        It       Normf           Normx         Damp.Fct.   New\n")
        end
    end
    if mPr >= 3 || niter == 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e                 %2i\n",niter,dlevf,dlevx,newt))
    end
    if mPr == 2 && niter != 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e      %7.5f    %2i\n",niter,dlevf,dlevx,fc,newt))
    end
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n")
    end
    flush(printIO)
    return nothing
end

"""
function nPrv1(dlevf::Float64, dlevx::Float64, fc::Float64, niter::Int64, newt::Int64,
    iRank::Int64, mPr::Int64, printIO, qMixIO::Bool, cond1::Float64)

Printing of intermediate values (Type 1 routine)

Parameters
-------------
For all the parameters check n2int
"""
function nPrv1(dlevf::Float64, dlevx::Float64, fc::Float64, niter::Int64, newt::Int64,
    iRank::Int64, mPr::Int64, printIO, qMixIO::Bool, cond1::Float64)
    # Begin
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
        if mPr >= 3
            write(printIO,"        It       Normf           Normx                     New      Rank        Cond\n")
        end
        if mPr == 2
            write(printIO,"        It       Normf           Normx         Damp.Fct.   New      Rank        Cond\n")
        end
    end
    if mPr >= 3 || niter == 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e",niter,dlevf,dlevx),
        @sprintf("                 %2i      %4i          %10.3e\n",newt,iRank,cond1))
    end
    if mPr == 2 && niter != 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e",niter,dlevf,dlevx),
        @sprintf("      %7.5f    %2i      %4i          %10.3e\n",fc,newt,iRank,cond1))
    end
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n")
    end
    flush(printIO)
    return nothing
end

"""
function nPrv2(dlevf::Float64, dlevx::Float64, fc::Float64, niter::Int64,
    printIO, qMixIO::Bool, cmark::AbstractString)

Printing of intermediate values (Type 2 routine)

## Input parameters
-------------------
| Variables | Description                                 |
|-----------|---------------------------------------------|
| dlevf     | Standard level function value               |
| dlevx     | Standard level value                        |
| fc        | Current damping factor                      |
| niter     | Current number of Newton iterations         |
| qMixIO    | Decision parameter for printing             |
| cmark     | Marker character to be printed before dlevx |
"""
function nPrv2(dlevf::Float64, dlevx::Float64, fc::Float64, niter::Int64,
    printIO, qMixIO::Bool, cmark::AbstractString)
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
        write(printIO,"        It       Normf           Normx         Damp.Fct.\n")
    end
    write(printIO,@sprintf("      %4i     %10.3e    %1s %10.3e      %7.5f\n",niter,dlevf,cmark,dlevx,fc))
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n")
    end
    flush(printIO)
    return nothing
end

"""
function nSout(n::Int64, x::Vector{Float64}, mode::Int64, mPr::Int64, printIO,
    nIter::Int64, dLevF::Float64, sumX::Float64)

Printing of iterate (user customizable routine)

## Input parameters
-------------------
| Variables | Description                                                     |
|-----------|-----------------------------------------------------------------|
| n         | Size of the problem                                             |
| x         | Iterate vector                                                  |
| mode = 1  | This routine is called before the first Newton iteration step   |
|      = 2  | This routine is called with an intermediate iterate x           |
|      = 3  | This is the last call with the solution vector x                |
|      = 4  | This is the last call with the final, but not solution vector x |
| mPr       | Decision parameter for printing                                 |
| printIO   | IO handle for printing                                          |
| nIter     | Current number of Newton iterations                             |
| dLevF     | Standard level function value                                   |
| sumX      | Scaled natural level function value                             |
"""
function nSout(n::Int64, x::Vector{Float64}, mode::Int64, mPr::Int64, printIO,
    nIter::Int64, dLevF::Float64, sumX::Float64)
    # Begin
    if mode == 1
        write(printIO,@sprintf("\n%s\n%s%5i\n\n%s\n","  Start data:","  N =",n,
            "  Format: iteration-number, (x(i),i=1,...N) , Normf , Normx "))
        write(printIO,@sprintf("%s\n","  Initial data:"))
    elseif mode == 3
        write(printIO,@sprintf("%s\n","  Solution data:"))
    elseif mode == 4
        write(printIO,@sprintf("%s\n","  Final data:"))
    end
    write(printIO,@sprintf(" %5i\n",nIter))
    l2 = 0
    for l1 = 1:n
        write(printIO,@sprintf("%18.10e ",x[l1]))
        l2 += 1
        if l2 == 3
            write(printIO," \n")
            l2 = 0
        end
    end
    write(printIO,@sprintf("%18.10e %18.10e \n",dLevF,
        sqrt(sumX/n)))
    if mode == 1 && mPr >= 2
        write(printIO,@sprintf("%s\n","  Intermediate data:"))
    elseif mode >= 3
        write(printIO,@sprintf("%s\n","  End data:"))
    end
    flush(printIO)
    return nothing
end

"""
function wnorm(n::Int64, z::Vector{Float64}, xw::Vector{Float64})

Return the norm to be used in exit (termination) criteria

## Input parameters
-------------------
| Variables | Description                                                     |
|-----------|-----------------------------------------------------------------|
| n         | Size of the problem                                             |
| z         | The vector of which norm is to be computed                      |
| xw        | Scaling value of z                                              |

## Output
---------
The mean square root norm of z subject to the scaling values in xw.
"""
function wnorm(n::Int64, z::Vector{Float64}, xw::Vector{Float64})
    # Begin
    return sqrt(sum((z./xw).^2)/n)
end
