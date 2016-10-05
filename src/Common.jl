"""
# Summary:
checkOptions : Checking of common input parameters and options.

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
    # TODO: Get the type of elements in x

    # Check whether warnings need to be printed
    printWarn = getOption!(opt,OPT_PRINTWARNING,0)
    printIOwarn = getOption!(opt,OPT_PRINTIOWARN,STDOUT)

    # Get the machine related constants
    great   = 1.0/small

    # Initialize return code to 0
    retCode = 0

    # Check dimensional parameter n
    if n <= 0
        retCode = 20
        write(printIOwarn,"ERROR: Bad input to dimensional parameter n supplied","\n",
            "Choose n positive, your input is: n = $n")
        return retCode
    end

    # Problem type specification by user
    nonLin = getOption(opt,OPT_NONLIN,0)
    if nonLin == 0
        nonLin = 3
        setOption!(opt,OPT_NONLIN,nonLin)
    end

    # Checking and conditional adaptation of user given RTOL
    # if RTOL is not set, set it to 1e-6
    rTol = getOption!(opt,OPT_RTOL,1e-6)
    if rTol <= 0.0
        retCode = 21
        write(printIOwarn,"ERROR: Nonpositive $OPT_RTOL supplied")
        return retCode
    else
        tolMin = epMach*10.0*n
        if rTol < tolMin
            rTol = tolMin
            setOption!(opt,OPT_RTOL,rTol)
            if printWarn == 1
                write(printIOwarn,"WARNING: User prescribed $OPT_RTOL increased to a reasonable smallest value RTOL = $rTol")
            end
        end

        tolMax = 1.0e-1
        if rTol > tolMax
            rTol = tolMax
            setOption!(opt,OPT_RTOL,rTol)
            if printWarn == 1
                write(printIOwarn,"WARNING: User prescribed $OPT_RTOL decreased to a reasonable largest value RTOL = $rTol")
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
            write(printIOwarn,"ERROR: Negative value in xScal[$i] supplied")
            return retCode
        end

        if xScal[i] == 0.0
            xScal[i] = defScal
        end
        # Avoid overflow due to division by xScal[i]
        if xScal[i] > 0.0 && xScal[i] < small
            if printWarn == 1
                write(printIOwarn,"WARNING: xScal[$i] = $xScal[i] too small, increased to $small")
            end
            xScal[i] = small
        end
        # Avoid underflow due to division by xScal[i]
        if xScal[i] > great
            if printWarn == 1
                write(printIOwarn,"WARNING: xScal[$i] = $xScal[i] too big, increased to $great")
            end
            xScal[i] = great
        end
    end

    # Assign the Jacobian depending on user input
    # Multiple dispatch calls the required function based on
    # the storage requirement of the user
    jacGen = getOption!(opt,OPT_JACGEN,0)
    if jacGen == 1
        jacFcn = getOption!(opt,OPT_JACFCN,0)
        if jacFcn == 0
            retCode = 99
            write(printIOwarn,"ERROR: The Jacobian function OPT_JACFCN is not supplied. ",
            "Please supply a Jacobian function or use OPT_JACGEN = 2 or 3 for numerical differentiation based jacobian evaluation.")
            return retCode
        end
    elseif jacGen == 0
        jacGen = 2
        opt.options[OPT_JACGEN] = jacGen
    end

    # TODO: Add check for IO

    return retCode
end

"""
# Summary:
initializeOptions : Initialization of options based on the solver input argument.

## Input parameters
-------------------
| Variables | Description                                                    |
|-----------|----------------------------------------------------------------|
| opt*      | Options set by the user                                        |
| wk*       | Internal workspace specific to the solver                      |
| n         | Size of the problem                                            |
| m1        | In full mode = n and in band mode = 2*ml+mu+1                  |
| nBroy     | Maximum number of possible consecutive iterative Broyden steps |
| qRank1    | Decision parameter for Rank-1 updates                          |
| solver = 1| Specifies that the solver is NLEQ1                             |
|        = 2| Specifies that the solver is NLEQ2                             |

(* marks inout parameters)
"""
function initializeOptions(opt::OptionsNLEQ, wk::OptionsNLEQ,
    n::Int64, m1::Int64, nBroy::Int64, qRank1::Bool, solver::Int64)
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

    return nothing
end

"""
# Summary:
printInitialization : Print a summary of the initialization.

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
    mStor::Int64, ml::Int64, mu::Int64, qNoRowScal::Bool, qRank1::Bool, nonLin::Int64,
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
# Summary:
printStats : Print a summary of the statistics.

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
# Summary :
nScal : To be used in connection with NLEQ1 and NLEQ2.
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
# Summary :
nScrf : Row scaling of a (m,n)-matrix in full storage mode

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
    end
    a = aout[:,:]
end

"""
# Summary :
nScrb : Row scaling of a (n,n)-matrix in band storage mode

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
    a[:] = aout
end

"""
# Summary :
nLvls : Provides descaled solution, error norm and level functions
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

"""
# Summary :
nPrv2 : Printing of intermediate values (Type 2 routine)

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
        "\n");
    end
    return nothing
end

"""
# Summary :
nSout : Printing of iterate (user customizable routine)
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
    return nothing
end

"""
# Summary :
wnorm : Return the norm to be used in exit (termination) criteria

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
