"""
# Summary:
checkOptions : Checking of input parameters and options for NLEQ1.
"""
function checkOptions(n::Int64, x::Vector{Float64}, xScal::Vector{Float64},
    opt::OptionsNLEQ)
    # Begin

    # Machine related constants
    epMach  = 1e-17
    small   = 1e-150
    great   = 1.0/small

    # Check the print stream variables
    printWarn   = getOption!(opt, OPT_PRINTWARNING, 0)
    printIOwarn = getOption!(opt, OPT_PRINTIOWARN, 0)
    if printIOwarn == 0
        printIOwarn = STDOUT
    elseif typeof(printIOwarn) != IOStream
        write(STDOUT, "ERROR: Please provide a file handle for writing warning messages",
            " \n or use the default option OPT_PRINTIOWARN = 0 for output to STDOUT.\n")
        return 99
    end
    printIOmon  = getOption!(opt, OPT_PRINTIOMON, 0)
    if printIOmon == 0
        printIOmon = STDOUT
    elseif typeof(printIOmon) != IOStream
        write(STDOUT, "ERROR: Please provide a file handle for writing iteration monitor",
            " \n or use the default option OPT_PRINTIOMON = 0 for output to STDOUT.\n")
        return 99
    end
    printIOsol  = getOption!(opt, OPT_PRINTIOSOL, 0)
    if printIOsol == 0
        printIOsol = STDOUT
    elseif typeof(printIOsol) != IOStream
        write(STDOUT, "ERROR: Please provide a file handle for writing the solution",
            " \n or use the default option OPT_PRINTIOSOL = 0 for output to STDOUT.\n")
        return 99
    end

    # Check dimensional parameter n
    if n <= 0
        write(printIOwarn, "ERROR: Bad input to dimensional parameter n supplied","\n",
            "Choose n positive, your input is: n = $n\n")
        return 20
    end

    # Check if x and xScal are of same length
    if length(x) != length(xScal)
        write(printIOwarn, "ERROR: Length of x and xScal do not match","\n")
        return 20
    end

    # Check for successive call
    qSucc = getOption!(opt, OPT_QSUCC, 0)
    if qSucc > 1 || qSucc < 0
        write(printIOwarn, "ERROR: Invalid option for OPT_QSUCC provided.\n")
        return 99
    end

    # Check the mode
    mode = getOption!(opt, OPT_MODE, 0)
    if mode > 1 || mode < 0
        write(printIOwarn, "ERROR: Invalid option for OPT_MODE provided.\n")
        return 99
    end

    # Check the Jacobian
    jacGen = getOption!(opt, OPT_JACGEN, 2)
    if jacGen == 1
        jacFcn = getOption!(opt, OPT_JACFCN, 0)
        if jacFcn == 0
            write(printIOwarn, "ERROR: The Jacobian function OPT_JACFCN is not supplied. ",
            "Please supply a Jacobian function or use OPT_JACGEN = 2 or 3 for numerical differentiation based jacobian evaluation.\n")
            return 99
        end
    elseif jacGen < 1 || jacGen > 3
        write(printIOwarn, "ERROR: Invalid option for OPT_JACGEN provided.\n")
        return 99
    end

    # Check the storage
    mStor = getOption!(opt, OPT_MSTOR, 0)
    if mStor < 0 || mStor > 1
        write(printIOwarn, "ERROR: Invalid option for OPT_MSTOR provided.\n")
        return 99
    end

    # Check the bandwidth if in banded mode
    if mStor == 1
        ml = getOption!(opt, OPT_ML, 0)
        mu = getOption!(opt, OPT_MU, 0)
        if ml < 0 || ml > 9999999
            write(printIOwarn, "ERROR: Invalid option for OPT_ML provided.\n")
            return 99
        end
        if mu < 0 || mu > 9999999
            write(printIOwarn, "ERROR: Invalid option for OPT_MU provided.\n")
            return 99
        end
    end

    # Check scaling strategy
    iScal = getOption!(opt, OPT_ISCAL, 0)
    if iScal < 0 || iScal > 1
        write(printIOwarn, "ERROR: Invalid option for OPT_ISCAL provided.\n")
        return 99
    end

    # Problem type specification by user
    nonLin = getOption!(opt, OPT_NONLIN, 3)
    if nonLin < 1 || nonLin > 4
        write(printIOwarn, "ERROR: Invalid option for OPT_NONLIN provided.\n")
        return 99
    end

    # Check rank-1 update strategy
    qRank1 = getOption!(opt, OPT_QRANK1, 0)
    if qRank1 < 0 || qRank1 > 1
        write(printIOwarn, "ERROR: Invalid option for OPT_QRANK1 provided.\n")
        return 99
    end

    # Check if ordinary Newton iteration is to be performed
    qOrdi = getOption!(opt, OPT_QORDI, 0)
    if qOrdi < 0 || qOrdi > 1
        write(printIOwarn, "ERROR: Invalid option for OPT_ORDI provided.\n")
        return 99
    end

    # Check if Simplified Newton iteration is to be performed
    qSimpl = getOption!(opt, OPT_QSIMPL, 0)
    if qSimpl < 0 || qSimpl > 1
        write(printIOwarn, "ERROR: Invalid option for OPT_SIMPL provided.\n")
        return 99
    end

    # Check if automatic row scaling is allowed ot inhibited
    qNScal = getOption!(opt, OPT_NOROWSCAL, 0)
    if qNScal < 0 || qNScal > 1
        write(printIOwarn, "ERROR: Invalid option for OPT_NOROWSCAL provided.\n")
        return 99
    end

    # Check bounded damping strategy
    iBDamp = getOption!(opt, OPT_BOUNDEDDAMP, 0)
    if iBDamp < 0 || iBDamp > 2
        write(printIOwarn, "ERROR: Invalid option for OPT_BOUNDEDDAMP provided.\n")
        return 99
    end

    # Check convergence order monitor
    iOrMon = getOption!(opt, OPT_IORMON, 2)
    if iOrMon < 1 || iOrMon > 3
        write(printIOwarn, "ERROR: Invalid option for OPT_IORMON provided.\n")
        return 99
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

    # If this point is reached means everything went on smoothly
    return 0
end

function initializeOptions(opt, wk, n, m1, qRank1)
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
    initOption!(wk, WK_QU, zeros(n))
    initOption!(wk, WK_T1, zeros(n))
    initOption!(wk, WK_T2, zeros(n))

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

function printInitialization(n, printIOmon, rTol, jacGen, mStor, ml, mu,
    qNoRowScal, qRank1, nonLin, qBDamp, fcBand, qOrdi, qSimpl, nItmax)
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

function printStats(stats, printIOmon)
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

# TODO: Get printIO from call rather than inside function
function nScal(n,x,xa,xScal,iScal,qIniSc,opt)
    small = getMachineConstants(6)
    # Begin
    xw = zeros(n)
    if iScal == 1
        xw[:] = xScal
    else
        for l1 = 1:n
            xw[l1] = max(xScal[l1],max((abs(x[l1])+abs(xa[l1]))*0.5,small))
        end
    end

    mPr = opt.options[OPT_PRINTITERATION]
    if mPr >= 6
        printIO = opt.options[OPT_PRINTIO]
        write(printIO,"\n\n",
        "+++++++++++++++++++++++++++++++++++++++++++++++++\n",
        "      x-components         scaling-components\n")
        for l1 = 1:n
            write(printIO,"  %18.10e   %18.10e\n",x[l1],xw[l1])
        end
        write(printIO,"+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
    end
    return xw
end

function nScrf(m,n,a)
    # Begin
    fw = zeros(n)
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
    return (aout,fw)
end

function nScrb(n,lda,ml,mu,a)
    # Begin
    fw = zeros(n)
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
    return (aout,fw)
end

function nLvls(n,dxq,dx1,xw,f,mPr,qdscal)
    # Begin
    if qdscal
        # ----------------------------------------------------------------------
        # 1.2 Descaling of solution dx1 (stored to dxq)
        dxq = dx1.*xw
    end
    # --------------------------------------------------------------------------
    # 2 Evaluation of scaled natural level function sumx and scaled maximum
    # error norm conv
    conv = maximum(abs(dx1))
    sumx = sum(dx1.^2)
    # --------------------------------------------------------------------------
    # 3 Evaluation of (scaled) standard level function dlevf
    dlevf = sqrt(sum(f.^2)/n)
    return (dxq,conv,sumx,dlevf)
end


function nPrv2(dlevf,dlevx,fc,niter,mPr,printIO,qMixIO,cmark)
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

function nSout(n,x,mode,opt,wkNLEQ1,mPr,printIO)
    # Begin
    qNorm = true
    if qNorm
        if mode == 1
            write(printIO,@sprintf("%s\n%s%5i\n\n%s\n","  Start data:","  N =",n, "  Format: iteration-number, (x(i),i=1,...N) , Normf , Normx "))
            write(printIO,@sprintf("%s\n","  Initial data:"))
        elseif mode == 3
            write(printIO,@sprintf("%s\n","  Solution data:"))
        elseif mode == 4
            write(printIO,@sprintf("%s\n","  Final data:"))
        end
        write(printIO,@sprintf(" %5i\n",wkNLEQ1.options[STATS_NITER]))
        l2 = 0
        for l1 = 1:n
            write(printIO,@sprintf("%18.10e ",x[l1]))
            l2 += 1
            if l2 == 3
                write(printIO,"%1s\n"," ")
                l2 = 0
            end
        end
        write(printIO,@sprintf("%18.10e %18.10e \n",wkNLEQ1.options[STATS_DLEVF],sqrt(wkNLEQ1.options[STATS_SUMX]/n)))
        if mode == 1 && mPr >= 2
            write(printIO,@sprintf("%s\n","  Intermediate data:"))
        elseif mode >= 3
            write(printIO,@sprintf("%s\n","  End data:"))
        end
    end
    return nothing
end

function wnorm(n,z,xw)
    # Begin
    return sqrt(sum((z./xw).^2)/n)
end
