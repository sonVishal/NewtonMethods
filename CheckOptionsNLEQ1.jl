include("Options.jl")
include("MachineConsts.jl")

# Summary :
# checkOptions : Checking of input parameters and options for NLEQ1.

function checkOptions(n,x,xScal,opt::OptionsNLEQ)

    # TODO: Get the type of elements in x

    # TODO: Logging and printing of Logs!!

    # Check whether warnings need to be printed
    warnFlag = getOption(opt,"OPT_PRINT_WARNINGS",0)

    # Get the machine related constants
    epMach  = getMachineConstants(3)
    great   = sqrt(getMachineConstants(2)/10.0)
    small   = getMachineConstants(6)

    # Initialize return code to 0
    retCode = 0

    # Check dimensional parameter n
    if n <= 0
        retCode = 20
        println("ERROR: Bad input to dimensional parameter n supplied","\n",
            "Choose n positive, your input is: n = $n")
        return retCode
    end

    # Problem type specification by user
    nonLin = getOption(opt,"OPT_NONLIN",0)
    if nonLin == 0
        nonLin = 3
        setOption!(opt,"OPT_NONLIN",nonLin)
    end

    # Checking and conditional adaptation of user given RTOL
    # if RTOL is not set, set it to 1e-6
    rTol = getOption!(opt,"OPT_RTOL",1e-6)
    if rTol <= 0.0
        retCode = 21
        println("ERROR: Nonpositive OPT_RTOL supplied")
        return retCode
    else
        tolMin = epMach*10.0*n
        if rTol < tolMin
            rTol = tolMin
            if warnFlag == 1
                warn("User prescribed OPT_RTOL increased to a reasonable smallest value RTOL = $rTol")
            end
        end

        tolMax = 1.0e-1
        if rTol > tolMax
            rTol = tolMax
            if warnFlag == 1
                warn("User prescribed OPT_RTOL decreased to a reasonable largest value RTOL = $rTol")
            end
        end
        setOption!(opt,"OPT_RTOL",rTol)
    end

    # Test user given accuracy and scaling on proper values
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
            println("ERROR: Negative value in xScal[$i] supplied")
            return retCode
        end

        if xScal[i] == 0.0
            xScal[i] = defScal
        end
        # Avoid overflow due to division by xScal[i]
        if xScal[i] > 0.0 && xScal[i] < small
            if warnFlag == 1
                warn("xScal[$i] = $xScal[i] too small, increased to $small")
            end
            xScal[i] = small
        end
        # Avoid underflow due to division by xScal[i]
        if xScal[i] > great
            if warnFlag == 1
                warn("xScal[$i] = $xScal[i] too big, increased to $great")
            end
            xScal[i] = great
        end
    end

    # Define lower and upper bounds for options
    L = [0; 0; 0; 0];

    # Check if the Jacobian is Dense/Sparse or Banded matrix
    mstor = getOption!(opt,"OPT_MSTOR",0)
    if mstor == 0   # Dense or Sparse
        m1 = 0
        m2 = 0
    elseif mstor == 1   # Banded
        m1 = n - 1
        m2 = n - 1
    else
        m1 = 9999999
        m2 = 9999999
    end

    # TODO: Consider having Jacobian option similar to Dr. Ludwig

    qSucc   = getOption!(opt,"OPT_QSUCC",0)
    mode    = getOption!(opt,"OPT_MODE",0)
    jacGen  = getOption!(opt,"OPT_JACGEN",0)

    # Check bounds of options
    if getOption!(opt,"OPT_QSUCC",0) < 0 || getOption!(opt,"OPT_QSUCC",0) > 1


    return retCode
end
