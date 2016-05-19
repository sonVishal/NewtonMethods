include("Options.jl")
include("Constants.jl")

# Summary :
# checkOptions : Checking of input parameters and options for NLEQ1.

function checkOptions(n::Number,x::Vector,xScal::Vector,opt::OptionsNLEQ)

    # TODO: Get the type of elements in x

    # Check whether warnings need to be printed
    warnFlag = getOption(opt,OPT_PRINTWARNING,0)

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
        println("ERROR: Nonpositive $OPT_RTOL supplied")
        return retCode
    else
        tolMin = epMach*10.0*n
        if rTol < tolMin
            rTol = tolMin
            setOption!(opt,OPT_RTOL,rTol)
            if warnFlag == 1
                warn("User prescribed $OPT_RTOL increased to a reasonable smallest value RTOL = $rTol")
            end
        end

        tolMax = 1.0e-1
        if rTol > tolMax
            rTol = tolMax
            setOption!(opt,OPT_RTOL,rTol)
            if warnFlag == 1
                warn("User prescribed $OPT_RTOL decreased to a reasonable largest value RTOL = $rTol")
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

    # TODO: Not in MATLAB file but not sure if the other options should be tested
    # # Define lower and upper bounds for options
    # optL = zeros(8);
    # optU = [1; 1; 3; 1; 9999999; 9999999; 1; 1];
    #
    # # Check if the Jacobian is Dense/Sparse or Banded matrix
    # mStor = getOption!(opt,OPT_MSTOR,0)
    # if mStor == 0   # Dense or Sparse
    #     optU[6] = 0
    #     optU[7] = 0
    # elseif mStor == 1   # Banded
    #     optU[6] = n - 1
    #     optU[7] = n - 1
    # end
    #
    # qSucc   = getOption!(opt,OPT_QSUCC,0)
    # mode    = getOption!(opt,OPT_MODE,0)
    # jacGen  = getOption!(opt,OPT_JACOBIFCN,0)
    # ml      = getOption!(opt,OPT_ML,0)
    # mu      = getOption!(opt,OPT_MU,0)
    # iScal   = getOption!(opt,OPT_ISCAL,0)
    #
    # # Check bounds of options
    # # TODO: find a nicer way to do this
    # if qSucc < optL[1] || qSucc > optU[1]
    #     retCode = 30
    #     error("Invalid option specified: OPT_QSUCC = $qSucc
    #     range of permitted value is $(optL[1]) to $(optU[1])")
    # end
    # if mode < optL[2] || mode > optU[2]
    #     retCode = 30
    #     error("Invalid option specified: OPT_MODE = $mode
    #     range of permitted value is $(optL[2]) to $(optU[2])")
    # end
    # if typeof(jacGen) == Number
    #     if jacGen < optL[3] || jacGen > optU[3]
    #         retCode = 30
    #         error("Invalid option specified: OPT_JACOBIFCN = $jacGen
    #         range of permitted value is $(optL[3]) to $(optU[3])")
    #     end
    # end
    # if mStor < optL[4] || qSucc > optU[4]
    #     retCode = 30
    #     error("Invalid option specified: OPT_MSTOR = $mStor
    #     range of permitted value is $(optL[4]) to $(optU[4])")
    # end
    # if ml < optL[5] || ml > optU[5]
    #     retCode = 30
    #     error("Invalid option specified: OPT_ML = $ml
    #     range of permitted value is $(optL[5]) to $(optU[5])")
    # end
    # if mu < optL[6] || mu > optU[6]
    #     retCode = 30
    #     error("Invalid option specified: OPT_MU = $mu
    #     range of permitted value is $(optL[6]) to $(optU[6])")
    # end
    # if iScal < optL[7] || iScal > optU[7]
    #     retCode = 30
    #     error("Invalid option specified: OPT_ISCAL = $iScal
    #     range of permitted value is $(optL[7]) to $(optU[7])")
    # end
    # if warnFlag < optL[8] || qSucc > optU[8]
    #     retCode = 30
    #     error("Invalid option specified: OPT_PRINTWARNING = $warnFlag
    #     range of permitted value is $(optL[8]) to $(optU[8])")
    # end
    return retCode
end
