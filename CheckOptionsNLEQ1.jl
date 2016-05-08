include("Options.jl");
include("MachineConsts.jl");
# Check options if they are correct
function checkOptions(n,x,xScal,opt::OptionsNLEQ)

    # TODO: Get the type of elements in x

    # Check for printing options
    # TODO: Take input (using options) on where to print
    # Currently printing default to stdout
    printFlag = getOption(opt,"PR_ERR",0);

    # Get the machine related constants
    epMach = getMachineConstants(3);
    great = sqrt(getMachineConstants(2)/10.0);
    small = getMachineConstants(6);

    # Initialize return code to 0
    retCode = 0;

    # Check dimensional parameter n
    if n <= 0
        if printFlag >= 1
            println("Error: bad input to dimensional parameter n supplied","\n",
                "Choose n positive, your input is: n = $n");
        end
        retCode = 20;
    end

    # Problem type specification by user
    nonLin = getOption(opt,"NONLIN",0);
    if nonLin == 0
        nonLin = 3;
    end
    setOption!(opt,"NONLIN",nonLin);

    # Checking and conditional adaptation of user given RTOL
    rTol = getOption(opt,"RTOL",0.0);
    # if RTOL is not set, set it to 1e-6
    if rTol == 0.0
        rTol = 1e-6;
    elseif rTol < 0
        if printFlag >= 1
            println("Error: Nonpositive RTOL supplied");
        end
        retCode = 21;
    else
        tolMin = epMach*10.0*n;
        if rTol < tolMin
            rTol = tolMin;
            if printFlag >= 2
                println("Warning: user prescribed RTOL increased to
                    a reasonable smallest value RTOL = $rtol");
            end
        end

        tolMax = 1.0e-1;
        if rTol > tolMax
            rTol = tolMax;
            if printFlag >= 2
                println("Warning: user prescribed RTOL decreased to
                    a reasonable largest value RTOL = $rtol");
            end
        end
    end

    # Test user given accuracy and scaling on proper values
    if nonLin >= 3
        defScal = rTol;
    else
        defScal = 1.0;
    end

    for i = 1:n
        # Scaling Values cannot be negative
        # Positive scaling values give scale invariance
        if xScal[i] < 0.0
            if printFlag >= 1
                println("Error: negative value in xScal[$i] supplied");
            end
            retCode = 22;
        end

        if xScal[i] == 0.0
            xScal[i] = defScal;
        end
        # Avoid overflow due to division by xScal[i]
        if xScal[i] > 0.0 && xScal[i] < small
            if printFlag >= 2
                println("Warning: xScal[$i] = $xScal[i] too small, increased to $small");
            end
            xScal[i] = small;
        end
        # Avoid underflow due to division by xScal[i]
        if xScal[i] > great
            if printFlag >= 2
                println("Warning: xScal[$i] = $xScal[i] too big, increased to $great");
            end
            xScal[i] = great;
        end
    end

    return retCode;
end
