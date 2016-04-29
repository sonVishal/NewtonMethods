include("CheckOptionsNLEQ1.jl");

# Bookkeeping
# TODO: Make everything type independent
# Currently everything is assumed to be Float64
function nleq1(f,x,xscal,opt::OptionsNLEQ)
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
        return retCode
    end
    return retCode
end
