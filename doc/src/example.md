# Example for NLEQ1 solver

Please check the example folder for the *.jl file.

```julia
# Import the NewtonMethods package
using NewtonMethods

# Include the Chebyshev polynomial function and Jacobian generating functions
include("ChebyPoly.jl")

# Open files for writing
fSol  = open("nleq1.dat","w")
fRest = open("nleq1.out","w")

# We will solve Chebyshev polynomials of dimensions 2 to 9
dimMax = 9
for dim = 2:dimMax
    n1 = dim + 1

    # Initialize the options
    opt = OptionsNLEQ(OPT_MODE              => 1, # Use stepwise mode
                      OPT_JACGEN            => 1, # Use Jacobian function provided
                      OPT_JACFCN            => chebyQuadJac, # User defined Jacobian function
                      OPT_MSTOR             => 0, # Full storage mode
                      OPT_NOROWSCAL         => 0, # Inhibit automatic row scaling
                      OPT_PRINTWARNING      => 1, # Print all warnings
                      OPT_PRINTITERATION    => 3, # Print intermediate iterations
                      OPT_PRINTSOLUTION     => 2, # Print solution
                      OPT_PRINTIOWARN       => fRest, # Where to print warnings?
                      OPT_PRINTIOMON        => fRest, # Where to print iterations?
                      OPT_PRINTIOSOL        => fSol,  # Where to print solution?
                      OPT_NITMAX            => 10, # Maximum number of iterations
                      OPT_RTOL              => 1e-5) # Required tolerance

    x0    = collect(1:dim)./n1 # Initial guess
    xScal = zeros(x0) # Initial scaling vector

    retCode = -1 # Initialize return code to -1.

    i = 1 # Counter for iterations (optional)

    println("Calling the while loop for solving the Cheby Quad equation of dimension $dim")
    # If return code is again -1 then repeat.
    while retCode == -1
        # Perform one Newton iteration
        (x0, stats, retCode) = nleq1(chebyQuad, x0, xScal, opt)
        write(fRest, @sprintf("Returned from call %4i of NLEQ1\n",i))
		flush(fSol)
	    flush(fRest)
        i += 1
    end
    # Print solution or final iterate
    println("Solution = $x0")

    println("DONE")
end

# Close the files that were opened for writing
close(fSol)
close(fRest)
```

Example for NLEQ2 is also included in the `examples` folder.
