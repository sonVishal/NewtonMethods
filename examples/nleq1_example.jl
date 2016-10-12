using NewtonMethods
include("ChebyPoly.jl")

fSol = open("nleq1.dat","w")
fRest = open("nleq1.out","w")
dimMax = 9
for dim = 2:dimMax
    n1 = dim + 1

    # Initialize the options
    opt = OptionsNLEQ(OPT_MODE              => 1,
                      OPT_JACGEN            => 1,
                      OPT_JACFCN            => chebyQuadJac,
                      OPT_MSTOR             => 0,
                      OPT_NOROWSCAL         => 0,
                      OPT_PRINTWARNING      => 1,
                      OPT_PRINTITERATION    => 3,
                      OPT_PRINTSOLUTION     => 2,
                      OPT_PRINTIOWARN       => fRest,
                      OPT_PRINTIOMON        => fRest,
                      OPT_PRINTIOSOL        => fSol,
                      OPT_NITMAX            => 10,
                      OPT_RTOL              => 1e-5)

    x0    = collect(1:dim)./n1
    xScal = zeros(x0)

    retCode = -1

    i = 1

    println("Calling the while loop for solving the Cheby Quad equation of dimension $dim")
    while retCode == -1
        (x0, _, retCode) = nleq1(chebyQuad, x0, xScal, opt)
        write(fRest, @sprintf("Returned from call %4i of NLEQ1\n",i))
		flush(fSol)
	    flush(fRest)
        i += 1
    end
    println("Solution = $x0")

    println("DONE")
end
close(fSol)
close(fRest)
