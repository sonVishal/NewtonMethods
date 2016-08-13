include("ChebyQuad.jl");

# Dimension of the problem 2 <= dim <= 9
dim = 2
n1 = dim + 1;

# Get the Cheby Quad function
# (name,fcn,jac,x) = chebyquad(dim)

fSol = open("solution.dat","w")
fRest = open("solution.out","w")

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
                  OPT_RTOL              => 1e-8)

x0    = collect(1:dim)./n1
xScal = zeros(x0)

retCode = -1
stats   = []

println("Calling the while loop for solving the Cheby Quad equation of dimension $dim");
while retCode == -1
    (x0, stats, retCode) = nleq1(chebyQuad,x0,xScal,opt);
end
println("Solution = $x0")
flush(fSol)
flush(fRest)
close(fSol)
close(fRest)

println("DONE")
