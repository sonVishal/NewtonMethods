include("ChebyQuad.jl");

# Dimension of the problem 2 <= dim <= 9
dim = 2;
n1 = dim+1;

# Get the Cheby Quad function
(name,fcn,jac,x) = chebyquad(dim);

# Initialize the options
opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_JACFCN            => jac,
                  OPT_MSTOR             => 0,
                  OPT_PRINTIO           => STDOUT,
                  OPT_PRINTWARNING      => 1,
                  OPT_PRINTITERATION    => 3,
                  OPT_PRINTSOLUTION     => 2,
                  OPT_NITMAX            => 10,
                  OPT_RTOL              => 1e-5);

x0    = collect(1:dim)./n1;
xScal = zeros(x0);

(sol, stats, retCode) = nleq1(fcn,x0,xScal,opt);

println("DONE")
