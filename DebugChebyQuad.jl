include("ChebyQuad.jl");

# Dimension of the problem 2 <= dim <= 9
dim = 2
n1 = dim + 1;

# Get the Cheby Quad function
# (name,fcn,jac,x) = chebyquad(dim)



# Initialize the options
opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_JACFCN            => chebyQuadJac,
                  OPT_MSTOR             => 0,
                  OPT_NOROWSCAL         => 0,
                  OPT_PRINTIO           => STDOUT,
                  OPT_PRINTWARNING      => 1,
                  OPT_PRINTITERATION    => 3,
                  OPT_PRINTSOLUTION     => 2,
                  OPT_NITMAX            => 10,
                  OPT_RTOL              => 1e-5)

x0    = collect(1:dim)./n1
xScal = zeros(x0)

# f = zeros(x0)
# J = zeros(length(x0),length(x0))
#
# chebyQuad(x0,f)
# chebyQuadJac(x0,J)
#
# println(f)
# println(J)

(sol, stats, retCode) = nleq1(chebyQuad,x0,xScal,opt);

println("DONE")
