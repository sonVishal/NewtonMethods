using NewtonMethods
include("ChebyQuad.jl")
dim = 8
n1 = dim + 1

fSol = open("nleq2_$dim.dat","w")
fRest = open("nleq2_$dim.out","w")

# Initialize the options
opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_PRINTWARNING      => 1,
                  OPT_PRINTITERATION    => 3,
                  OPT_PRINTSOLUTION     => 2,
                  OPT_PRINTIOWARN       => fRest,
                  OPT_PRINTIOMON        => fRest,
                  OPT_PRINTIOSOL        => fSol,
                  OPT_JACFCN            => chebyQuadJac,
                  OPT_MSTOR             => 0,
                  OPT_NOROWSCAL         => 0,
                  OPT_NITMAX            => 200,
                  OPT_RTOL              => 1e-5)

x0    = collect(1:dim)./n1
xScal = zeros(x0)

retCode = -1
stats   = []

i = 1

while retCode == -1
    (x0, stats, retCode) = nleq2(chebyQuad, x0, xScal, opt)
    write(fRest, @sprintf("Returned from call %4i of NLEQ2\n",i))
    flush(fSol)
    flush(fRest)
    i += 1
end

close(fSol)
close(fRest)
