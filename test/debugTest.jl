using NewtonMethods
include("ChebyQuad.jl")
clearWorkspace("nleq2")
dim = 2
n1 = dim + 1

fSol = open("nleq2.dat","w")
fRest = open("nleq2.out","w")

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
    println("returned from call $i")
    i += 1
end

flush(fSol)
flush(fRest)

close(fSol)
close(fRest)
