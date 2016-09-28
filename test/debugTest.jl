using NewtonMethods
include("ChebyQuad.jl")
dim = 2
n1 = dim + 1

# Initialize the options
opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_JACFCN            => chebyQuadJac,
                  OPT_MSTOR             => 0,
                  OPT_NOROWSCAL         => 0,
                  OPT_NITMAX            => 200,
                  OPT_RTOL              => 1e-5)

x0    = collect(1:dim)./n1
xScal = zeros(x0)

retCode = -1
stats   = []

count = 1

while retCode == -1
    (x0, stats, retCode) = nleq2(chebyQuad, x0, xScal, opt)
    println("returned from call $count")
    count += 1
end

f = zero(x0)

chebyQuad(f,x0)

println(f)
