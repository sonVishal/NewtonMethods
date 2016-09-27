dim = 2
n1 = dim + 1

# Initialize the options
opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_JACFCN            => chebyQuadJac,
                  OPT_MSTOR             => 0,
                  OPT_NOROWSCAL         => 0,
                  OPT_NITMAX            => 10,
                  OPT_RTOL              => 1e-5)

x0    = collect(1:dim)./n1
xScal = zeros(x0)

retCode = -1
stats   = []

while retCode == -1
    (x0, stats, retCode) = nleq1(chebyQuad,x0,xScal,opt)
end

refSol = [0.211324865405173;
            0.788675134594827]

relNormDiff = norm(x0-refSol)/norm(refSol)

println(relNormDiff)
