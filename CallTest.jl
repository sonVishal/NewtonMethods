# Initialize the options
function f(x,y)
    y[1] = sum(x.^4) - 32;
    y[2] = sum(x.^5) - 186;
    y[3] = sum(x.^6) - 803;
    return nothing;
end

function Df(x,J)
    J[1,1:3] = 4*x.^3;
    J[2,1:3] = 5*x.^4;
    J[3,1:3] = 6*x.^5;
    return nothing
end

opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_JACFCN            => Df,
                  OPT_MSTOR             => 0,
                  OPT_NOROWSCAL         => 0,
                  OPT_NITMAX            => 10,
                  OPT_RTOL              => 1e-3);

x0    = ones(3);
xScal = zeros(x0)

retCode = -1
stats   = []

println("Calling the while loop for solving the Test equation");
while retCode == -1
    (x0, stats, retCode) = nleq1(f,x0,xScal,opt);
end
println("Solution = $x0")
println("DONE")
