function f(x,y)
    y = x^2-2*x+1;
end

function Df(x,J)
    J = 2*x-2;
end

opt = OptionsNLEQ(OPT_MODE              => 1,
                  OPT_JACGEN            => 1,
                  OPT_JACFCN            => Df,
                  OPT_MSTOR             => 0,
                  OPT_NOROWSCAL         => 0,
                  OPT_NITMAX            => 10,
                  OPT_RTOL              => 1e-3);

x0    = 0.0
xScal = 0.0

retCode = -1
stats   = []

println("Calling the while loop for solving the Test equation");
while retCode == -1
    (x0, stats, retCode) = nleq1(f,x0,xScal,opt);
end
println("Solution = $x0")
println("DONE")
