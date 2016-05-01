include("NLEQ1.jl")

opt = OptionsNLEQ();
wk  = OptionsNLEQ();

function f(x)
    return x;
end

function Jac(x)
    return x;
end

x = ones(4);
xScal = ones(4);

(stats,retCode) = nleq1(f,x,xScal,opt,wk);

println("retCode = $retCode","\n");
println("Options = $opt","\n");
println("WK = $wk");
