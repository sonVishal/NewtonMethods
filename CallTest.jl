include("NLEQ1.jl")

function f(x)
    return x;
end

function Jac(x)
    return x;
end

opt = OptionsNLEQ(OPT_PRINTWARNING => 1, OPT_RTOL => 1);

x = ones(4);
xScal = ones(4);

(stats,retCode) = nleq1(f,x,xScal,opt);

println("retCode = $retCode","\n");
println("Options = $opt","\n");
