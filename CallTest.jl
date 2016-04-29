include("NLEQ1.jl")

opt = OptionsNLEQ("RTOL" => 1e-6,"PR_ERR" => 2);

function f(x)
    return x;
end

function Jac(x)
    return x;
end

x = ones(4);
xScal = ones(4);

retCode = nleq1(f,x,xScal,opt);

println("$retCode");
