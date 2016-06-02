include("NLEQ1.jl")
function fcn(x,f)
    f = zeros(x)
    f[1] = x[1]
    f[2] = x[2]
    return nothing
end

opt = OptionsNLEQ(OPT_PRINTWARNING => 1, OPT_PRINTITERATION => 5)

x = ones(2);
xScal = ones(2);

(sol, stats, retCode) = nleq1(fcn,x,xScal,opt);

println("retCode = $retCode","\n");
println("Options = $opt","\n");
println("SOlution = $sol");

# throw(EvaluationError(Jac));
