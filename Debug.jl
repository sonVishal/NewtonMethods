include("Error.jl")

include("Helper.jl")

include("Jacobian.jl")

include("Options.jl")

include("CheckOptionsNLEQ1.jl")

include("NLEQ1Main.jl")

include("NLEQ1.jl")

println("Done including all the files")

function chebyquad(n::Integer)
    name = "chebyquad"
    const tk = 1/n

    function f!(x, fvec)
        fill!(fvec, 0)
        for j = 1:n
            temp1 = 1.0
            temp2 = 2x[j]-1
            temp = 2temp2
            for i = 1:n
                fvec[i] += temp2
                ti = temp*temp2 - temp1
                temp1 = temp2
                temp2 = ti
            end
        end
        iev = -1.0
        for k = 1:n
            fvec[k] *= tk
            if iev > 0
                fvec[k] += 1/(k^2-1)
            end
            iev = -iev
        end
    end

    function g!(x, fjac)
        for j = 1:n
            temp1 = 1.
            temp2 = 2x[j] - 1
            temp = 2*temp2
            temp3 = 0.0
            temp4 = 2.0
            for k = 1:n
                fjac[k,j] = tk*temp4
                ti = 4*temp2 + temp*temp4 - temp3
                temp3 = temp4
                temp4 = ti
                ti = temp*temp2 - temp1
                temp1 = temp2
                temp2 = ti
            end
        end
    end
    x = collect(1:n)/(n+1)
    return name, f!, g!, x
end

function test(x,f)
    f[:] = x;
end

#opt = OptionsNLEQ(OPT_PRINTWARNING => 1, OPT_PRINTITERATION => 5);
optTest = OptionsNLEQ(OPT_RTOL => 1e-3);

#(name,fcn,jac,x) = chebyquad(2);

x = ones(3);
xScal = ones(3);
println("Done setting up the problem")
println("Calling the function")
(sol, stats, retCode) = nleq1(test,x,xScal,optTest);
println("Done calling the function")
println("retCode = $retCode","\n");
println("Options = $optTest","\n");
println("Answer = $sol")
