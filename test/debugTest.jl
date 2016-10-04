using NewtonMethods
include("ChebyQuad.jl")
fSol = open("nleq2.dat","w")
fRest = open("nleq2.out","w")

refSol = Dict{Int64,Vector}()

refSol[2] = [0.2113248654, 0.7886751346]
refSol[3] = [0.1464466086, 0.5000000000, 0.8535533914]
refSol[4] = [0.1026727638, 0.4062037628, 0.5937962372, 0.8973272362]
refSol[5] = [0.8375125650e-1, 0.3127292952, 0.5000000000, 0.6872707048,
                0.9162487435]
refSol[6] = [0.6687658775e-1, 0.2887407796, 0.3666821896, 0.6333178104,
                0.7112592204, 0.9331234122]
refSol[7] = [0.5806914962e-1, 0.2351716124, 0.3380440947, 0.5000000000,
                0.6619559053, 0.7648283876, 0.9419308504]
refSol[8] = [0.5473028465e-1, 0.2096506534, 0.2670502974, 0.5000286281,
                0.5000338603, 0.7329693255, 0.7903332923, 0.9451666568]
refSol[9] = [0.4420534614e-1, 0.1994906723, 0.2356191085, 0.4160469079,
                0.5000000000, 0.5839530921, 0.7643808915, 0.8005093277,
                0.9557946539]

for dim = 2:9
    n1 = dim + 1

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

    i = 1

    while retCode == -1
        (x0, _, retCode) = nleq2(chebyQuad, x0, xScal, opt)
        write(fRest, @sprintf("Returned from call %4i of NLEQ2\n",i))
        flush(fSol)
        flush(fRest)
        i += 1
    end
    err = norm(x0-refSol[dim],Inf)/norm(refSol[dim],Inf)
    println("Relative error in Inf norm = $err")
end
close(fSol)
close(fRest)
