using NewtonMethods
include("ChebyQuad.jl")
fSol = open("nleq1.dat","w")
fRest = open("nleq1.out","w")

refSol = Dict{Int64,Vector}()

refSol[2] = [0.21132486540517281, 0.78867513459482719]
refSol[3] = [0.14644660863174364, 0.5, 0.85355339136825636]
refSol[4] = [0.10267276384825177, 0.40620376278203324, 0.59379623721796682,
                0.89732723615174825]
refSol[5] = [0.083751256499240576, 0.31272929522379189, 0.5, 0.68727070477620811,
                0.91624874350075947]
refSol[6] = [0.066876587754941208, 0.28874077956671557, 0.36668218955839393,
                0.63331781044160607, 0.71125922043328438, 0.93312341224505879]
refSol[7] = [0.058069149620871681, 0.23517161236154796, 0.33804409473596703, 0.5,
                0.66195590526403303, 0.76482838763845207, 0.94193085037912838]
refSol[8] = [0.047550008453679384, 0.2380258672311448, 0.24314594839706133,
                0.50202780086466359, 0.49797219913533641, 0.75685405160281116,
                0.76197413276898263, 0.95244999154632048]
refSol[9] = [0.044205346135779103, 0.19949067230774983, 0.23561910847320636,
                0.41604690789257154, 0.50000000000000011, 0.58395309210742818,
                0.76438089152679367, 0.80050932769225003, 0.95579465386422091]

# for dim = 2:9
dim = 2
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
                  OPT_NITMAX            => 10,
                  OPT_RTOL              => 1e-5)

x0    = collect(1:dim)./n1
xScal = zeros(x0)

retCode = -1

i = 1

while retCode == -1
    (x0, _, retCode) = nleq1(chebyQuad, x0, xScal, opt)
    write(fRest, @sprintf("Returned from call %4i of NLEQ1\n",i))
    flush(fSol)
    flush(fRest)
    i += 1
end
err = norm(x0-refSol[dim],Inf)/norm(refSol[dim],Inf)
println("Dimension $dim, relative error in Inf norm = $err")
# end

close(fSol)
close(fRest)
