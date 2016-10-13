using NewtonMethods
include("ChebyQuad.jl")
fSol = open("nleq1.dat","w")
fRest = open("nleq1.out","w")

refSol = Dict{Int64,Vector}()

refSol[2] = [0.21132486540517284, 0.78867513459482719];

refSol[3] = [0.14644660863174364, 0.5, 0.85355339136825636];

refSol[4] = [0.10267276384825176, 0.40620376278203324, 0.5937962372179667,
                0.89732723615174825];

refSol[5] = [0.083751256499240576, 0.31272929522379189, 0.5, 0.68727070477620811,
                0.91624874350075947];

refSol[6] = [0.066876587754941222, 0.28874077956671557, 0.36668218955839399,
                0.63331781044160607, 0.71125922043328438, 0.93312341224505879];

refSol[7] = [0.058069149620871653, 0.23517161236154782, 0.33804409473596703, 0.5,
                0.66195590526403303, 0.76482838763845207, 0.94193085037912838];

refSol[8] = [0.0475500084536796, 0.238025867229826, 0.243145948398381,
                0.502027800864952, 0.497972199135048, 0.756854051601463,
                0.76197413277033,0.95244999154632];

refSol[9] = [0.044205346135779117, 0.19949067230774983, 0.23561910847320633,
                0.41604690789257159, 0.5, 0.5839530921074283, 0.76438089152679367,
                0.80050932769225014, 0.95579465386422091];

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
msg = "Solution rejected"
if err < 1e-12
    msg = "Solution accepted"
end
println(msg)

close(fSol)
close(fRest)
