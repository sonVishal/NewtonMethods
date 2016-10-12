function testNLEQ1()
    refSol = Dict{Int32,Any}();
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
    refSol[8] = [0.0475500084536796, 0.238025867229826, 0.243145948398381,
                    0.502027800864952, 0.497972199135048, 0.756854051601463,
                    0.76197413277033,0.95244999154632]
    refSol[9] = [0.044205346135779103, 0.19949067230774983, 0.23561910847320636,
                    0.41604690789257154, 0.50000000000000011, 0.58395309210742818,
                    0.76438089152679367, 0.80050932769225003, 0.95579465386422091]

    testResult = false;
    for dim = 2:9
        n1 = dim + 1

        # Initialize the options
        opt = OptionsNLEQ(OPT_JACGEN            => 1,
                          OPT_JACFCN            => chebyQuadJac,
                          OPT_MSTOR             => 0,
                          OPT_NOROWSCAL         => 0,
                          OPT_NITMAX            => 10,
                          OPT_RTOL              => 1e-5)

        x0    = collect(1:dim)./n1
        xScal = zeros(x0)

        (x0, _, retCode) = nleq1(chebyQuad,x0,xScal,opt)

        relNormDiff = norm(x0-refSol[dim],Inf)/norm(refSol[dim],Inf)

        # Here we check < 1e-12 due to the fact that we use the Julia LU Decomposition
        # rather than the same code for factorization and solution as with NLEQ1
        if dim == 2
            testResult = relNormDiff < 1e-12
        else
            testResult &= relNormDiff < 1e-12
        end
    end
    return testResult
end

testNLEQ1()
