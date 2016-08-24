function chebyQuad(f,x)
    n = length(x)
    f[:] = zeros(n);
    for i = 2:2:n
        f[i] = n/(i*i-1)
    end
    for l = 1:n
        factt = 4.0*x[l] - 2.0
        ti2 = 1.0
        ti1 = 0.5*factt
        f[1] = ti1 + f[1]
        for i = 2:n
            ti = factt*ti1-ti2
            f[i] = ti + f[i]
            ti2 = ti1
            ti1 = ti
        end
    end
    return nothing
end

function chebyQuadJac(x,J)
    n = length(x)
    for j = 1:n
        factt = 4.0*x[j] - 2.0
        ti2 = 1.0
        ti1 = 0.5*factt
        tabli2 = 0.0
        tabli1 = 2.0
        J[1,j] = tabli1
        for i = 2:n
            ti = factt*ti1-ti2
            tabli = 4.0*ti1+factt*tabli1-tabli2
            J[i,j] = tabli
            ti2 = ti1
            ti1 = ti
            tabli2 = tabli1
            tabli1 = tabli
        end
    end
    return nothing
end

function testNLEQ1()
    refSol = Dict{Int32,Any}();
    refSol[2] = [0.211324865405173;
                0.788675134594827]

    refSol[3] = [0.146446608631744;
                0.5;
                0.853553391368256]

    refSol[4] = [0.102672763848252;
                0.406203762782033;
                0.593796237217967;
                0.897327236151748]

    refSol[5] = [0.0837512564992406;
                0.312729295223792;
                0.5;
                0.687270704776208;
                0.916248743500759]

    refSol[6] = [0.0668765877549412;
                0.288740779566716;
                0.366682189558394;
                0.633317810441606;
                0.711259220433284;
                0.933123412245059]

    refSol[7] = [0.0580691496208717;
                0.235171612361548;
                0.338044094735967;
                0.5;
                0.661955905264033;
                0.764828387638452;
                0.941930850379128]

    refSol[8] = [0.0475500084536796;
                0.238025867229826;
                0.243145948398381;
                0.502027800864952;
                0.497972199135048;
                0.756854051601463;
                0.76197413277033;
                0.95244999154632]

    refSol[9] = [0.0442053461357791;
                0.19949067230775;
                0.235619108473206;
                0.416046907892572;
                0.5;
                0.583953092107428;
                0.764380891526794;
                0.80050932769225;
                0.955794653864221]

    testResult = true;
    for dim = 2:9
        n1 = dim + 1

        # Initialize the options
        opt = OptionsNLEQ(OPT_MODE              => 1,
                          OPT_JACGEN            => 1,
                          OPT_JACFCN            => chebyQuadJac,
                          OPT_MSTOR             => 0,
                          OPT_NOROWSCAL         => 0,
                          OPT_NITMAX            => 10,
                          OPT_RTOL              => 1e-5)

        x0    = collect(1:dim)./n1
        xScal = zeros(x0)

        retCode = -1
        stats   = []

        while retCode == -1
            (x0, stats, retCode) = nleq1(chebyQuad,x0,xScal,opt)
        end

        relNormDiff = norm(x0-refSol[dim])/norm(refSol[dim])

        testResult &= relNormDiff <= 1e-12
    end
    return testResult
end

testNLEQ1()
