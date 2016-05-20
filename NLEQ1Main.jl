function n1int(n, fcn, jac, x, xScal, rTol, nItmax, nonLin, opt, retCode, wk,
    m1, m2, nBroy, xIter, sumXall, dLevFall, sumXQall, tolAll, fcAll,
    a, dxSave, dx, dxQ, xa, xwa, f, fa, eta, xw, fw, dxQa, sumxa0, sumxa1, fcmon,
    fcStart, fcMin, sigma, sigma2, fcA, fcKeep, fcPri, dMyCor, conv, sumXs,
    dLevF, mStor, printWarning, printIterationMonitor, printSolution, printIO,
    nIter, nCorr, nFcn, nJac, nRejR1, newt, iConv, qBDamp, stats)

    epMach  = getMachineConstants(3)
    small   = getMachineConstants(6)
    n       = length(x)
    # Begin
    # --------------------------------------------------------------------------
    # 1 Initialization
    # --------------------------------------------------------------------------
    # 1.1 Control variables
    qSucc       = opt.options[OPT_QSUCC]        == 1
    qScale      = opt.options[OPT_NOROWSCAL]    != 1
    qOrdi       = opt.options[OPT_QORDI]        == 1
    qSimpl      = opt.options[OPT_QSIMPL]       == 1
    qRank1      = opt.options[OPT_QRANK1]       == 1
    iOrMon      = getOption!(opt, OPT_IORMON, 0)
    if iOrMon == 0
        iOrMon = 2
    end
    iScal       = getOption!(opt, OPT_ISCAL, 0)
    mode        = getOption!(opt, OPT_MODE,  0)
    jacGen      = opt.options[OPT_JACGEN]
    qMixIO      = printIterationMonitor != 0 && printSolution != 0
    qLU         = !qSimpl
    # --------------------------------------------------------------------------
    # 1.2 Derived dimensional parameters
    if mStor == 0
        ml = 0
        mu = 0
    elseif mStor == 1
        ml = m1 - m2
        mu = m2 - 1 - ml
    end
    # --------------------------------------------------------------------------
    # 1.3 Derived internal parameters
    fcMin2  = fcMin*fcMin
    rSmall  = sqrt(10.0*rTol)
    # --------------------------------------------------------------------------
    # 1.4 Adaptation of input parameters, if necessary
    if fc < fcMin
        fc = fcMin
    end
    if fc > 1.0
        fc = 1.0
    end
    # --------------------------------------------------------------------------
    # 1.5 Initial preparations
    qJcRfr              = false
    qLInit              = false
    qIter               = true
    iFail               = 0
    fcBand              = 0.0
    if qBDamp == 1
        fcBand = fcband
    end
    # --------------------------------------------------------------------------
    # 1.5.1 Numerical differentiation related initializations
    if jacGen == 2
        aJdel = getOption!(opt,OPT_AJDEL,0.0)
        if aJdel <= small
            aJdel = sqrt(epMach*10.0)
        end
        aJmin = getOption!(opt,OPT_AJMIN,0.0)
    elseif jacGen == 3
        etaDif = getOption!(opt,OPT_ETADIF,0.0)
        if etaDif <= small
            etaDif = 1.0e-6
        end
        etaIni = getOption!(opt,OPT_ETAINI,0.0)
        if etaIni <= small
            etaIni = 1.0e-6
        end
        epDiff = sqrt(epMach*10.0)
        etaMax = sqrt(epDiff)
        etaMin = epDiff*etaMax
    end
    # --------------------------------------------------------------------------
    # 1.5.2 Miscellaneous preparations of the first iteration step
    if !qSucc
        push!(xIter,x)

        nIter  = 0
        stats[STATS_NITER] = nIter;
        nCorr  = 0
        nRejR1 = 0
        nFcn   = 0
        nJac   = 0
        iConv  = 0
        conv   = 0.0

        fcKeep  = fc
        fcA     = fc
        fcPri   = fc
        fcMon   = fc
        sumxa0  = 0.0
        sumxa1  = 0.0

        qGenJ   = true
        qIniSc  = true
        fcK2    = fc

        if jacGen == 3
            eta = etaIni*ones(n)
        end

        xa[:] = x[:]

        alphaE  = 0.0
        cAlpha  = 0.0
        cLin0   = 0.0
        qMStop  = false
        # ----------------------------------------------------------------------
        # 1.6 Print monitor header
        if printIterationMonitor >= 2 && !qMixIO
            write(printIO,
            "\n",
            "*******************************************************************",
            "\n",
            "        It       Normf           Normx         Damp.Fct.   New\n")
        end
        # ----------------------------------------------------------------------
        # 1.7 Startup step
        # ----------------------------------------------------------------------
        # 1.7.1 Computation of residual vector
        try
            fcn(x,f)
        catch err
            retCode = 82
            # qIter   = false
            throw(EvaluationError(fcn,err))
        end
        nFcn += 1
        if length(f) ~= n
            retCode = 22
            # qIter   = false
            throw(DimensionMismatch("Dimension of the function output does not
            match the input dimension. Please check the function \"$fcn\" again."))
        end
    else
        qIniSc = false
    end
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Main iteration loop

    # Repeat
    while qIter
        # ----------------------------------------------------------------------
        # 2 Startup of iteration step
        if !qJcRfr
            # ------------------------------------------------------------------
            # 2.1 Scaling of variables x(n)
            xw = n1scal(n , x, xa, xScal, iScal, qIniSc, opt)
            qIniSc = false
            if nIter != 0
                # --------------------------------------------------------------
                # 2.2 Aposteriori estimate of damping factor
                dxqa = dxq
                if !qOrdi
                    fcNumP = sum((dx./xw).^2)
                    th = fc - 1.0
                    fcDnm = sum(((dxqa+th*dx)./xw).^2)
                    # ----------------------------------------------------------
                    # 2.2.2 Decision criterion for Jacobian update technique
                    # qGenJ == true   numerical differentation,
                    # qGenJ == false  rank1 updating
                    qGenJ = true
                    if fc == fcPri
                        qGenJ = fc < 1.0 || fcA < 1.0 || dMyCor <= fc*sigma ||
                        !qRank1 || newt + 2 > nBroy

                        fcA = fc
                    else
                        dMyCor = fcA*fcA*0.5*sqrt(fcNump/fcDnm)
                        if nonLin <= 3
                            fcCor = min(1.0,dMyCor)
                        else
                            fcCor = min(1.0,0.5*dMyCor)
                        end
                        fcA = max(min(fc,fcCor),fcMin)

                        if printIterationMonitor >= 5
                            write(printIO,"\n",
                            @sprintf("+++  aposteriori estimate  +++\n"),
                            @sprintf(" fcCor  = %18.10e\n",fcCor),
                            @sprintf(" fc     = %18.10e\n",fc),
                            @sprintf(" dMyCor = %18.10e\n",dMyCor),
                            @sprintf(" fcNumP = %18.10e\n",fcNumP),
                            @sprintf(" fcDnm  = %18.10e\n",fcDnm),
                            @sprintf("++++++++++++++++++++++++++++++\n"))
                        end
                    end
                    fck2 = fcA
                    # ----------------------------------------------------------
                    # 2.2.1 Computation of the numerator of damping
                    # factor predictor
                    fcNmp2 = sum((dxqa./xw).^2)
                    fcNump = fcNump*fcNmp2
                end
            end
        end
        qJcRfr = false
        # ----------------------------------------------------------------------
        # 2.3 Jacobian matrix
        # ----------------------------------------------------------------------
        # 2.3.1 Jacobian generation by routine jac or difference approximation
        # if qGenJ == true
        # - or -
        # Rank-1 update of Jacobian if qGenJ == false
        if qGenJ && (!qSimpl || nIter == 0)
            newt = 0
            if jacGen == 1
                try
                    jac(x,a)
                catch err
                    throw(EvaluationError(jac,err))
                end
            else
                if mStor == 0
                    if jacGen == 3
                        #TODO: Jacobian with num diff feedback
                    end
                    if jacGen == 2
                        #TODO: Jacobian with num diff without feedback
                    end
                elseif mStor == 1
                    if jacGen == 3
                        #TODO: banded num diff with feedback
                    end
                    if jacGen == 2
                        #TODO: banded num diff without feedback
                    end
                end
            end
            nJac += 1

            #TODO: break if jacobian computation fails
            # Use try catch method
            if jacGen == 1 && iFail < 0
                retCode = 83
                break
            end

            if jacGen != 1 && iFail != 0
                retCode = 82
                break
            end

        elseif ~qSimpl
            newt += 1
        end

        if newt == 0 && (qLU || nIter == 0)
            # ------------------------------------------------------------------
            # 2.3.2.1 Save scaling values
            xwa = xw[1:n]
            if issparse(a)
                nza = nnz(a)
                (row,col) = findn(a);
                for k = 1:nza
                    a[row[k],col[k]] = -a[row[k],col[k]]*xw[col[k]]
                end
            else
                # --------------------------------------------------------------
                # 2.3.2.2 Prepare Jacobian for use by band-solver
                if mStor == 1
                    for l1 = 1:n
                        for l2 = m2:-1:1
                            a[l2+ml,l1] = a[l2,l1]
                        end
                    end
                end
                # --------------------------------------------------------------
                # 2.4 Prepare solution of the linear system
                # --------------------------------------------------------------
                # 2.4.1 Internal column scaling of matrix A
                if mStor == 0
                    for k = 1:n
                        a[1:n,k] = -a[1:n,k]*xw[k]
                    end
                elseif mStro == 1
                    for k = 1:n
                        l2 = max(1+m2-k,ml+1)
                        l3 = max(n+m2-k,m1)
                        a[l2:l3,k] = -a[l2:l3,k]*xw[k]
                    end
                end
            end
            # ------------------------------------------------------------------
            # 2.4.2 Row scaling of matrix A
            if qScale
                if mStor == 0
                    (a,fw) = n1scrf(n,n,a)
                elseif mStro == 1
                    (a,fw) = n1scrb(n,m1,ml,mu,a)
                end
            else
                fw = ones(n);
            end
        end
        # ----------------------------------------------------------------------
        # 2.4.3 Save and scale values of F(n)
        fa[:] = f;
        t1 = f.*fw;
        # ----------------------------------------------------------------------
        # 3 Central part of iteration step
        # ----------------------------------------------------------------------
        # 3.1 Solution of the linear system
        # ----------------------------------------------------------------------
        # 3.1.1 Decomposition of (n,n) matrix A
        if newt == 0 && (qLU || nIter == 0)
            (l,u,p,iFail) = n1fact(n,m1,ml,mu,a,opt);
            if iFail != 0
                if iFail == 1
                    retCode = 1
                else
                    retCode = 80
                end
                break
            end
        end
        qLInit = 1;
        # ----------------------------------------------------------------------
        # 3.1.2 Solution of (n,n) system
        if newt == 0
            (t1,iFail) = n1solv(n,m1,ml,mu,l,u,p,t1,opt);
            if iFail != 0
                retCode = 81
                break
            end
        else
            alfa1 = sum(dx.*dxq./xw.^2)
            alfa2 = sum(dx.^2./xw.^2)
            alfa = alfa1/alfa2
            beta = 1.0 - alfa
            t1 = (dxq+(fca-one)*alfa*dx)/beta
            if newt == 1
                dxSave[1:n,1] = dx
            end
            dxSave[1:n,newt+1] = t1
            dx = t1
            t1 = t1./xw
        end
        # ----------------------------------------------------------------------
        # 3.2 Evaluation of scaled natural level function sumX
        # scaled maximum error norm conv
        # evaluation of (scaled) standard level function
        # and computation of ordinary Newton corrections dx[n]
        if ~qSimpl
            (dx,conv,sumx,dlevf) = n1lvls(n,dx,t1,xw,f,newt == 0)
        else
            (dx,conv,sumx,dlevf) = n1lvls(n,dx,t1,xwa,f,newt == 0)
        end



    end
