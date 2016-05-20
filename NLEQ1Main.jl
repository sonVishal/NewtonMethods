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
                dxQa = dxQ
                if !qOrdi
                    fcNumP = sum((dx./xw).^2)
                    th = fc - 1.0
                    fcDnm = sum(((dxQa+th*dx)./xw).^2)
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
                    fcNmp2 = sum((dxQa./xw).^2)
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
            alfa1 = sum(dx.*dxQ./xw.^2)
            alfa2 = sum(dx.^2./xw.^2)
            alfa = alfa1/alfa2
            beta = 1.0 - alfa
            t1 = (dxQ+(fca-one)*alfa*dx)/beta
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
        # evaluation of (scaled) standard level function dlevf
        # and computation of ordinary Newton corrections dx[n]
        if !qSimpl
            (dx,conv,sumX,dLevF) = n1lvls(n,dx,t1,xw,f,newt == 0)
        else
            (dx,conv,sumX,dLevF) = n1lvls(n,dx,t1,xwa,f,newt == 0)
        end
        stats[STATS_SUMX]   = sumX
        stats[STATS_DLEVF]  = dLevF
        xa[:]    = x
        sumXa[:] = sumX
        dLevXa   = sqrt(sumXa/n)
        conva    = conv
        dxANrm   = wnorm(n,dx,xw)
        push!(sumXall,dLevXa)
        push!(dLevFall,dLevF)

        # ----------------------------------------------------------------------
        # 3.3 A - priori estimate of damping factor fc
        if nIter != 0 && nonLin != 1 && newt == 0 && !qOrdi
            # 3.3.1 Computation of the denominatior of a-priori estimate
            fcDnm = sum(((dx-dxQa)./xw).^2)*sumX
            # ------------------------------------------------------------------
            # 3.3.2 New damping factor
            if fcDnm > fcNumP*fcMin2 || (nonLin == 4 && fcA^2*fcNumP < 4.0*fcDnm)
                dMyPri = fcA*sqrt(fcNumP/fcDnm)
                fcPri  = min(dMyPri,1.0)
                if nonLin == 4
                    fcPri = min(0.5*dMyPri,1.0)
                end
            else
                fcPri = 1.0

                dMyPri = -1.0
            end

            if printIterationMonitor >= 5
                write(printIO,"\n",
                @sprintf("+++++  apriori estimate  +++++\n"),
                @sprintf(" fcPri  = %18.10e\n",fcPri),
                @sprintf(" fc     = %18.10e\n",fc),
                @sprintf(" fcA    = %18.10e\n",fcA),
                @sprintf(" dMyPri = %18.10e\n",dMyPri),
                @sprintf(" fcNumP = %18.10e\n",fcNumP),
                @sprintf(" fcDnm  = %18.10e\n",fcDnm),
                @sprintf("++++++++++++++++++++++++++++++\n"))
            end

            fc = max(fcPri,fcMin)
            if qBDamp
                fcbh = fcA*fcBand;
                if fc > fcbh
                    fc = fcbh
                    if printIterationMonitor >= 4
                        write(printIO, "*** Increase rest. act. (a priori)\n")
                    end
                end
                fcbh = fcA/fcBand
                if fc < fcbh
                    fc = fcbh
                    if printIterationMonitor >= 4
                        write(printIO, "*** Decrease rest. act. (a priori)\n")
                    end
                end
            end
        end

        if iOrMon >= 2
            sumxa2 = sumxa1
            sumxa1 = sumxa0
            sumxa0 = dLevXa
            if sumxa0 == 0.0
                sumxa0 = small
            end
            # Check convergence rate (linear or superlinear)
            # iconv : Convergence idicator
            #           = 0: No convergence indicated yet
            #           = 1: Damping factor is 1.0e0
            #           = 2: Superlinear convergence detected (alpha >=1.2)
            #           = 3: Quadratic convergence detected (alpha > 1.8)
            fcMon = min(fc,fcMon)
            if fcMon < 1.0
                iConv  = 0
                alphaE = 0.0
            end
            if fcMon == 1.0 && iConv == 0
                iConv = 1
            end
            if nIter >= 1
                cLin1 = cLin0
                cLin0 = sumxa0/sumxa1
            end
            if iConv >= 1 && nIter >= 2
                alphaK = alphaE
                alphaE = 0.0
                if cLin1 <= 0.95
                    alphaE = log(cLin0)/log(cLin1)
                end
                if alphaK != 0.0
                    alphaK = 0.5*(alphaE+alphaK)
                end
                alphaA = min(alphaK,alphaE)
                cAlphaK = cAlpha
                cAlpha = 0.0
                if alphaE != 0.0
                    cAlpha = sumxa1/sumxa2^alphaE
                end
                sumXte = sqrt(cAlpha*cAlphaK)*sumxa1^alphaK-sumxa0
                if alphaA >= 1.2 && iConv == 1
                    iConv = 2
                end
                if alphaA > 1.8
                    iConv = 3
                end
                if printIterationMonitor >= 4
                    write(printIO,"\n",
                    @sprintf(" ** iConv: %1i",iConv),
                    @sprintf("  alpha:       %9.2e",alphaE),
                    @sprintf("  const-alpha: %9.2e",cAlpha),
                    @sprintf("  const-lin:   %9.2e **\n",cLin0),
                    @sprintf(" ** alpha-post: %9.2e",alphaK),
                    @sprintf("  check:       %9.2e",sumXte),
                    @sprintf("                    **\n"))
                end
                if iConv >= 2 && alphaA < 0.9
                    if iOrMon == 3
                        retCode = 4
                        break;
                    else
                        qmStop = 1
                    end
                end
            end
        end
        fcMon = fc
        # ----------------------------------------------------------------------
        # 3.4 Save natural level for later computations of corrector
        # and print iterate
        fcNumK = sumX
        if printIterationMonitor >= 2
            n1prv1(dLevF,dLevXa,fcKeep,nIter,newt,printIterationMonitor,printIO,qMixIO)
        end
        nRed    = 0
        qNext   = 0
        qRep    = 0
        qRed    = 0
        iCnv    = 0

        # Damping-factor reduction loop
        # ======================================================================
        while qRed
            # ------------------------------------------------------------------
            # 3.5 Preliminary new iterate
            x = xa + dx*fc
            push!(fcAll,fc)
            # ------------------------------------------------------------------
            # 3.5.2 Exit, if problem is specified as being linear
            if nonLin == 1
                retCode = 0
                break
            end
            #-------------------------------------------------------------------
            # 3.6.1 Computation of the residual vector
            try
                fcn(x,f)
            catch err
                retCode = 82
                throw(EvaluationError(fcn,err))
            end
            nFcn += 1
            # TODO: Understand what is happening here
            # and handle the failure properly
            # What does iFail = 1 and iFail = 2 mean??
            if iFail == 1 || iFail == 2
                if qOrdi
                    retCode = 82
                    break
                end
                if iFail == 1
                    fcRedu = 0.5
                else
                    fcRedu = f[1]

                    if fcRedu <= 0 || fcRedu >= 1
                        retCode = 82
                        break
                    end
                end
                if printIterationMonitor >= 2
                    write(printIO,
                    @sprintf("        %2i",nIter),
                    @sprintf(" %s could not be evaluated     %7.5f    %2i\n",fcn,fc,newt))
                end
                fch = fc
                fc  = fcRedu*fc
                if fch > fcMin
                    fc = max(fc,fcMin)
                end
                if qBDamp
                    fcbh = fch/fcBand
                    if fc < fcbh
                        fc = fcbh
                        if printIterationMonitor >= 4
                            write(printIO," *** Decrease rest. act. (fcn redu.) ***\n")
                        end
                    end
                end
                if fc < fcMin
                    retCode = 3
                    break
                end
            else
                if qOrdi
                    # ----------------------------------------------------------
                    # 3.6.2 Convergence test for ordinary Newton iteration
                    push!(tolAll,dxANrm)
                    if dxANrm <= rTol
                        retCode = 0
                        iCnv    = 1
                        break
                    end
                else
                    t1 = f.*fw
                    # ------------------------------------------------------
                    # 3.6.3 Solution of linear (n,n) system
                    (t1,iFail) = n1solv(n,m1,ml,mu,l,u,p,t1,opt)
                    if iFail != 0
                        retCode = 81
                        break
                    end
                    if newt > 0
                        dxQ = t1.*xwa
                        for iLoop = 1:newt
                            sum1 = sum((dxQ.*dxSave[1:n,iLoop])./xw.^2)
                            sum2 = sum((dxSave[1:n,iLoop]./xw).^2)
                            beta = sum1/sum2
                            dxQ = dxQ + beta*dxSave[1:n,iLoop+1]
                            t1 = dxQ./xw
                        end
                    end
                    # ------------------------------------------------------
                    # 3.6.4 Evaluation of scaled natural level function
                    #       sumX
                    #       scaled maximum error norm conv and evaluation
                    #       of (scaled) standard level function dLevFn
                    if !qSimpl
                        (dxQ,conv,sumX,dLevFn) =
                        n1lvls(n,dxQ,t1,xw,f,printIterationMonitor,newt==0)
                    else
                        (dxQ,conv,sumX,dLevFn) =
                        n1lvls(n,dxQ,t1,xwa,f,printIterationMonitor,newt==0)
                    end
                    push!(sumXQall,sqrt(sumX/n))
                    dxNrm = wnorm(n,dxQ,xw)
                    # ------------------------------------------------------
                    # 3.6.5 Convergence test
                    push!(tolAll,dxNrm)
                    if dxNrm <= rTol && dxANrm <= rSmall && fc == 1.0
                        retCode = 0
                        iCnv = 1
                        break
                    end

                    fcA = fc
                    # ------------------------------------------------------
                    # 3.6.6 Evaluation of reduced damping factor
                    th = fcA - 1.0
                    fcDnm = sum(((dxQ+th*dx)./xw).^2)
                    if fcDnm != 0.0
                        dMyCor = fcA*fcA*0.5*sqrt(fcNumK/fcDnm)
                    else
                        dMyCor = 1.0e+35
                    end
                    if nonLin <= 3
                        fcCor = min(1.0,dMyCor)
                    else
                        fcCor = min(1.0,0.5*dMyCor)
                    end

                    if printIterationMonitor >= 5
                        write(printIO,
                        @sprintf(" +++ corrector computation +++\n"),
                        @sprintf("  fcCor    = %18.10e\n",fcCor),
                        @sprintf("  fc       = %18.10e\n",fc),
                        @sprintf("  dMyCor   = %18.10e\n",dMyCor),
                        @sprintf("  fcNumK   = %18.10e\n",fcNumK),
                        @sprintf("  fcDnm    = %18.10e\n",fcDnm),
                        @sprintf("  fcA      = %18.10e\n",fcA),
                        @sprintf("++++++++++++++++++++++++++++++\n"))
                    end
                end
                # ----------------------------------------------------------
                # 3.7 Natural monotonicity test
                if sumX > sumXa && !qOrdi
                    # ------------------------------------------------------
                    # 3.8 Output of iterate
                    if printIterationMonitor >= 3
                        n1prv2(dLevFn,sqrt(sumX/n),fc,niter,printIterationMonitor,
                        printIO,qMixIO,"*")
                    end
                    if qmStop
                        retCode = 4
                        break
                    end
                    fch = min(fcCor,0.5*fc)
                    if fc > fcMin
                        fc = max(fch,fcMin)
                    else
                        fc = fch
                    end
                    if qBDamp
                        fcbh = fcA/fcBand
                        if fc < fcbh
                            fc = fcbh
                            if printIterationMonitor >= 4
                                write(printIO,
                                " *** Decrese rest. sct. (a posteriori) ***\n")
                            end
                        end
                    end
                    fcMon = fc

                    if printIterationMonitor >= 5
                        write(printIO,
                        " +++ corrector setting 1 +++\n",
                        @sprintf("fc    = %18.10e\n",fc),
                        " +++++++++++++++++++++++++++\n")
                    end

                    qRep = 1
                    nCorr += 1
                    nRed += 1
                    # ------------------------------------------------------
                    # 3.10 If damping factor is too small:
                    #      Refreash Jacobian, if current Jacobian was computed
                    #      by a Rank-1 update, else fail and exit
                    qJcRfr = fc < fcMin || (newt > 0 && nRed > 1)

                    if qJcRfr && newt == 0
                        retCode = 3
                        break
                    end
                else
                    if !qOrdi && !qRep && fcCor > sigma2*fc
                        if printIterationMonitor >= 3
                            n1prv2(dLevFn,sqrt(sumX/n),fc,nIter,
                            printIterationMonitor,printIO,qMixIO,"+")
                        end
                        fc = fcCor

                        if printIterationMonitor >= 5
                            write(printIO,
                            " +++ corrector setting 2 +++\n",
                            @sprintf("fc    = %18.10e\n",fc),
                            " +++++++++++++++++++++++++++\n")
                        end

                        qRep = 1
                    else
                        qNext = 1
                    end
                end
            end
            qRed = !(qNext||qJcRfr)
        end
        # End of damping factor reduction loop
        if nonLin == 1 || iCnv == 1 || (retCode != 0 && retCode != -1)
            break
        end
        # ======================================================================
        if qJcRfr
            # ------------------------------------------------------------------
            # 3.11 Restore former values for repeating iteration step
            nRejR1 += 1
            x[:] = xa
            f[:] = fa
            if printIterationMonitor >= 2
                write(printIO,
                @sprintf("        %2i not accepted damping factor %7.5f",nIter,fc),
                @sprintf("    %2i\n",newt))
            end
            fc  = fcKeep
            fcA = fck2
            if nIter == 0
                fc = fcMin
            end
            qGenJ = 1
        else
            # ------------------------------------------------------------------
            # 4 Preparations to start the following iteration step
            # ------------------------------------------------------------------
            # Print values
            if printIterationMonitor >= 3 && !qOrdi
                n1prv2(dLevFn,sqrt(sumX/n),fc,nIter+1,printIterationMonitor,
                printIO,qMixIO,"*")
            end
            # Print the natural level of the current iterate and
            # return it in one-step mode
            sumXs = sumX
            sumX = sumXa
            if printSolution >= 2 && nIter != 0
                n1sout(n,xa,2,opt,wk,printSolution,printIO)
            elseif printSolution >= 1 && nIter == 0
                n1sout(n,xa,1,opt,wk,printSolution,printIO)
            end
            nIter += 1
            stats[STATS_NITER] = nIter
            push!(xIter,x)
            dLevF = dLevFn
            if nIter >= nItmax
                retCode = 2
                break
            end
            fcKeep = fc
            # ------------------------------------------------------------------
            # 4.2 Return if in one-step mode
            if mode == 1
                qSucc = 1
                return
            end
        end
    end
    # End repeat
    # End of main iteration loop
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # 9 Exits
    # --------------------------------------------------------------------------
    # 9.1 Solution exit
end
