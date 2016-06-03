#include("Helper.jl")
#include("Jacobian.jl")
function n1int(n, fcn, x, xScal, rTol, nItmax, nonLin, opt, retCode, wk,
    m1, m2, nBroy, xIter, sumXall, dLevFall, sumXQall, tolAll, fcAll,
    a, dxSave, dx, dxQ, xa, xwa, f, fa, eta, xw, fw, dxQa, sumxa0, sumxa1, fcMon,
    fc, fcMin, sigma, sigma2, fcA, fcKeep, fcPri, dMyCor, conv, sumXs,
    dLevF, mStor, printWarning, mPr, mPrSol, printIO,
    nIter, nCorr, nFcn, nFcnJ, nJac, nRejR1, newt, iConv, qBDamp, stats)

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
    qMixIO      = mPr != 0 && mPrSol != 0
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
    if qBDamp
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
        stats[STATS_NITER] = nIter
        nCorr  = 0
        nRejR1 = 0
        nFcn   = 0
        nFcnJ  = 0
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
        if mPr >= 2 && !qMixIO
            write(printIO,
            "\n",
            "*******************************************************************",
            "\n",
            "\tIt\tNormf\t\tNormx\t\tDamp.Fct.\tNew\n")
        end
        # ----------------------------------------------------------------------
        # 1.7 Startup step
        # ----------------------------------------------------------------------
        # 1.7.1 Computation of residual vector
        try
            fcn(x,f)
        catch
            retCode = 82
            qIter   = false
            #throw(EvaluationError(fcn,err))
        end
        nFcn += 1
        if length(f) != n
            retCode = 22
            qIter   = false
            #=throw(DimensionMismatch("Dimension of the function output does not
            match the input dimension. Please check the function \"$fcn\" again."))=#
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

                        if mPr >= 5
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
        # TODO: Multiple dispatch for the Jacobian evaluation
        # ----------------------------------------------------------------------
        # 2.3.1 Jacobian generation by routine jac or difference approximation
        # if qGenJ == true
        # - or -
        # Rank-1 update of Jacobian if qGenJ == false
        if qGenJ && (!qSimpl || nIter == 0)
            newt = 0
            if jacGen == 1
                jac = getOption(opt,OPT_JACFCN,0)
                try
                    jac(x,a)
                catch
                    retCode = 82
                    iFail   = -1
                    #throw(EvaluationError(jac,err))
                end
            else
                if mStor == 0
                    if jacGen == 3
                        (a,eta,nFcnJ,iFail) = n1jcf(fcn,n,n,x,f,xw,eta,etaMin,etaMax,etaDif,conv,nFcnJ)
                    end
                    if jacGen == 2
                        (a,nFcnJ,iFail) = n1jac(fcn,n,n,x,f,xw,aJdel,aJmin,nFcnJ)
                    end
                elseif mStor == 1
                    if jacGen == 3
                        (a,eta,nFcnJ,iFail) = n1jcfb(fcn,n,m1,ml,x,f,xw,eta,etaMin,etaMax,etaDif,conv,nFcnJ)
                    end
                    if jacGen == 2
                        (a,nFcnJ,iFail) = n1jacb(fcn,n,m1,ml,x,f,xw,aJdel,aJmin,nFcnJ)
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

        elseif !qSimpl
            newt += 1
        end

        if newt == 0 && (qLU || nIter == 0)
            # ------------------------------------------------------------------
            # 2.3.2.1 Save scaling values
            xwa = xw[1:n]
            if issparse(a)
                nza = nnz(a)
                (row,col) = findn(a)
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
                fw = ones(n)
            end
        end
        # ----------------------------------------------------------------------
        # 2.4.3 Save and scale values of F(n)
        fa[:] = f
        t1 = f.*fw
        # ----------------------------------------------------------------------
        # 3 Central part of iteration step
        # ----------------------------------------------------------------------
        # 3.1 Solution of the linear system
        # ----------------------------------------------------------------------
        # 3.1.1 Decomposition of (n,n) matrix A
        if newt == 0 && (qLU || nIter == 0)
            (l,u,p,iFail) = n1fact(n,m1,ml,mu,a,opt)
            if iFail != 0
                if iFail == 1
                    retCode = 1
                else
                    retCode = 80
                end
                break
            end
        end
        qLInit = true
        # ----------------------------------------------------------------------
        # 3.1.2 Solution of (n,n) system
        if newt == 0
            (t1,iFail) = n1solv(n,m1,ml,mu,l,u,p,t1,opt)
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
            (dx,conv,sumX,dLevF) = n1lvls(n,dx,t1,xw,f,mPr,newt == 0)
        else
            (dx,conv,sumX,dLevF) = n1lvls(n,dx,t1,xwa,f,mPr,newt == 0)
        end
        stats[STATS_SUMX]   = sumX
        stats[STATS_DLEVF]  = dLevF
        xa[:]    = x
        sumXa = sumX
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

            if mPr >= 5
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
                fcbh = fcA*fcBand
                if fc > fcbh
                    fc = fcbh
                    if mPr >= 4
                        write(printIO, "*** Increase rest. act. (a priori)\n")
                    end
                end
                fcbh = fcA/fcBand
                if fc < fcbh
                    fc = fcbh
                    if mPr >= 4
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
                if mPr >= 4
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
                        break
                    else
                        qmStop = true
                    end
                end
            end
        end
        fcMon = fc
        # ----------------------------------------------------------------------
        # 3.4 Save natural level for later computations of corrector
        # and print iterate
        fcNumK = sumX
        if mPr >= 2
            n1prv1(dLevF,dLevXa,fcKeep,nIter,newt,mPr,printIO,qMixIO)
        end
        nRed    = 0
        qNext   = false
        qRep    = false
        qRed    = true
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
            catch
                iFail = -1
                retCode = 82
                # throw(EvaluationError(fcn,err))
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
                if mPr >= 2
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
                        if mPr >= 4
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
                        n1lvls(n,dxQ,t1,xw,f,mPr,newt==0)
                    else
                        (dxQ,conv,sumX,dLevFn) =
                        n1lvls(n,dxQ,t1,xwa,f,mPr,newt==0)
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

                    if mPr >= 5
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
                    if mPr >= 3
                        n1prv2(dLevFn,sqrt(sumX/n),fc,niter,mPr,
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
                            if mPr >= 4
                                write(printIO,
                                " *** Decrese rest. sct. (a posteriori) ***\n")
                            end
                        end
                    end
                    fcMon = fc

                    if mPr >= 5
                        write(printIO,
                        " +++ corrector setting 1 +++\n",
                        @sprintf("fc    = %18.10e\n",fc),
                        " +++++++++++++++++++++++++++\n")
                    end

                    qRep = true
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
                        if mPr >= 3
                            n1prv2(dLevFn,sqrt(sumX/n),fc,nIter,
                            mPr,printIO,qMixIO,"+")
                        end
                        fc = fcCor

                        if mPr >= 5
                            write(printIO,
                            " +++ corrector setting 2 +++\n",
                            @sprintf("fc    = %18.10e\n",fc),
                            " +++++++++++++++++++++++++++\n")
                        end

                        qRep = true
                    else
                        qNext = true
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
            if mPr >= 2
                write(printIO,
                @sprintf("        %2i not accepted damping factor %7.5f",nIter,fc),
                @sprintf("    %2i\n",newt))
            end
            fc  = fcKeep
            fcA = fck2
            if nIter == 0
                fc = fcMin
            end
            qGenJ = true
        else
            # ------------------------------------------------------------------
            # 4 Preparations to start the following iteration step
            # ------------------------------------------------------------------
            # Print values
            if mPr >= 3 && !qOrdi
                n1prv2(dLevFn,sqrt(sumX/n),fc,nIter+1,mPr,
                printIO,qMixIO,"*")
            end
            # Print the natural level of the current iterate and
            # return it in one-step mode
            sumXs = sumX
            sumX = sumXa
            if mPrSol >= 2 && nIter != 0
                n1sout(n,xa,2,opt,wk,mPrSol,printIO)
            elseif mPrSol >= 1 && nIter == 0
                n1sout(n,xa,1,opt,wk,mPrSol,printIO)
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
                qSucc = true
                return
            end
        end
    end
    # End repeat
    # End of main iteration loop
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------

    # TODO: Check which exit messages are not requried due to error()

    # --------------------------------------------------------------------------
    # 9 Exits
    # --------------------------------------------------------------------------
    # 9.1 Solution exit
    aprec = -1.0

    if retCode == 0 || retCode == 4
        if nonLin != 1
            if !qOrdi
                if retCode == 0
                    aprec = sqrt(sumX/n)
                    x += dxQ
                    push!(xIter,x)
                else
                    aprec = sqrt(sumXa/n)
                    if alphaA > 0.0 && iOrMon == 3
                        x += dx
                        push!(xIter,x)
                    end
                end
                # Print final monitor output
                if mPr >= 2
                    if retCode == 0
                        n1prv2(dLevFn,sqrt(sumX/n),fc,nIter+1,
                        mPr,printIO,qMixIO,"*")
                    elseif iOrMon == 3
                        n1prv1(dLevFn,sqrt(sumXa/n),fc,nIter,newt,
                        mPr,printIO,qMixIO)
                    end
                end
                if iOrMon >= 2
                    if iConv <= 1 && alphaE != 0.0 && alphaK != 0.0
                        retCode = 5
                    end
                end
            else
                # if qOrdi is true
                aprec = sqrt(sumXa/n)
            end
            if mPr >= 1
                if qOrdi || retCode == 4
                    nOut = nIter
                else
                    nOut = nIter + 1
                end
                write(printIO,"\n\n\n ",
                @sprintf("Solution of nonlinear system of equations obtained within "),
                @sprintf("%3i iteration steps\n\n",nOut),
                @sprintf("Achieved relative accuracy %10.3e\n",aprec))
            end
        else
            if mPr >= 1
                write(printIO,"\n\n\n ",
                @sprintf("Solution of nonlinear system of equations obtained by NLEQ1\n"),
                @sprintf("No estimate available for the achieved relative accuracy\n"))
            end
        end
    end
    # --------------------------------------------------------------------------
    # 9.2 Fail exit messages
    # --------------------------------------------------------------------------
    # 9.2.1 Termination, since Jacobian matrix became singular
    if retCode == 1 && mPr >= 1
        write(printIO,"\nIteration terminated due to singular Jacobian matrix\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.2 Termination after more than nItmax iterations
    if retCode == 2 && mPr >= 1
        write(printIO,"\n",
        @sprintf("Iteration terminates after nItmax %3i iteration steps\n",nItmax))
    end
    # --------------------------------------------------------------------------
    # 9.2.3 Damping factor fc became too small
    if retCode == 3 && mPr >= 1
        write(printIO,"\n",
        @sprintf("Damping factor has become too small: lambda = %10.3e\n",fc))
    end
    # --------------------------------------------------------------------------
    # 9.2.4.1 Superlinear convergence slowed down
    if retCode == 4 && mPr >= 1
        if iConv == 2
            ctyp = "superlinear"
        end
        if iConv == 3
            ctyp = "quadratic"
        end
        # TODO: write a wrapper for error(), warn(), info() such that it prints
        # either on STDOUT or in file
        if qmStop
            write(printIO,"\nWARNING: Monotonicity test failed after ",ctyp,
            " convergence was already checked\nrTol requirement may be too",
            " stringent\n")
        else
            write(printIO,"\nWARNING: ",ctyp, " convergence slowed down\n",
            "rTol requirement may be too stringent\n")
        end
    end
    # --------------------------------------------------------------------------
    # 9.2.4.2 Convergence criterion satisfied before superlinear
    #         convergence has been established
    if retCode == 5 && dLevFn == 0.0
        retCode = 0
    end
    if retCode == 5 && mPr >= 1
        write(printIO,"\n",
        "WARNING: No quadratic or superlinear convergence established yet\n",
        "         your solution may perhaps be less accurate\n",
        "         as indicated by the standard error estimate\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.5
    if retCode == 22 && mPr >= 1
        write(printIO,"ERROR: dimensions of startvector and problem function ",
        "output differ:\n",
        @sprintf("      length(x0) = %5i,     length(fcn(x0)) = %5i\n",n,length(f)))
    end
    # --------------------------------------------------------------------------
    # 9.2.6 Error exit due to linear solver routine n1fact
    if retCode == 80 && mPr >= 1
        write(printIO,@sprintf("ERROR: %5i",iFail)," signalled by linear solver n1fact\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.7 Error exit due to linear solver routine n1solv
    if retCode == 81 && mPr >= 1
        write(printIO,@sprintf("ERROR: %5i",iFail)," signalled by linear solver n1solv\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.8 Error exit due to fail of user function fcn
    if retCode == 82 && mPr >= 1
        write(printIO,@sprintf("ERROR: %5i",iFail)," signalled by user function fcn\n")
    end
    # --------------------------------------------------------------------------
    # 9.2.9 Error exit due to fail of user function Jacobian
    if retCode == 83 && mPr >= 1
        write(printIO,@sprintf("ERROR: %5i",iFail)," signalled by user function jac\n")
    end
    if (retCode == 82 || retCode == 83) && nIter <= 1 && mPr >= 1
        write(printIO,"Try to find a better initial guess for the solution\n")
    end
    # --------------------------------------------------------------------------
    # 9.3 Common exit
    if mPr >= 3 && retCode != 0 && retCode != 4 && nonLin != 1
        write(printIO,"\n",@sprintf("Achieved relative accuracy %10.3e\n",conva))
        aprec = conva
    end
    rTol = aprec
    sumX = sumXa
    if mPrSol >= 2 && nIter != 0
        mode = 2
        if qOrdi
            mode = 3
        end
        n1sout(n,xa,mode,opt,wk,mPrSol,printIO)
    elseif mPrSol >= 1 && nIter == 0
        n1sout(n,xa,1,opt,wk,mPrSol,printIO)
    end
    if !qOrdi
        if retCode != 4
            nIter += 1
        end
        dLevF = dLevFn
        # Print solution or final iteration vector
        if mPrSol >= 1
            if retCode == 0
                modefi = 3
            else
                modefi = 4
            end
            n1sout(n,x,modefi,opt,wk,mPrSol,printIO)
        end
    end
    # End of exits

    # --------------------------------------------------------------------------
    # 10 Prepare all the variables for returning
    xScal = xw

    opt.options[OPT_QSUCC]  = Int(qSucc)

    stats[STATS_NITER]  = nIter
    stats[STATS_NCORR]  = nCorr
    stats[STATS_NFCN]   = nFcn
    stats[STATS_NFCNJ]  = nFcnJ
    stats[STATS_NJAC]   = nJac
    stats[STATS_NREJR1] = nRejR1
    stats[STATS_NEW]    = newt
    stats[STATS_ICONV]  = iConv

    return (x, xScal, retCode, stats)
    # End of function n1int
end
