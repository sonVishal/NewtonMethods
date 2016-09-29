function n2prjn(n, iRank, u, d, qe, p, v)
    # Begin
    v[:] = u[p]
    iRk1 = iRank + 1
    del = 0.0
    for i = iRk1:n
        sh = sum(qe[1:i-1,i].*v[1:i-1])
        s = (v[i] - sh)/d[i]
        del += s*s
        v[i] = s
    end
    return del
end

function deccon(a, nRow, nCol, mCon, m, n, iRankC, iRank, cond, d, pivot,
    kRed, ah)
    # Begin
    # --------------------------------------------------------------------------
    # Machine constants
    epMach  = 1e-17
    small   = 1e-150
    # --------------------------------------------------------------------------
    # 1 Initialization
    level = 1
    k     = 1
    dd    = 1
    jj    = 1
    i1    = 1
    hMax  = 0.0
    iFail = 0
    v     = zeros(n)

    if iRank > n
        iRank = n
    end
    if iRank > m
        iRank = m
    end
    # --------------------------------------------------------------------------
    # 1.1 Special case m = 1 and n = 1
    if m == 1 && n == 1
        pivot[1] = 1
        d[1] = a[1,1]
        cond = 1.0
        return (iRankC, iRank, cond, cond, iFail)
    end
    if kRed >= 0
        # ----------------------------------------------------------------------
        # 1.2 Initialize pivot-array
        pivot[1:n] = collect(1:n)
        # ----------------------------------------------------------------------
        # 2. Constrained Householder triangularization
        jd = 1
        iRanC1 = iRankC + 1
        mh = mCon
        iRankH = iRankC
        iData = 0
        if mh == 0
            iRankH = iRank
            mh = m
            iData = 1
        end
        iRk1 = iRank
        for k = 1:iRk1
            level = 1
            if k != n
                k1 = k+1
                qLoop = true
                while qLoop
                    if jd != 0
                        for j = k:n
                            d[j] = sum(a[k:mh,j].^2)
                        end
                    end
                    # ----------------------------------------------------------
                    # 2.1 Column pivoting
                    s1 = d[k]
                    jj = k
                    for l1 = k:n
                        if d[l1] > s1
                            s1 = d[l1]
                            jj = l1
                        end
                    end
                    h = d[jj]
                    if jd == 1
                        hMax = h/max(1.0e1,cond*0.05)
                    end
                    jd = 0
                    if h < hMax
                        jd = 1
                    end
                    qLoop = !(h >= hMax)
                end
                if jj != k
                    # ----------------------------------------------------------
                    # Column interchange
                    i = pivot[k]
                    pivot[k] = pivot[jj]
                    pivot[jj] = i
                    d[jj] = d[k]
                    for l1 = 1:m
                        s1 = a[l1,jj]
                        a[l1,jj] = a[l1,k]
                        a[l1,k] = s1
                    end
                end
            end
            # end for k != n case
            h = sum(a[k:mh,k].^2)
            t = sqrt(h)
            # ------------------------------------------------------------------
            # 2.3.0 A-priori test on pseudo-rank
            if k == 1 || k == iRanC1
                dd = t/cond
            end
            if t <= dd || k > iRankH
                # --------------------------------------------------------------
                # Rank reduction
                iRankH = k-1
                if mh != mCon || iData == 1
                    iRank = iRankH
                    if iRankC == iRank
                        level = 4
                    else
                        level = 3
                    end
                else
                    iRankC = iRankH
                    if iRankC != mCon
                        mh = m
                        iRankH = iRank
                        jd = 1
                        iData = 1
                        continue
                    else
                        error("Internal Error of deccon")
                    end
                end
            end
            if level == 1
                # --------------------------------------------------------------
                # 2.4 Householder step
                s = a[k,k]
                t = -sign(s)*abs(t)
                d[k] = t
                a[k,k] = s-t
                if k != n
                    t = 1.0/(h-s*t)
                    for j = k1:n
                        s = sum(a[k:mh,k].*a[k:mh,j])
                        s *= t
                        s1 = -s
                        if s != 0.0
                            # Update the sub columns
                            a[k:m,j] = a[k:m,j]+a[k:m,k]*s1
                        end
                        # update sub column norms
                        d[j] = d[j]-a[k,j]^2
                    end
                    if k == iRankC
                        mh = m
                        jd = 1
                        iRankH = iRank
                    end
                    if k == iRk1
                        level = 3
                    end
                else
                    level = 4
                end
                # end for Householder step
            end
            if level > 1
                break
            end
        end
    else
        k = -1
        level = 3
    end
    # --------------------------------------------------------------------------
    # 3 Rank-deficient psuedo-inverse
    if level == 3
        iRk1 = iRank+1
        for j = iRk1:n
            for ii = 1:iRank
                i = iRk1 - ii
                s = a[i,j]
                if ii != 1
                    sh = sum(a[i,i1:iRank]'.*v[i1:iRank])
                    s = s - sh
                end
                i1 = i
                v[i] = s/d[i]
                ah[i,j] = v[i]
            end
            for i = iRk1:j
                s = sum(ah[1:i-1,i].*v[1:i-1])
                if i != j
                    v[i] = -s/d[i]
                    ah[i,j] = -v[i]
                end
            end
            if s > -1.0
                d[j] = sqrt(s+1.0)
            else
                iFail = -2
                return (iRankC, iRank, cond, v[1], iFail)
            end
        end
    end
    # --------------------------------------------------------------------------
    # 9 Exit
    if iRankC != 0
        sh = d[iRankC]
        if sh != 0.0
            sh = abs(d[1]/sh)
        end
    else
        sh = 0.0
    end
    v[1] = sh
    if k == iRank
        t = d[iRank]
    end
    if iRankC+1 <= iRank && t != 0.0
        s = abs(d[iRankC+1]/t)
    else
        s = 0.0
    end
    cond = s
    iFail = 0
    return (iRankC, iRank, cond, v[1], iFail)
end

function solcon(a, nRow, nCol, mCon, m, n, x, b, iRankC, iRank, d, pivot,
    kRed, ah)
    # Begin
    v = zeros(n)
    s = 0.0
    # --------------------------------------------------------------------------
    # 1 Solution for pseudo-rank zero
    if iRank == 0
        x[1:n] = zeros(n)
        return (iRankC, iRank)
    end
    if iRank <= iRankC && iRank != n
        iRanC1 = iRankC + 1
    end
    if kRed >= 0 && (m != 1 || n != 1)
        # ----------------------------------------------------------------------
        # 2 Constrained householder transformations of right-hand side
        mh = mCon
        if iRankC == 0
            mh = m
        end
        for j = 1:iRank
            s = sum(a[j:mh,j].*b[j:mh])
            s = s/(d[j]*a[j,j])
            b[j:m] += a[j:m,j]*s
            if j == iRankC
                mh = m
            end
        end
    end
    # --------------------------------------------------------------------------
    # 3 Solution of upper triangular system
    iRk1 = iRank + 1
    for ii = 1:iRank
        i = iRk1 - ii
        i1 = i + 1
        s = b[i]
        if ii != 1
            sh = sum(a[i,i1:iRank]'.*v[i1:iRank])
            s = s-sh
        end
        v[i] = s/d[i]
    end
    if iRank != n && iRank != iRankC
        # ----------------------------------------------------------------------
        # 3.2 Computation of the best constrained least squares-solution
        for j = iRk1:n
            s = sum(ah[1:j-1].*v[1:j-1])
            v[j] = -s/d[j]
        end
        j1 = 1
        for jj = 1:n
            j = n-jj+1
            if jj != 1
                s = sum(ah[j,j1:n]'.*v[j1:n])
            end
            if jj != 1 && j <= iRank
                v[j] -= s
            else
                j1 = j
                v[j] = -(s+v[j])/d[j]
            end
        end
    end
    # --------------------------------------------------------------------------
    # 4 Back-permutation of solution components
    x[:] = v[pivot]
    return (iRankC, iRank)
end

function n1fact(n,lda,ml,mu,a,opt)
    # Begin
    mStor = opt.options[OPT_MSTOR]
    if mStor == 0
        try
            (l,u,p) = lu(a)
            iFail = 0
        catch
            iFail = 1
        end
    elseif mStor == 1
        # Band mode: l holds the complete lu-factorization of p*a
        l = a[:,:]
        # Use LINPACK function to compute the LU in place
        # p stores the pivot vectors and not the permuted identity
        (l,p,iFail) = dgbfa(a,lda,n,ml,mu)
        u = []
    end
    if iFail != 0
        iFail = 1
    end
    return (l,u,p,iFail)
end

function n2fact(n,lda,ldaInv,ml,mu,a,aInv,cond,iRank,opt,p,d,iRepeat)
    # Begin
    mPrWarn = opt.options[OPT_PRINTWARNING]
    printIO = opt.options[OPT_PRINTIOWARN]
    mCon = 0
    iRepeat = -iRepeat
    # TODO: This next line is not clear tmp = IWK(2) which is the variable iRankC from deccon
    tmp = n
    if iRepeat == 0
        tmp = mCon
    end
    (iRankC, iRank, cond, subCond, iFail) = deccon(a, lda, n, mCon, n, n,
                                            tmp, iRank, cond, d, p,
                                            iRepeat, aInv)
    if iFail == -2 && mPrWarn == 1
        write(printIO,"\n","Deccon failed to compute Rank-deficient QR-Decomposition\n")
    end
    if iRank != 0
        cond = abs(d[1]/d[iRank])
        setOption!(wkNLEQ2, WK_SENS1, abs(d[1]))
    else
        cond = 1.0
        setOption!(wkNLEQ2, WK_SENS1, 0.0)
    end
    # TODO: Reassign output variables either here or where the function is called
    return (cond, iFail)
end

function n1solv(n,lda,ml,mu,l,u,p,b,opt)
    # Begin
    mStor = opt.options[OPT_MSTOR]
    if mStor == 0
        x = b[p]
        x = l\x
        x = u\x
    elseif mStor == 1
        x = dgbsl(l,lda,n,ml,mu,p,b,0)
    end
    iFail = 0
    return (x,iFail)
end

function n2solv(n,lda,ldaInv,ml,mu,a,aInv,b,z,iRank,opt,iRepeat,d,pivot)
    # Begin
    mCon = 0
    iRepeat = -iRepeat
    # TODO: Don't know this yet
    tmp = 0
    (iRankC,iRank) = solcon(a,lda,n,mCon,n,n,z,b,tmp,iRank,d,pivot,
        iRepeat,aInv)
    iFail = 0
    # TODO: Reassign output variables either here or where the function is called
    return iFail
end

function n1Prv1(dlevf,dlevx,fc,niter,newt,mPr,printIO,qMixIO)
    # Begin
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
        if mPr >= 3
            write(printIO,"        It       Normf           Normx                     New\n")
        end
        if mPr == 2
            write(printIO,"        It       Normf           Normx         Damp.Fct.   New\n")
        end
    end
    if mPr >= 3 || niter == 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e                 %2i\n",niter,dlevf,dlevx,newt))
    end
    if mPr == 2 && niter != 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e      %7.5f    %2i\n",niter,dlevf,dlevx,fc,newt))
    end
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
    end
    return nothing
end

function n2Prv1(dlevf,dlevx,fc,niter,newt,iRank,mPr,printIO,qMixIO,cond1)
    # Begin
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
        if mPr >= 3
            write(printIO,"        It       Normf           Normx                     New      Rank        Cond\n")
        end
        if mPr == 2
            write(printIO,"        It       Normf           Normx         Damp.Fct.   New      Rank        Cond\n")
        end
    end
    if mPr >= 3 || niter == 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e",niter,dlevf,dlevx),
        @sprintf("                 %2i      %4i          %10.3e\n",newt,iRank,cond1))
    end
    if mPr == 2 && niter != 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e",niter,dlevf,dlevx),
        @sprintf("      %7.5f    %2i      %4i          %10.3e\n",fc,newt,iRank,cond1))
    end
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
    end
    return nothing
end
