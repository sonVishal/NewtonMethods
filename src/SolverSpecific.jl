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

function n1prv1(dlevf,dlevx,fc,niter,newt,mPr,printIO,qMixIO)
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

function n2prv1(dlevf,dlevx,fc,niter,newt,mPr,printIO,qMixIO,cond1,iRank)
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
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e                 %2i      %4i  %10.3e\n",niter,dlevf,dlevx,newt,iRank,cond1))
    end
    if mPr == 2 && niter != 0
        write(printIO,@sprintf("      %4i     %10.3e      %10.3e      %7.5f    %2i      %4i  %10.3e\n",niter,dlevf,dlevx,fc,newt,iRank,cond1))
    end
    if qMixIO
        write(printIO,"  ******************************************************************",
        "\n");
    end
    return nothing
end

function n2prjn()
end

function n2fact()
end

function n2solv()
end


function deccon()
end

function solcon()
end
