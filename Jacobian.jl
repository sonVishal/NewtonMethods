# TODO: write the functions starting at 1357 in f file
# TODO: Multiple dispatch between banded and dense for each => only 2 functions
# TODO: Wrap user given function to give iFail using try catch mechanism

# function numDiffJac(x,J)
#     Evaluation of a dense Jacobian matrix using finite difference
#     approximation adapted for use in nonlinear systems solver NLEQ1

function n1jac(fcn,n,lda,x,fx,yscal,ajdel,ajmin,nFcn)
    # Begin
    iFail = 0
    nFcnNew = nFcn
    for k = 1:n
        w = x[k]
        su = sign(x[k])
        if su == 0
            su = 1
        end
        u = max(max(abs(x[k]),ajmin),yscal[k])*ajdel*su;
        x[k] = w+u
        try
            fcn(x,fu)
            iFail = 0
        catch
            iFail = -1
        end
        nFcnNew += 1
        if iFail != 0
            break
        end
        x[k] = w
        a[1:n,k] = (fu-fx)/u
    end
    return (a,nFcnNew,iFail)
end

# function numDiffJacBanded(x,J)
#     Evaluation of a banded Jacobian matrix using finite difference
#     approximation adapted for use in nonlinear systems solver NLEQ1
function n1jacb(fcn,n,lda,ml,x,fx,yscal,ajdel,ajmin,nFcn)
    # Begin
    nFcnNew = nFcn
    iFail = 0
    mu = lda - 2*ml -1
    ldab = ml + mu + 1
    a = zeros(ldab,n)
    w = zeros(n)
    u = zeros(n)

    for jj = 1:ldab
        for k = jj:ldab:n
            w[k] = x[k]
            su = sign(x[k])
            if su == 0
                su = 1
            end
            u[k] = max(max(abs(x[k]),ajmin),yscal[k])*ajdel*su
            x[k] = w[k] + u[k]
        end
        try
            fcn(x,fu)
            iFail = 0
        catch
            iFail = -1
        end
        if iFail != 0
            break;
        end

        for k = jj:ldab:n
            x[k] = w[k]
            i1 = max(1,k-mu)
            i2 = min(n,k+ml)
            mh = mu + 1 - k
            a[mh+i1:mh+i2,k] = (fu[i1:i2]-fx[i1:i2])/u[k]
        end
    end
    return (a,nFcnNew,iFail)
end

# function numDiffJacFB(x,J)
#     Approximation of dense Jacobian matrix for nonlinear systems solver
#     NLEQ1 with feed-back control of discretization and rounding errors

function n1jcf(fcn,n,lda,fx,yscal,eta,etamin,etamax,etadif,conv,nFcn)
    # Constant
    small2 = 0.1
    # Begin
    iFail = 0
    nFcnNew = nFcn
    a = zeros(n,n)
    for k = 1:n
        is = 0
        qFine = false
        qExit = false
        while !qFine
            w = x[k]
            su = sign(x[k])
            if su == 0
                su = 1
            end
            u = eta[k]*yscal[k]*su
            x[k] = w + u
            try
                fcn(x,fu)
                iFail = 0
            catch
                iFail = -1
            end
            nFcnNew += 1
            if iFail != 0
                qExit = true
                break;
            end
            x[k] = w
            sumd = 0.0
            for i = 1:n
                hg = max(abs(fx[i]),abs(fu[i]))
                fhi = fu[i] - fx[i]
                if hg != 0.0
                    sumd = sumd + (fhi/hg)^2
                end
                a[i,k] = fhi/u
            end
            sumd = sqrt(sumd/n)
            qFine = true
            if sumd != 0.0 && is == 0
                eta[k] = min(etamax,max(etamin,sqrt(etadif/sumd)*eta[k]))
                is = 1
                qFine = conv < small2 || sumd >= etamin
            end
        end
        if qExit
            break;
        end
    end
    return (a,eta,nFcnNew,iFail)
end

# function numDiffJac(x,J)
#     Approximation of dense Jacobian matrix for nonlinear systems solver
#     NLEQ1 with feed-back control of discretization and rounding errors
function n1jcfb(fcn,n,lda,ml,x,fx,yscal,eta,etamin,etamax,etadif,conv,nFcn)
    # Constants
    small2 = 0.1
    # Begin
    nFcnNew = nFcn
    mu = lda - 2*ml - 1
    ldab = ml + mu + 1
    a = zeros(ldab,n)
    w = zeros(n)
    u = zeros(n)

    for jj = 1:ldab
        is = 0
        qFine = false
        qExit = false
        while !qFine
            for k = jj:ldab:n
                w[k] = x[k]
                su = sign(x[k])
                if su == 0
                    su = 1
                end
                u[k] = eta[k]*yscal[k]*su
                x[k] = w[k] + u[k]
            end

            try
                fcn(x,fu)
                iFail = 0
            catch
                iFail = -1
            end
            nFcnNew += 1
            if iFail != 0
                qExit = true
                break;
            end

            for k = jj:ldab:n
                x[k] = w[k]
                sumd = 0.0
                i1 = max(1,k-mu)
                i2 = min(n,k+ml)
                mh = mu + 1 - k
                for i = i1:i2
                    hg = max(abs(fx[i]),abs(fu[i]))
                    fhi = fu[i] - fx[i]
                    if hg != 0.0
                        sumd += (fhi/hg)^2
                    end
                    a[mh+i,k] = fhi/u[k]
                end
                sumd = sqrt(sumd/n)
                qFine = true
                if sumd != 0.0 && is == 0
                    eta[k] = min(etamax,max(etamin,sqrt(etadif/sumd)*eta[k]))
                    is = 1
                    qFine = conv < small2 || sumd >= etamin
                end
            end
        end
        if qExit
            break;
        end
    end
    return (a,eta,nFcnNew,iFail)
end
