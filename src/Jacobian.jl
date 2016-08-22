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
    a = zeros(n,n)
    fu = zero(x)
    for k = 1:n
        w = x[k]
        su = sign(x[k])
        if su == 0
            su = 1
        end
        u = max(max(abs(x[k]),ajmin),yscal[k])*ajdel*su;
        x[k] = w + u
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

# Factors a double precision band matrix by elimination
function dgbfa(abd,lda,n,ml,mu)
    m = ml + mu + 1
    info = 0
    ipvt = zeros(n)
    # zero intial fill-in columns
    j0 = mu + 2
    j1 = min(n,m) - 1
    if j1 >= j0
        for jz = j0:j1
            i0 = m + 1 - jz
            abd[i0:ml,jz] = zeros(ml-i0+1,1)
        end
    end
    jz = j1
    ju = 0
    # Gaussian elimination with partial pivoting
    nm1 = n - 1
    if nm1 >= 1
        for k = 1:nm1
            kp1 = k + 1
            # zero next fill-in column
            jz += 1
            if jz <= n && ml >= 1
                abd[1:ml,jz] = zeros(ml)
            end
            # find l = pivot index
            lm = min(ml,n-k)
            l = idamax(lm+1,abd[m:m+lm,k],1) + m - 1
            ipvt[k] = l + k - m
            # zero pivot implies this column already triangularized
            if abd[l,k] != 0.0
                # interchange if necessary
                if l != m
                    t = abd[l,k]
                    abd[l,k] = abd[m,k]
                    abd[m,k] = t
                end
                # compute multipliers
                t = -1.0/abd[m,k]
                abd[m+1:m+lm,k] = dscal(lm,t,abd[m+1:m+lm,k],1)
                # row elimintaion with column indexing
                ju = min(max(ju,mu+ipvt[k]),n)
                mm = m
                if ju >= kp1
                    for j = kp1:ju
                        l = l - 1
                        mm = mm - 1
                        t = abd[l,j]
                        if l!= mm
                            abd[l,j] = abd[mm,j]
                            abd[mm,j] = t
                        end
                        abd[mm+1:mm+lm,j] = daxpy(lm,t,abd[m+1:m+lm,k],1,abd[mm+1:mm+lm,j],1)
                    end
                end
            else
                info = k
            end
        end
    end
    ipvt[n] = n
    if abs[m,n] == 0.0
        info = n
    end
    return (abd,ipvt,info)
end

# Solves the double precision band system
# a*x = b or transpose(a)*x = b depending on flag
function dgbsl(abd,lda,n,ml,mu,ipvt,b,flag)
    m = mu + ml + 1
    nm1 = n - 1
    if flag == 0
        # solve a*x = b
        # first solve l*y = b
        if ml != 0 && nm1 >= 1
            for k = 1:nm1
                lm = min(ml,n-k)
                l = ipvt[k]
                t = b[l]
                if l != k
                    b[l] = b[k]
                    b[k] = t
                end
                b[k+1:k+lm] = daxpy(lm,t,abd[m+1:m+lm,k],1,b[k+1:k+lm],1)
            end
        end
        # now solve u*x = y
        for kb = 1:n
            k = n + 1 - kb
            b[k] = b[k]/abd[m,k]
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b[k]
            b[lb:lb+lm-1] = daxpy(lm,t,abd[la:la+lm-1,k],1,b[lb:lb+lm-1],1)
        end
    else
        # flag = nonzero => solve transpose(a)*x = b
        # first solve transpose(u)*y = b
        for k = 1:n
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = Base.LinAlg.BLAS.dot(lm,abd[la:la+lm-1,k],1,b[lb:lb+lm-1],1)
            b[k] = (b[k]-t)/abd[m,k]
        end
        # now solve transpose(l)*x = y
        if ml != 0 && nm1 >= 1
            for kb = 1:n
                k = n - kb
                lm = min(ml,n-k)
                b[k] += Base.LinAlg.BLAS.dot(lm,abd[m+1:m+lm,k],1,b[k+1:k+lm],1)
                l = ipvt[k]
                if l != k
                    t = b[l]
                    b[l] = b[k]
                    b[k] = t
                end
            end
        end
    end
    x = b
    return x
end

# TODO: replace with BLAS functions from Julia
# a*x+y (x and y vectors but a is scalar const)
function daxpy(n,da,dx,incx,dy,incy)
    if n > 0
        if incx != 1
            dyout = zeros(dy)
        end
        if incx > 0 && incy > 0
            dyout[1:incy:(n-1)*incy+1] = dy[1:incy:(n-1)*incy+1]+da*dx[1:incx:(n-1)*incx+1]
        else
            error("daxpy: incx and incy must be positive")
        end
    else
        dyout = dy
    end
    return dyout
end

# given vector x and scalar a - output a*x
function dscal(n,da,dx,incx)
    if n > 0
        if incx != 1
            dxout = zeros(dx)
        end
        if incx > 0
            dxout[1:incx:n*incx] = da*dx[1:incx:n*incx]
        else
            error("dscal: incx must be positive")
        end
    end
    return dxout
end

# find the index of element having max abs value
function idamax(n,dx,incx)
    if n < 1
        return 0
    end
    if n == 1
        return 1
    end

    # general increment
    ix = 1
    idamx = 1
    dmax = abs(dx[1])
    ix = ix + incx
    for i = 2:n
        if abs(dx[ix]) > dmax
            idamx = i
            dmax = abs(dx[ix])
        end
        ix = ix + incx
    end
    return idamx
end
