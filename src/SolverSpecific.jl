"""
# Summary:
function dgbfa{T}(abd::Array{T,2}, lda::Int64, n::Int64, ml::Int64, mu::Int64)

This function is similar to the function gbtrf! provided in Base.LinAlg.LAPACK.
Julia documentation mentions the following:

*Note that the LAPACK API provided by Julia can and will change in the future.*
*Since this API is not user-facing, there is no commitment to support/deprecate*
*this specific set of functions in future releases.*

Due to this reason dgbfa function has been written as per the one provided in the
original version of these solvers at: [NewtonLib](http://elib.zib.de/pub/elib/codelib/NewtonLib/index.html)
based on LINPACK written for Matlab by Cleve Moler, University of New Mexico, Argonne National Lab.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable    | Description                                                  |
|:------------|:-------------------------------------------------------------|
| abd[lda,n]* | Contains the matrix in band storage. See the comments below. |
| lda         | Leading dimension of the matrix. lda >= 2*ml + mu + 1        |
| n           | Order of the original matrix                                 |
| ml          | Number of diagonals below the main diagonal. 0 <= ml < n     |
| mu          | Number of diagonals above the main diagonal. 0 <= mu < n     |

(* marks inout parameters)

## Output parameters
--------------------
| Variable    | Description                                                  |
|:------------|:-------------------------------------------------------------|
| abd[lda,n]* | An upper triangular matrix in band storage and the multipliers which were used to obtain it. Factorization can be written as a = l*u where l is a product of permutation and unit lower triangular matrices and u is upper triangular. |
| ipvt        | Vector containing the pivot indices                          |
| info = 0    | Normal value                                                 |
| info = k    | if u(k,k) = 0.0. This only indicates that dgbsl will divide by zero if called. |

## Band storage
If a is a band matrix, the following program segment will set up input matrix abd
```
ml = (bandwidth below the diagonal)
mu = (bandwidth above the diagonal)
m = ml + mu + 1
for j = 1:n
    i1 = max(1,j-mu)
    i2 = min(n,j+ml)
    for i = i1:i2
        k = i - j + m
        abd[k,j] = a[i,j]
    end
end
```
"""
function dgbfa{T}(abd::Array{T,2}, lda::Int64, n::Int64, ml::Int64, mu::Int64)
    # Begin
    # Copy the matrix so that it remains unchanged at exit
    ipvt = zeros(Int64,n)
    a = abd[:,:]
    m = ml + mu + 1
    info = 0
    # zero intial fill-in columns
    j0 = mu + 2
    j1 = min(n,m) - 1
    if j1 >= j0
        for jz = j0:j1
            i0 = m + 1 - jz
            a[i0:ml,jz:jz] = zeros(T,ml-i0+1,1)
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
                a[1:ml,jz:jz] = zeros(T,ml)
            end
            # find l = pivot index
            lm = min(ml,n-k)
            l = maxabs(a[m:m+lm,k:k]) + m - 1
            ipvt[k] = l + k - m
            # zero pivot implies this column already triangularized
            if a[l,k] != 0.0
                # interchange if necessary
                if l != m
                    t = a[l,k]
                    a[l,k] = a[m,k]
                    a[m,k] = t
                end
                # compute multipliers
                t = -1.0/a[m,k]
                a[m+1:m+lm,k:k] *= t
                # row elimintaion with column indexing
                ju = min(max(ju,mu+ipvt[k]),n)
                mm = m
                if ju >= kp1
                    for j = kp1:ju
                        l = l - 1
                        mm = mm - 1
                        t = a[l,j]
                        if l!= mm
                            a[l,j] = a[mm,j]
                            a[mm,j] = t
                        end
                        a[mm+1:mm+lm,j:j] += t*a[m+1:m+lm,k:k]
                    end
                end
            else
                info = k
            end
        end
    end
    ipvt[n] = n
    if a[m,n] == 0.0
        info = n
    end
    return (a,ipvt,info)
end

"""
function dgbsl{T}(abd::Array{T,2}, lda::Int64, n::Int64, ml::Int64, mu::Int64,
    ipvt::Vector{Int64}, b::Vector{T}, flag::Int64)

Solves the double precision band system using the factors computed by dgbfa.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable    | Description                                                  |
|:------------|:-------------------------------------------------------------|
| abd[lda,n]  | Matrix output from dgbfa                                     |
| lda         | Leading dimension of the matrix abd                          |
| n           | Order of the original matrix                                 |
| ml          | Number of diagonals below the main diagonal                  |
| mu          | Number of diagonals above the main diagonal                  |
| ipvt        | Vector containing the pivot indices from dgbfa               |
| b           | Right hand side vector                                       |
| info = 0    | Solve a*x = b                                                |
| info != 0   | Solve trans(a)*x = b                                         |

## Output parameters
--------------------
| Variable    | Description                                                  |
|:------------|:-------------------------------------------------------------|
| x           | Solution vector                                              |
"""
function dgbsl{T}(abd::Array{T,2}, lda::Int64, n::Int64, ml::Int64, mu::Int64,
    ipvt::Vector{Int64}, b::Vector{T}, flag::Int64)
    # Begin
    x = b[:]
    m = mu + ml + 1
    nm1 = n - 1
    if flag == 0
        # solve a*x = b
        # first solve l*y = b
        if ml != 0 && nm1 >= 1
            for k = 1:nm1
                lm = min(ml,n-k)
                l = ipvt[k]
                t = x[l]
                if l != k
                    x[l] = x[k]
                    x[k] = t
                end
                x[k+1:k+lm] += t*abd[m+1:m+lm,k:k]
            end
        end
        # now solve u*x = y
        for kb = 1:n
            k = n + 1 - kb
            x[k] = x[k]/abd[m,k]
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -x[k]
            x[lb:lb+lm-1] += t*abd[la:la+lm-1,k:k]
        end
    else
        # flag = nonzero => solve transpose(a)*x = b
        # first solve transpose(u)*y = b
        for k = 1:n
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = abd[la:la+lm-1,k:k].*b[lb:lb+lm-1]
            x[k] = (x[k]-t)/abd[m,k]
        end
        # now solve transpose(l)*x = y
        if ml != 0 && nm1 >= 1
            for kb = 1:n
                k = n - kb
                lm = min(ml,n-k)
                x[k] += lm,abd[m+1:m+lm,k:k].*b[k+1:k+lm]
                l = ipvt[k]
                if l != k
                    t = x[l]
                    x[l] = x[k]
                    x[k] = t
                end
            end
        end
    end
    return x
end

"""
function n2prjn{T}(n::Int64, iRank::Int64, u::Vector{T}, d::Vector{T},
    qe::Array{T,2}, p::Vector{Int64}, v::Vector{T})

Provides the projection to the appropriate subspace in case
    of rank - reduction. To be used in connection with DECCON/SOLCON.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable | Description                                      |
|:---------|:-------------------------------------------------|
| n        | Number of parameters to be estimated             |
| iRank    | Pseudo rank of decomposed Jacobian matrix        |
| u[n]     | Scaled Newton correction                         |
| d[n]     | Diagonal elements of upper triangular matrix     |
| qe[n]    | Part of pseudoinverse Jacobian matrix            |
| p[n]     | Pivot vector resulting from matrix decomposition |
| v[n]     | Temporary work array                             |

## Output parameters
--------------------
| Variable | Description                 |
|:---------|:----------------------------|
| del      | Defect                      |
"""
function n2prjn{T}(n::Int64, iRank::Int64, u::Vector{T}, d::Vector{T},
    qe::Array{T,2}, p::Vector{Int64}, v::Vector{T})
    # Begin
    v[:] = u[p]
    iRk1 = iRank + 1
    del = 0.0
    for i = iRk1:n
        sh = sum(qe[1:i-1,i:i].*v[1:i-1])
        s = (v[i] - sh)/d[i]
        del += s*s
        v[i] = s
    end
    return del
end

"""
function deccon{T}(a::Array{T,2}, nRow::Int64, nCol::Int64, mCon::Int64,
    m::Int64, n::Int64, iRankC::Int64, iRank::Int64, cond::T,
    d::Vector{T}, pivot::Vector{Int64}, kRed::Int64, ah::Array{T,2})

Constrained QR-decomposition of (m,n)-system  with computation of pseudoinverse
in case of rank-defeciency. First mcon rows belong to equality constraints.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable      | Description                                                   |
|:--------------|:--------------------------------------------------------------|
| a[nRow,nCol]* | Array holding the (m,n) matrix to be decomposed               |
| nRow          | Declared number of rows of array a                            |
| nCol          | Declared number of columns of array a and rows and columns of ah |
| mCon          | Number of equality constraints. mCon <= n. Internally reduced if equality constraints are linearly dependent |
| m             | Current number of rows of array a                             |
| n             | Current number of columns of array a                          |
| iRankC*       | Prescribed maximum pseudo-rank of constrained part of a. iRankC <= mCon |
| iRank*        | Prescribed maximum pseudo-rank of matrix a. iRank <= n        |
| cond*         | Permitted upper bound for the subcondition of least squares part of a |
| kRed >= 0     | Householder triangularization. Build up pseudo-inverse if iRank < n |
| kRed < 0      | Reduction of pseudo-rank of matrix a, skipping Householder triangularization, build-up new pseudo-inverse |

(* marks inout parameters)

## Output parameters
--------------------
| Variable      | Description                                                   |
|:--------------|:--------------------------------------------------------------|
| a[nRow,nCol]* | Upper triangular part consisting of transformed matrix and lower triangular part consisting of Householder transformations |
| iRankC*       | New pseudo-rank of constrained part of a such that abs(d[1]/d[iRankC]) < 1/epMach |
| iRank*        | New pseudo-rank of matrix a such that abs(d[iRankC+1]/d[iRank]) < cond |
| cond*         | The subcondition number belonging to the least squares part of A. In case of rank reduction: subcondition number which led to rank reduction |
| d             | Diagonal elements of the upper triangular matrix              |
| pivot         | Index vector storing permutation of columns due to pivoting   |
| ah[nCol,nCol] | In case of rank-defect used to compute the psuedo-inverse     |
| iFail = 0     | deccon computation was successful                             |
| iFail = -2    | Numerically negative diagonal element encountered during computation of pseudo-inverse due to extremely bad conditioned matrix a. deccon is unable to continue rank-reduction |
"""
function deccon{T}(a::Array{T,2}, nRow::Int64, nCol::Int64, mCon::Int64,
    m::Int64, n::Int64, iRankC::Int64, iRank::Int64, cond::T,
    d::Vector{T}, pivot::Vector{Int64}, kRed::Int64, ah::Array{T,2})
    # Begin
    # --------------------------------------------------------------------------
    # 1 Initialization
    level = 1
    k     = 1
    dd    = 1
    jj    = 1
    i1    = 1
    hMax  = zero(T)
    iFail = 0
    v     = zeros(T,n)
    t     = zero(T)

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
        cond = one(T)
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
                            d[j] = sum(a[k:mh,j:j].^2)
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
            h = sum(a[k:mh,k:k].^2)
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
                        s = sum(a[k:mh,k:k].*a[k:mh,j:j])
                        s *= t
                        s1 = -s
                        if s != 0.0
                            # Update the sub columns
                            a[k:m,j:j] = a[k:m,j:j]+a[k:m,k:k]*s1
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
                    sh = sum(a[i:i,i1:iRank]'.*v[i1:iRank])
                    s = s - sh
                end
                i1 = i
                v[i] = s/d[i]
                ah[i,j] = v[i]
            end
            for i = iRk1:j
                s = sum(ah[1:i-1,i:i].*v[1:i-1])
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
        sh = zero(T)
    end
    v[1] = sh
    if k == iRank
        t = d[iRank]
    end
    if iRankC+1 <= iRank && t != 0.0
        s = abs(d[iRankC+1]/t)
    else
        s = zero(T)
    end
    cond = s
    iFail = 0
    return (iRankC, iRank, cond, v[1], iFail)
end

"""
function solcon{T}(a::Array{T,2}, nRow::Int64, nCol::Int64, mCon::Int64, m::Int64,
    n::Int64, x::Vector{T}, b::Vector{T}, iRankC::Int64, iRank::Int64,
    d::Vector{T}, pivot::Vector{Int64}, kRed::Int64, ah::Array{T,2})

Best constrained linear least squares solution of (m,n)-system. First mcon rows
comprise mcon equality constraints. To be used in connection with subroutine deccon.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable | Description             |
|:---------|:------------------------|
| a[m,n], nRow, nCol, m, n, mCon, iRankC, iRank, d[n], pivot[n], ah[n,n], kRed | Refer to input and output parameters of deccon     |
| b[m]*    | Right hand side of the linear system if kRed >= 0. Right hand side of the upper linear system if kRed < 0 |

(* marks inout parameters)

## Output parameters
--------------------
| Variable | Description                        |
|:---------|:-----------------------------------|
| x[n]     | Best LSQ-solution of linear system |
| b[m]*    | Right-hand side of upper triangular system (transformed right-hand side of linear system) |
"""
function solcon{T}(a::Array{T,2}, nRow::Int64, nCol::Int64, mCon::Int64, m::Int64,
    n::Int64, x::Vector{T}, b::Vector{T}, iRankC::Int64, iRank::Int64,
    d::Vector{T}, pivot::Vector{Int64}, kRed::Int64, ah::Array{T,2})
    # Begin
    v = zeros(T,n)
    s = zero(T)
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
            s = sum(a[j:mh,j:j].*b[j:mh])
            s = s/(d[j]*a[j,j])
            b[j:m] += a[j:m,j:j]*s
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
            sh = sum(a[i:i,i1:iRank]'.*v[i1:iRank])
            s = s-sh
        end
        v[i] = s/d[i]
    end
    if iRank != n && iRank != iRankC
        # ----------------------------------------------------------------------
        # 3.2 Computation of the best constrained least squares-solution
        for j = iRk1:n
            s = sum(ah[1:j-1,j:j].*v[1:j-1])
            v[j] = -s/d[j]
        end
        j1 = 1
        for jj = 1:n
            j = n-jj+1
            s = zero(T)
            if jj != 1
                s = sum(ah[j:j,j1:n]'.*v[j1:n])
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
    x[pivot] = v[:]
    return iRank
end

"""
function nFact{T}(n::Int64, lda::Int64, ml::Int64, mu::Int64, a::Array{T,2},
    mStor::Int64, l::Array{T,2}, u::Array{T,2}, p::Vector{Int64})

Call linear algebra subprogram for factorization of a (n,n)-matrix.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable  | Description                                             |
|:----------|:--------------------------------------------------------|
| n         | Order of the linear system                              |
| lda       | Leading dimension of the matrix a                       |
| ml        | Lower bandwidth of the matrix (only for banded systems) |
| mu        | Upper bandwidth of the matrix (only for banded systems) |
| a[lda,n]  | Matrix to be factorized                                 |
| mStor = 0 | Full storage mode for matrix                            |
| mStor = 1 | Band storage mode for matrix                            |

## Output parameters
--------------------
| Variable | Description                                        |
|:---------|:---------------------------------------------------|
| iFail    | Exit code in case of errors                        |
| = 0      | Matrix decomposition successful                    |
| = 1      | Decomposition failed - matrix numerically singular |
| l[lda,n] | Lower triangular part of decomposed matrix in case of full mode. The complete LU-decomposition in band mode |
| u[lda,n] | Upper triangular part of decomposed matrix in case of full mode. Unused in case of band mode. |
| p        | Vector of pivot indices                            |
"""
function nFact{T}(n::Int64, lda::Int64, ml::Int64, mu::Int64, a::Array{T,2},
    mStor::Int64, l::Array{T,2}, u::Array{T,2}, p::Vector{Int64})
    # Begin
    iFail = 0
    if mStor == 0
        try
            (l[:,:],u[:,:],p[:]) = lu(a)
        catch
            iFail = 1
        end
    elseif mStor == 1
        # Band mode: l holds the complete lu-factorization of p*a
        l[:,:] = a
        # Use LINPACK function to compute the LU in place
        # p stores the pivot vectors and not the permuted identity
        (l[:,:],p[:],iFail) = dgbfa(a,lda,n,ml,mu)
        if iFail != 0
            iFail = 1
        end
    end
    return iFail
end

"""
function nFact{T}(n::Int64, lda::Int64, ldaInv::Int64, ml::Int64, mu::Int64,
    a::Array{T,2}, aInv::Array{T,2}, cond::T, iRank::Int64,
    opt::OptionsNLEQ, p::Vector{Int64}, d::Vector{T}, iRepeat::Int64, iRankC::Int64)

Call linear algebra subprogram for factorization of a (n,n)-matrix
with rank decision and casual computation of the rank deficient
pseudo-inverse matrix.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable  | Description                                             |
|:----------|:--------------------------------------------------------|
| n         | Order of the linear system                              |
| lda       | Leading dimension of the matrix a                       |
| ldaInv    | Leading dimension of the matrix aInv                    |
| ml        | Lower bandwidth of the matrix (only for banded systems) |
| mu        | Upper bandwidth of the matrix (only for banded systems) |
| a[lda,n]  | Matrix to be factorized                                 |
| cond*     | Maximum permitted subcondition for the prescribed rank  |
| opt       | User prescribed options                                 |
| iRankC, iRank, cond | Refer to input and output parameters of deccon |
| iRank*    | Prescribed maximum pseudo-rank of matrix a              |
| cond*     | Permitted upper bound for the subcondition of leastsquares part of a |

## Output parameters
--------------------
| Variable       | Description                                        |
|:---------------|:---------------------------------------------------|
| aInv[ldaInv,n] | If matrix a is rank deficient this array holds the pseudo-inverse of a |
| iFail = 0      | Matrix decomposition successful                    |
| p              | Vector of pivot indices                            |
| d              | Refer to deccon                                    |
"""
function nFact{T}(n::Int64, lda::Int64, ldaInv::Int64, ml::Int64, mu::Int64,
    a::Array{T,2}, aInv::Array{T,2}, cond::T, iRank::Int64,
    opt::OptionsNLEQ, p::Vector{Int64}, d::Vector{T}, iRepeat::Int64, iRankC::Int64)
    # Begin
    mPrWarn = opt.options[OPT_PRINTWARNING]
    printIO = opt.options[OPT_PRINTIOWARN]
    mCon = 0
    iRepeat = -iRepeat
    if iRepeat == 0
        iRankC = mCon
    end
    (iRankC, iRank, cond, subCond, iFail) = deccon(a, lda, n, mCon, n, n,
                                            iRankC, iRank, cond, d, p,
                                            iRepeat, aInv)
    if iFail == -2 && mPrWarn == 1
        write(printIO,"\n","Deccon failed to compute Rank-deficient QR-Decomposition\n")
    end
    if iRank != 0
        cond = abs(d[1]/d[iRank])
        setOption!(wkNLEQ2, WK_SENS1, abs(d[1]))
    else
        cond = one(T)
        setOption!(wkNLEQ2, WK_SENS1, zero(T))
    end
    return (cond, iRankC, iFail)
end

"""
function nSolv{T}(n::Int64, lda::Int64, ml::Int64, mu::Int64, l::Array{T,2},
    u::Array{T,2}, p::Vector{Int64}, b::Vector{T}, mStor::Int64)

Call linear algebra subprogram for solution of the linear system a*z = b

T = Float64 or BigFloat

## Parameters
-------------------
| Variable | Description             |
|:---------|:------------------------|
| n, lda, ml, mu, l, u, p, b, mStor, iFail | Refer nFact |
"""
function nSolv{T}(n::Int64, lda::Int64, ml::Int64, mu::Int64, l::Array{T,2},
    u::Array{T,2}, p::Vector{Int64}, b::Vector{T}, mStor::Int64)
    # Begin
    if mStor == 0
        b[:] = b[p]
        b[:] = l\b
        b[:] = u\b
    elseif mStor == 1
        b[:] = dgbsl(l,lda,n,ml,mu,p,b,0)
    end
    iFail = 0
    return iFail
end

"""
function nSolv{T}(n::Int64, lda::Int64, ldaInv::Int64, ml::Int64, mu::Int64,
    a::Array{T,2}, aInv::Array{T,2}, b::Vector{T},
    z::Vector{T}, iRank::Int64, iRepeat::Int64, d::Vector{T},
    pivot::Vector{Int64}, iRankC::Int64)

Call linear algebra subprogram for solution of the linear system a*z = b.
To be used with the factorization routine of a (n,n)-matrix
with rank decision and casual computation of the rank deficient
pseudo-inverse matrix.

T = Float64 or BigFloat

## Input parameters
-------------------
| Variable | Description                          |
|:---------|:-------------------------------------|
| b[n]*    | Right hand side of the linear system |
| n, lda, ldaInv, ml, mu, a, aInv, iRank, iRepeat, d, pivot, iRankC | Refer nFact |

(* marks inout parameters)

## Output parameters
--------------------
| Variable | Description                                                       |
|:---------|:------------------------------------------------------------------|
| b[n]*    | RHS transformed to the upper triangular part of the linear system |
"""
function nSolv{T}(n::Int64, lda::Int64, ldaInv::Int64, ml::Int64, mu::Int64,
    a::Array{T,2}, aInv::Array{T,2}, b::Vector{T},
    z::Vector{T}, iRank::Int64, iRepeat::Int64, d::Vector{T},
    pivot::Vector{Int64}, iRankC::Int64)
    # Begin
    mCon = 0
    iRepeat = -iRepeat
    iRank = solcon(a,lda,n,mCon,n,n,z,b,iRankC,iRank,d,pivot,
                iRepeat,aInv)
    iFail = 0
    return iFail
end
