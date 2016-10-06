using ForwardDiff

# Evaluation of Jacobian using Automatic differentiation
# TODO: this is not working
"""
function nJacFAD(fcn, x::Vector{Float64}, a::Array{Float64,2}, nFcn::Int64,
    chunkSize = ForwardDiff.pickchunk(x))

Evaluation of Jacobian matrix using forward mode automatic differentiation
for use in nonlinear systems solver.

## Input parameters
-------------------
| Variable  | Description                                                     |
|-----------|-----------------------------------------------------------------|
| fcn       | Function of the form fcn(f, x) to provide right-hand side       |
| x[n]      | Current scaled vector                                           |
| nFcn*     | fcn evaluation count                                            |
| chunkSize | Chunk size for computing the Jacobian. Refer doc of ForwardDiff |

(* marks inout parameters)

## Output parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| a        | Array containing the approximated Jacobian                      |
| nFcn*    | fcn evaluation count adjusted                                   |
| iFail    | Return code non-zero if Jacobian could not be computed          |
"""
function nJacFAD(fcn, x::Vector{Float64}, a::Array{Float64,2}, nFcn::Int64,
    chunkSize = ForwardDiff.pickchunk(x))
    # Begin
    y = zeros(x)
    iFail = 0
    try
        fcn(y,x)
    catch
        iFail = -1
    end
    nFcn += 1
    if iFail == 0
        try
            a = ForwardDiff.jacobian(fcn,y,x,chunk;usecache=true)
        catch
            iFail = -1
        end
    end
    return (nFcn, iFail)
end

"""
function nJac(fcn, n::Int64, lda::Int64, x::Vector{Float64}, fx::Vector{Float64},
    yscal::Vector{Float64}, ajdel::Float64, ajmin::Float64, nFcn::Float64,
    a::Array{Float64,2})

Evaluation of a dense Jacobian matrix using finite difference approximation
adapted for use in nonlinear systems solver.

## Input parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| fcn      | Function of the form fcn(f, x) to provide right-hand side       |
| n        | Number of rows and columns of the Jacobian                      |
| lda      | Leading dimension of the array "a"                              |
| x[n]     | Current scaled vector                                           |
| fx[n]    | Vector containing fcn(x)                                        |
| yscal[n] | Vector containing scaling factors                               |
| ajdel    | Perturbation of component k: abs(y(k))*ajdel                    |
| ajmin    | Minimum perturbation is ajmin*ajdel                             |
| nFcn*    | fcn evaluation count                                            |

(* marks inout parameters)

## Output parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| a[lda,n] | Array containing the approximated Jacobian                      |
| nFcn*    | fcn evaluation count adjusted                                   |
| iFail    | Return code non-zero if Jacobian could not be computed          |
"""
function nJac(fcn, n::Int64, lda::Int64, x::Vector{Float64}, fx::Vector{Float64},
    yscal::Vector{Float64}, ajdel::Float64, ajmin::Float64, nFcn::Float64,
    a::Array{Float64,2})
    # Begin
    iFail = 0
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
            fcn(fu,x)
        catch
            iFail = -1
        end
        nFcn += 1
        if iFail != 0
            break
        end
        x[k] = w
        a[1:n,k] = (fu-fx)/u
    end
    return (nFcn, iFail)
end

"""
function nJacb(fcn, n::Int64, lda::Int64, ml::Int64, x::Vector{Float64},
    fx::Vector{Float64}, yscal::Vector{Float64}, ajdel::Int64, ajmin::Int64,
    nFcn::Int64, a::Array{Float64,2})

Evaluation of a banded Jacobian matrix using finite difference approximation
adapted for use in nonlinear systems solver

## Input parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| fcn      | Function of the form fcn(f, x) to provide right-hand side       |
| n        | Number of rows and columns of the Jacobian                      |
| lda      | Leading dimension of the array "a"                              |
| ml       | Lower bandwidth of the Jacobian matrix                          |
| x[n]     | Current scaled vector                                           |
| fx[n]    | Vector containing fcn(x)                                        |
| yscal[n] | Vector containing scaling factors                               |
| ajdel    | Perturbation of component k: abs(y(k))*ajdel                    |
| ajmin    | Minimum perturbation is ajmin*ajdel                             |
| nFcn*    | fcn evaluation count                                            |

(* marks inout parameters)

## Output parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| a[lda,n] | Array containing the approximated Jacobian                      |
| nFcn*    | fcn evaluation count adjusted                                   |
| iFail    | Return code non-zero if Jacobian could not be computed          |
"""
function nJacb(fcn, n::Int64, lda::Int64, ml::Int64, x::Vector{Float64},
    fx::Vector{Float64}, yscal::Vector{Float64}, ajdel::Int64, ajmin::Int64,
    nFcn::Int64, a::Array{Float64,2})
    # Begin
    iFail = 0
    mu = lda - 2*ml -1
    ldab = ml + mu + 1
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
            fcn(fu,x)
        catch
            iFail = -1
        end
        nFcn += 1
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
    return (nFcn,iFail)
end

"""
function nJcf(fcn, n::Int64, lda::Int64, x::Vector{Float64}, fx::Vector{Float64},
    yscal::Vector{Float64}, eta::Vector{Float64}, etamin::Float64, etamax::Float64,
    etadif::Float64, conv::Float64, nFcn::Int64, a::Array{Float64,2})

Approximation of dense Jacobian matrix for nonlinear systems solver with
feed-back control of discretization and rounding errors

## Input parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| fcn      | Function of the form fcn(f, x) to provide right-hand side       |
| n        | Number of rows and columns of the Jacobian                      |
| lda      | Leading dimension of the array "a"                              |
| x[n]     | Current scaled vector                                           |
| fx[n]    | Vector containing fcn(x)                                        |
| yscal[n] | Vector containing scaling factors                               |
| eta[n]*  | Vector of scaled denominator differences                        |
| etamin   | Minimum allowed scaled denominator                              |
| etamax   | Maximum allowed scaled denominator                              |
| etadif   | = sqrt(1.1*epMach)                                              |
| conv     | Maximum norm of last (unrelaxed) Newton correction              |
| nFcn*    | fcn evaluation count                                            |

(* marks inout parameters)

## Output parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| a[lda,n] | Array containing the approximated Jacobian                      |
| eta[n]*  | Vector of scaled denominator differences adjusted               |
| nFcn*    | fcn evaluation count adjusted                                   |
| iFail    | Return code non-zero if Jacobian could not be computed          |
"""
function nJcf(fcn, n::Int64, lda::Int64, x::Vector{Float64}, fx::Vector{Float64},
    yscal::Vector{Float64}, eta::Vector{Float64}, etamin::Float64, etamax::Float64,
    etadif::Float64, conv::Float64, nFcn::Int64, a::Array{Float64,2})
    # Constant
    small2 = 0.1
    # Begin
    iFail = 0
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
                fcn(fu,x)
            catch
                iFail = -1
            end
            nFcn += 1
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
    return (nFcn,iFail)
end

"""
function nJcfb(fcn, n::Int64, lda::Int64, ml::Int64, x::Vector{Float64},
    fx::Vector{Float64}, yscal::Vector{Float64}, eta::Vector{Float64},
    etamin::Float64, etamax::Float64, etadif::Float64, conv::Float64,
    nFcn::Int64, a::Array{Float64,2})

Approximation of banded Jacobian matrix for nonlinear systems solver with
feed-back control of discretization and rounding errors

## Input parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| fcn      | Function of the form fcn(f, x) to provide right-hand side       |
| n        | Number of rows and columns of the Jacobian                      |
| lda      | Leading dimension of the array "a"                              |
| ml       | Lower bandwidth of the Jacobian matrix                          |
| x[n]     | Current scaled vector                                           |
| fx[n]    | Vector containing fcn(x)                                        |
| yscal[n] | Vector containing scaling factors                               |
| eta[n]*  | Vector of scaled denominator differences                        |
| etamin   | Minimum allowed scaled denominator                              |
| etamax   | Maximum allowed scaled denominator                              |
| etadif   | = sqrt(1.1*epMach)                                              |
| conv     | Maximum norm of last (unrelaxed) Newton correction              |
| nFcn*    | fcn evaluation count                                            |

(* marks inout parameters)

## Output parameters
-------------------
| Variable | Description                                                     |
|----------|-----------------------------------------------------------------|
| a[lda,n] | Array containing the approximated Jacobian                      |
| eta[n]*  | Vector of scaled denominator differences adjusted               |
| nFcn*    | fcn evaluation count adjusted                                   |
| iFail    | Return code non-zero if Jacobian could not be computed          |
"""
function nJcfb(fcn, n::Int64, lda::Int64, ml::Int64, x::Vector{Float64},
    fx::Vector{Float64}, yscal::Vector{Float64}, eta::Vector{Float64},
    etamin::Float64, etamax::Float64, etadif::Float64, conv::Float64,
    nFcn::Int64, a::Array{Float64,2})
    # Constants
    small2 = 0.1
    # Begin
    mu = lda - 2*ml - 1
    ldab = ml + mu + 1
    w = zeros(n)
    u = zeros(n)
    iFail = 0
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
                fcn(fu,x)
            catch
                iFail = -1
            end
            nFcn += 1
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
    return (nFcn,iFail)
end
