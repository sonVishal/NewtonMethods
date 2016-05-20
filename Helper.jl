function n1scal(n,x,xa,xScal,iScal,qIniSc,opt)
    small = getMachineConstants(6)
    # Begin
    if iScal == 1
        xw = xScal
    else
        xw = zeros(n)
        for l1 = 1:n
            xw[l1] = max(xScal[l1],max((abs(x[l1])+abs(xa[l1])*0.5),small))
        end
    end

    mPr = opt.options[OPT_PRINTITERATION]
    if mPr >= 6
        printIO = opt.options[OPT_PRINTIO]
        write(printIO,"\n\n",
        "+++++++++++++++++++++++++++++++++++++++++++++++++\n",
        "      x-components         scaling-components\n")
        for l1 = 1:n
            write(printIO,"  %18.10e   %18.10e\n",x[l1],xw[l1])
        end
        write(printIO,"+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
    end
    return xw
end

# TODO: There is a smarter way to write scaling using in built scale!(A,b) (col)
# and scale!(b,A) (row) functions
function n1scrf(m,n,a)
    # Begin
    fw = zeros(n)
    if issparse(a)
        nza = nnz(a)
        (row,col) = findn(a)
        aout = sparse(row,col,zeros(nza),m,n)
        for j = 1:nza
            rj = row[j]
            aaj = abs(a[rj,col[j]])
            fw[rj] = max(fw[rj],aaj)
        end
        for k = 1:m
            if fw[k] > 0.0
                fw[k] = 1.0/fw[k]
            else
                fw[k] = 1.0
            end
        end
        for j = 1:nza
            aout[row[j],col[j]] = a[row[j],col[j]]*fw[row[j]]
        end
    else
        aout = zeros(a)
        for k = 1:m
            s1 = max(abs(a[k,1:n]))
            if s1 > 0.0
                s1 = 1.0/s1
                aout[k,1:n] = a[k,1:n]*s1
            else
                fw[k] = 1.0
                aout[k,1:n] = a[k,1:n]
            end
        end
    end
    return (aout,fw)
end

function n1scrb(n,lda,ml,mu,a)
    # Begin
    fw = zeros(n)
    # TODO: zeros(a) or zeros(size(a)). The former creates a sparse if a is sparse
    aout = zeros(size(a))
    m2 = ml + mu + 1
    for k = 1:n
        s1 = 0.0
        l2 = max(1,k-ml)
        l3 = min(n,k+mu)
        k1 = m2 + k
        for l1 = l2:l3
            s1 = max(s1,abs(a[k1-l1,l1]))
        end
        if s1 > 0.0
            s1 = 1.0/s1
            fw[k] = s1
            for l1 = l2:l3
                aout[k1-l1,l1] = a[k1-l1,l1]*s1
            end
        else
            fw[k] = 1.0
        end
    end
    return (aout,fw)
end
