function chebyQuad(x,f)
    n = length(x)
    for i = 2:2:n
        i1 = i - 1
        f[i1] = 0.0
        f[i]  = n/(i^2-1)
    end
    if n%2 == 1
        f[n] = 0.0
    end
    for l = 1:n
        factt = 4.0*x[l] - 2.0
        ti2 = 1.0
        ti1 = 0.5*factt
        f[1] = ti1 + f[1]
        for i = 2:n
            ti = factt*ti1-ti2
            f[i] = ti + f[i]
            ti2 = ti1
            ti1 = ti
        end
    end
    return nothing
end

function chebyQuadJac(x,J)
    n = length(x)
    for j = 1:n
        factt = 4.0*x[j] - 2.0
        ti2 = 1.0
        ti1 = 0.5*factt
        tabli1 = 0.0
        tabli2 = 2.0
        J[1,j] = tabli1
        for i = 2:n
            ti = factt*ti1-ti2
            tabli = 4.0*ti1+factt*tabli1-tabli2
            J[i,j] = tabli
            ti2 = ti1
            ti1 = ti
            tabli2 = tabli1
            tabli1 = tabli
        end
    end
    return nothing
end
