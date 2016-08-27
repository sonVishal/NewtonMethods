function chebyQuad(f,x)
    n = length(x)
    f[:] = zeros(n);
    for i = 2:2:n
        f[i] = n/(i*i-1)
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

function chebyQuadJac(J,x)
    n = length(x)
    for j = 1:n
        factt = 4.0*x[j] - 2.0
        ti2 = 1.0
        ti1 = 0.5*factt
        tabli2 = 0.0
        tabli1 = 2.0
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
