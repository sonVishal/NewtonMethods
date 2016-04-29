function getMachineConstants(i)
# 1 = realmin
# 2 = realmax
# 3 = eps
# 4 = 2.0*eps
# 5 = log10(2.0)
# 6 = sqrt(realmin/eps)
    if i < 1 || i > 6
        error("getMachineConstants - i out of bounds");
    end

    if i == 1
        return realmin(Float64);
    elseif i == 2
        return realmax(Float64);
    elseif i == 3
        return eps(Float64);
    elseif i == 4
        return 2.0*eps(Float64);
    elseif i == 5
        return log10(2.0);
    elseif i == 6
        return sqrt(getMachineConstants(1)/getMachineConstants(3));
    end
end
