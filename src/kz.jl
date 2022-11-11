"""
compute moving average value for Vector `v` at index `col` with window `w`
skipping any `missing` values in `v`
"""
function mavg1d(v::Vector{T}, col::Int, w::Int) where T
    s::T = 0.0
    z = 0

    startcol = col-w>1 ? col-w : 1
    endcol = col+w<length(v) ? col+w : length(v)
    
    for i in startcol:endcol
        if !ismissing(v[i])
            z += 1
            s += v[i]
        end
    end
    if (z == 0)
        return missing
    else
        return s/z
    end
end


"""
Apply Kolmogorov-Zurbenko filter to Vector `x` with `window` over `iterations`
"""
function kz(x::Vector{T}, window::Int; iterations::Int=3) where T
    ans = Vector{T}(undef, length(x))
    tmp = copy(x)
    
    for k in 1:iterations
        for i in 1:length(x)
            ans[i] = mavg1d(tmp, i, window)
        end
        tmp[:] = ans[:]
    end
    return ans
end