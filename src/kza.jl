"""
source: https://github.com/cran/kza/blob/master/src/kzf.c#L13
"""
function R_differenced(y::Vector{T}, q::Int) where T
    n = length(y)
    d = Vector{T}(undef, n)
    dprime = Vector{T}(undef, n)
    
	# calculate d = |Z(i+q) - Z(i-q)|
    #R: for (i=0; i<q; i++) {REAL(d)[i] = fabs(REAL(y)[i+q] - REAL(y)[0]);}
	for i in 1:q
        try
            d[i] = abs(y[i+q] - y[1])
        catch
            @info i q
            throw("lol")
        end
    end
    #R: for (i=q; i<n-q; i++) {REAL(d)[i] = fabs(REAL(y)[i+q] - REAL(y)[i-q]);}
	for i in q+1:n-q
        d[i] = abs(y[i+q] - y[i-q])
    end
    #R: for (i=n-q; i<n; i++) {REAL(d)[i] = fabs(REAL(y)[n-1] - REAL(y)[i-q]);}
	for i in n-q+1:n
        d[i] = abs(y[n] - y[i-q])
    end

	# d'(t) = d(i+1)-d(i)
    #R: for(i=0; i<n-1; i++) REAL(dprime)[i] = REAL(d)[i+1]-REAL(d)[i];
	for i in 1:n-1
        dprime[i] = d[i+1]-d[i]
    end
	dprime[n] = dprime[n-1]
    return d, dprime
end


"""
compute moving average value for array `v` at index `col` with window `w`
skipping any `missing` values in `v`
"""
function mavg1d(v::Vector{T}, col::Int, w::Int) where T
    s::T = 0.0
    z = 0

    startcol = col-w>1 ? col-w : 1
    endcol = col+w<length(v) ? col+w : length(v)

    #@show col w startcol endcol
    
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
adaptive curvature function
"""
function adaptive(d::T, m::T) where T
	return 1 - d/m
end


"""
Apply Adaptive Kolmogorov-Zurbenko filter to 1D Vector `x` over window of width
`window` using `iterations` number of iterations

v: A vector of the time series
window: The window for the filter.
iterations: The number of iterations.
min_size:  Minimum size of window q.
tolerance: The smallest value to accept as nonzero.

source: https://github.com/cran/kza/blob/master/src/kza.c#L60
"""
function kza(v::Vector{T}, window::Integer; iterations=3, minimum_window_length=Int(round(0.05*window)), tolerance=1.0e-5) where T
    n = length(v)
    eps = tolerance
	q = window
	min_window_length = minimum_window_length

    y = kz(v, q; iterations=iterations)
    d, dprime = R_differenced(y, q)
 
    m = maximum(skipmissing(d))

    tmp = copy(v)
    ans = Vector{T}(undef, n)

    #for(i=0; i<INTEGER_VALUE(iterations); i++) {
    for i in 1:iterations
    	#for (t=0; t<n; t++) {
        for t in 1:n
		    #if (fabs(REAL(dprime)[t]) < eps) { /* dprime[t] = 0 */
            if abs(dprime[t]) < eps  # dprime[t] = 0
			    qh = Int(floor(q*adaptive(d[t], m)))
			    qt = Int(floor(q*adaptive(d[t], m)))
            elseif dprime[t] < 0.0
                qh = q
    		    qt = Int(floor(q*adaptive(d[t], m)))
            else
		    	qh = Int(floor(q*adaptive(d[t], m)))
			    qt = q;
            end
			qt = (qt < min_window_length) ? min_window_length : qt
			qh = (qh < min_window_length) ? min_window_length : qh
            
            #@info "step" d[t] adaptive(d[t], m) q t qt qh
	
	        # /* check bounds */
        	qh = (qh > n-t) ? n-t : qh; # head past end of serie
            qt = (qt >= t) ? t-1 : qt;  	        		
            #@info "bounds" t-qt t+qh t
   		    ans[t] = mean(skipmissing(tmp[t-qt:t+qh]))
        end
	    # copyVector(tmp, ans);
        tmp[:] = ans[:]
    end
	return ans
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