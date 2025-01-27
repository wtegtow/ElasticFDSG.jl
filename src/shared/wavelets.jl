function ricker(t, ts, fdom)  
    wavelet = @. (1 - (2*π^2 * fdom^2) * (t - ts)^2) * (exp(  (-π^2 * fdom^2 ) * (t - ts)^2 ) ) 
    wavelet ./= maximum(abs.(wavelet))
    return wavelet 
end

function gauss1d(t, ts, fdom) 
    wavelet = @. -2 * π^2 * fdom^2 * (t .- ts) * exp.(-π^2 * fdom^2 * (t - ts)^2)
    wavelet ./= maximum(abs.(wavelet))
    return wavelet 
end