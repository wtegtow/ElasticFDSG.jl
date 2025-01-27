"""
_________________________________________________________________________
This file contains recycled code from: 
https://github.com/geodynamics/seismic_cpml
_________________________________________________________________________
Roland Martin and Dimitri Komatitsch and Stephen D. Gedney,
A variational formulation of a stabilized unsplit convolutional perfectly
matched layer for the isotropic or anisotropic seismic wave equation.
_________________________________________________________________________
"""

function cmpl(N, 
              npoints_pml, 
              use_pml_start,
              use_pml_end,
              domain_, 
              vmax, 
              fdom, 
              dt, 
              rcoef,
              FLOAT = Float64)
    
    # params
    K_max = 1
    α_max = π * fdom  
    NPOW = 3

    n = length(domain_)
    h = abs(domain_[2]-domain_[1])
    small_float = h / 1e3

    # output arrays
    K_evn  = ones(FLOAT,n);
    K_odd  = ones(FLOAT,n);
    a_evn  = zeros(FLOAT,n);
    a_odd  = zeros(FLOAT,n);
    b_evn  = zeros(FLOAT,n);
    b_odd  = zeros(FLOAT,n);
    
    # intermediate arrays
    α_evn = zeros(FLOAT,n)
    α_odd = zeros(FLOAT,n);
    d_evn = zeros(FLOAT,n);
    d_odd = zeros(FLOAT,n);

    # some usefull variables
    pml_width = npoints_pml * h
    ghost_width = N*h
    d0 = -(NPOW-1) * vmax * log(rcoef) / (2 * pml_width)

    # pml regions (start = left, end = right of domain)
    pml_start = domain_[begin] + ghost_width : h : domain_[begin] + ghost_width + pml_width - h
    pml_end = domain_[end] - ghost_width - pml_width + h : h : domain_[end] - ghost_width 
   
    # even grid points 
    for i in eachindex(domain_)

        val_i = domain_[i]
        val_in_pml = nothing
        in_pml = false

        if any(isapprox.(val_i, pml_start, atol=small_float)) && use_pml_start == true
            val_in_pml = pml_start[end] - val_i
            in_pml = true

        elseif any(isapprox.(val_i, pml_end, atol=small_float)) && use_pml_end == true
            val_in_pml = val_i - pml_end[1]
            in_pml = true
        end

        if  in_pml == true
            val_in_pml /= pml_width
            d_evn[i] = d0 * val_in_pml^NPOW
            K_evn[i] = 1 + (K_max - 1) * val_in_pml^NPOW
            α_evn[i] = α_max * (1 - val_in_pml)
        
        if α_evn[i] < 0
           α_evn[i] = 0
        end
        end
    end

    # odd grid points
    for i in eachindex(domain_)
        val_i = domain_[i]
        val_in_pml = nothing
        in_pml = false

        if any(isapprox.(val_i, pml_start, atol=small_float)) && use_pml_start == true
            val_in_pml = pml_start[end] - (val_i + h/2)
            in_pml = true

        elseif any(isapprox.(val_i, pml_end, atol=small_float)) && use_pml_end == true
            val_in_pml = (val_i + h/2) - pml_end[1]
            in_pml = true
        end

        if in_pml == true
            val_in_pml /= pml_width 
            d_odd[i] = d0 * val_in_pml^NPOW
            K_odd[i] = 1 + (K_max - 1) * val_in_pml^NPOW
            α_odd[i] = α_max * (1 - val_in_pml)

        if α_odd[i] < 0
            α_odd[i] = 0
        end
        end
    end

    # fill values 
    for i in eachindex(domain_)
        b_evn[i] = exp(- (d_evn[i] / K_evn[i] + α_evn[i]) * dt)
        b_odd[i] = exp(- (d_odd[i] / K_odd[i] + α_odd[i]) * dt)

        if d_evn[i] != 0
            a_evn[i] = d_evn[i] * (b_evn[i] - 1) / (K_evn[i] * (d_evn[i] + K_evn[i] * α_evn[i]))
        end

        if d_odd[i] != 0
            a_odd[i] = d_odd[i] * (b_odd[i] - 1) / (K_odd[i] * (d_odd[i] + K_odd[i] * α_odd[i]))
        end
    end    

    return K_evn, K_odd, a_evn, a_odd, b_evn , b_odd
end