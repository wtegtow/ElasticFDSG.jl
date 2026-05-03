"""
    config_template_2d(; device, precision, fd_order, verbose, output_file,
                         t_start, t_end, dt,
                         fdom, wavelet, wavelet_center, seismic_moment,
                         src_x, src_z, Mxx, Mxz, Mzz, anisotropic,
                         xstart, xend, zstart, zend, pml_layer,
                         geophones, das_x_aligned, das_z_aligned,
                         snapshot_times, snapshot_fields) -> Dict

Build a 2D simulation configuration dictionary.

All keyword arguments are required. The returned `Dict` can be passed directly
to [`runsim`](@ref) or serialised to YAML by the user.

```
"""
function config_template_2d(;
    device::String,
    precision::String,
    fd_order::Int,
    verbose::Bool,
    output_file,
    t_start::Real,
    t_end::Real,
    dt::Real,
    fdom::Real,
    wavelet::String,
    wavelet_center::Real,
    seismic_moment::Real,
    src_x::Real,
    src_z::Real,
    Mxx::Real,
    Mxz::Real,
    Mzz::Real,
    anisotropic::Bool,
    xstart::String,
    xend::String,
    zstart::String,
    zend::String,
    pml_layer::Int,
    geophones,
    das_x_aligned,
    das_z_aligned,
    snapshot_times,
    snapshot_fields,
)
    # validate enum fields (mirrors _REQUIRED_SETTINGS / _REQUIRED_SOURCE_2D)
    precision  in ["Float32", "Float64"] ||
        error("config_template_2d: precision must be \"Float32\" or \"Float64\"")
    wavelet    in ["ricker", "gauss1d"]  ||
        error("config_template_2d: wavelet must be \"ricker\" or \"gauss1d\"")
    for (side, val) in [("xstart",xstart),("xend",xend),("zstart",zstart),("zend",zend)]
        val in ["absorbing", "free", "none"] ||
            error("config_template_2d: boundary \"$side\" must be \"absorbing\", \"free\", or \"none\"")
    end

    return Dict(
        "settings" => Dict(
            "device"                   => device,
            "precision"                => precision,
            "spatial_derivative_order" => fd_order,
            "verbose"                  => verbose,
            "output_file"              => output_file,
        ),
        "time" => Dict(
            "start"    => t_start,
            "end"      => t_end,
            "timestep" => dt,
        ),
        "source" => Dict(
            "dominant_frequency" => fdom,
            "wavelet_type"       => wavelet,
            "wavelet_center"     => wavelet_center,
            "seismic_moment"     => seismic_moment,
            "location"           => Dict(
                "x" => src_x,
                "z" => src_z,
            ),
            "moment_tensor"      => Dict(
                "Mxx"        => Mxx,
                "Mxz"        => Mxz,
                "Mzz"        => Mzz,
                "anisotropic" => anisotropic,
            ),
        ),
        "boundaries" => Dict(
            "xstart"    => xstart,
            "xend"      => xend,
            "zstart"    => zstart,
            "zend"      => zend,
            "pml_layer" => pml_layer,
        ),
        "receivers" => Dict(
            "geophones" => geophones,
            "das"       => Dict(
                "x_aligned" => das_x_aligned,
                "z_aligned" => das_z_aligned,
            ),
            "snapshots" => Dict(
                "times"  => snapshot_times,
                "fields" => snapshot_fields,
            ),
        ),
    )
end


"""
    config_template_3d(; device, precision, fd_order, verbose, output_file,
                         t_start, t_end, dt,
                         fdom, wavelet, wavelet_center, seismic_moment,
                         src_x, src_y, src_z,
                         Mxx, Mxy, Mxz, Myy, Myz, Mzz, anisotropic,
                         xstart, xend, ystart, yend, zstart, zend, pml_layer,
                         geophones, das_x_aligned, das_y_aligned, das_z_aligned,
                         snapshot_positions, snapshot_times, snapshot_fields) -> Dict

Build a 3D simulation configuration dictionary.

All keyword arguments are required. The returned `Dict` can be passed directly
to [`runsim`](@ref) or serialised to YAML by the user.

```
"""
function config_template_3d(;
    device::String,
    precision::String,
    fd_order::Int,
    verbose::Bool,
    output_file,
    t_start::Real,
    t_end::Real,
    dt::Real,
    fdom::Real,
    wavelet::String,
    wavelet_center::Real,
    seismic_moment::Real,
    src_x::Real,
    src_y::Real,
    src_z::Real,
    Mxx::Real,
    Mxy::Real,
    Mxz::Real,
    Myy::Real,
    Myz::Real,
    Mzz::Real,
    anisotropic::Bool,
    xstart::String,
    xend::String,
    ystart::String,
    yend::String,
    zstart::String,
    zend::String,
    pml_layer::Int,
    geophones,
    das_x_aligned,
    das_y_aligned,
    das_z_aligned,
    snapshot_positions,
    snapshot_times,
    snapshot_fields,
)
    precision in ["Float32", "Float64"] ||
        error("config_template_3d: precision must be \"Float32\" or \"Float64\"")
    wavelet in ["ricker", "gauss1d"] ||
        error("config_template_3d: wavelet must be \"ricker\" or \"gauss1d\"")
    for (side, val) in [("xstart",xstart),("xend",xend),("ystart",ystart),
                        ("yend",yend),("zstart",zstart),("zend",zend)]
        val in ["absorbing", "free", "none"] ||
            error("config_template_3d: boundary \"$side\" must be \"absorbing\", \"free\", or \"none\"")
    end

    return Dict(
        "settings" => Dict(
            "device"                   => device,
            "precision"                => precision,
            "spatial_derivative_order" => fd_order,
            "verbose"                  => verbose,
            "output_file"              => output_file,
        ),
        "time" => Dict(
            "start"    => t_start,
            "end"      => t_end,
            "timestep" => dt,
        ),
        "source" => Dict(
            "dominant_frequency" => fdom,
            "wavelet_type"       => wavelet,
            "wavelet_center"     => wavelet_center,
            "seismic_moment"     => seismic_moment,
            "location"           => Dict(
                "x" => src_x,
                "y" => src_y,
                "z" => src_z,
            ),
            "moment_tensor"      => Dict(
                "Mxx"        => Mxx,
                "Mxy"        => Mxy,
                "Mxz"        => Mxz,
                "Myy"        => Myy,
                "Myz"        => Myz,
                "Mzz"        => Mzz,
                "anisotropic" => anisotropic,
            ),
        ),
        "boundaries" => Dict(
            "xstart"    => xstart,
            "xend"      => xend,
            "ystart"    => ystart,
            "yend"      => yend,
            "zstart"    => zstart,
            "zend"      => zend,
            "pml_layer" => pml_layer,
        ),
        "receivers" => Dict(
            "geophones" => geophones,
            "das"       => Dict(
                "x_aligned" => das_x_aligned,
                "y_aligned" => das_y_aligned,
                "z_aligned" => das_z_aligned,
            ),
            "snapshots" => Dict(
                "plane_positions" => snapshot_positions,
                "times"           => snapshot_times,
                "fields"          => snapshot_fields,
            ),
        ),
    )
end