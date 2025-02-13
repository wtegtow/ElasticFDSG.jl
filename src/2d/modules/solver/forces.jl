function get_forces(fdsg2d::FDSG2D)

    use_point = fdsg2d.settings.config["source"]["point_source"]["use"]
    use_double = fdsg2d.settings.config["source"]["double_couple"]["use"]
    @assert use_point != use_double "Only one source type can be used at a time (point or double couple)"

    if use_point

        force_x = fdsg2d.settings.float.(zeros(fdsg2d.time.nt))
        force_y = fdsg2d.settings.float.(zeros(fdsg2d.time.nt))

        # only on vx (scaled by dx,dy for analytical solution)
        if fdsg2d.settings.config["source"]["point_source"]["act_on"]["on_vx"] == true 
            force_x = fdsg2d.source.STF .* fdsg2d.time.dt/(fdsg2d.source.rhosrc * fdsg2d.domain.dx * fdsg2d.domain.dy);
        end 

        # only on vy (scaled by dx,dy for analytical solution)
        if fdsg2d.settings.config["source"]["point_source"]["act_on"]["on_vy"] == true 
            force_y = fdsg2d.source.STF .* fdsg2d.time.dt/(fdsg2d.source.rhosrc * fdsg2d.domain.dx * fdsg2d.domain.dy);
        end 

        # directed force acting on vx and vy
        if fdsg2d.settings.config["source"]["point_source"]["act_on"]["on_vx"] == true  &&
           fdsg2d.settings.config["source"]["point_source"]["act_on"]["on_vy"] == true 

           force_angle = fdsg2d.settings.config["source"]["point_source"]["act_on"]["on_vx_and_vy"]["force_angle"]
           force_angle = deg2rad(force_angle)

           force_x = sin(force_angle) .* fdsg2d.source.STF .* fdsg2d.time.dt/fdsg2d.source.rhosrc;
           force_y = cos(force_angle) .* fdsg2d.source.STF .* fdsg2d.time.dt/fdsg2d.source.rhosrc;

        end 

        return [force_x,  force_y]

    elseif use_double

        strike = fdsg2d.settings.config["source"]["double_couple"]["strike"]
        strike = deg2rad(strike)

        Mxx = sin(strike)^2 - cos(strike)^2
        Mxy = 2 * cos(strike) * sin(strike)
        Myy = cos(strike)^2 - sin(strike)^2

        return fdsg2d.settings.float.([Mxx, Mxy, Myy] .* fdsg2d.time.dt / fdsg2d.source.rhosrc)

    end
end


function apply_forces!(forces, vx, vy, yid, xid, dx, dy, STF, ti) 

    # point source
    if length(forces) == 2
        force_x, force_y = forces
        vx[yid, xid] += force_x[ti]
        vy[yid, xid] += force_y[ti]
    
    # double couple
    elseif length(forces) == 3
    
        Mxx, Mxy, Myy = forces

        # Force X
        # Mxx
        vx[yid, xid+1] += (Mxx * STF[ti]) / (dx^2 * dy)
        vx[yid, xid]   -= (Mxx * STF[ti]) / (dx^2 * dy)
    
        # Mxy
        vx[yid+1, xid+1] += (Mxy * STF[ti]) / (4 * dx * dy^2)
        vx[yid+1, xid]   += (Mxy * STF[ti]) / (4 * dx * dy^2)
        vx[yid-1, xid+1] -= (Mxy * STF[ti]) / (4 * dx * dy^2)
        vx[yid-1, xid]   -= (Mxy * STF[ti]) / (4 * dx * dy^2)
    
        # Force Y
        # Myy
        vy[yid+1, xid] += (Myy * STF[ti]) / (dx * dy^2)
        vy[yid, xid]   -= (Myy * STF[ti]) / (dx * dy^2)
    
        # Mxy
        vy[yid+1, xid+1] += (Mxy * STF[ti]) / (dx^2 * 4 * dy)
        vy[yid, xid+1]   += (Mxy * STF[ti]) / (dx^2 * 4 * dy)
        vy[yid+1, xid-1] -= (Mxy * STF[ti]) / (dx^2 * 4 * dy)
        vy[yid, xid-1]   -= (Mxy * STF[ti]) / (dx^2 * 4 * dy)
        end
end