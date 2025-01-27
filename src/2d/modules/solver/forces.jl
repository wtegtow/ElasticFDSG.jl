function get_forces(fdsg2d::FDSG2D)

    force_type = fdsg2d.settings.config["source"]["point_source"] == true ? "point" : "double"
    force_angle = deg2rad(fdsg2d.settings.config["source"]["force_angle"])

    if force_type == "point"

        force_x = sin(force_angle) .* fdsg2d.source.STF .* fdsg2d.time.dt/fdsg2d.source.rhosrc;
        force_y = cos(force_angle) .* fdsg2d.source.STF .* fdsg2d.time.dt/fdsg2d.source.rhosrc;

        force_x = fdsg2d.settings.float.(force_x)
        force_y = fdsg2d.settings.float.(force_y)

        return [force_x,  force_y]

    elseif force_type == "double"

        Mxx = sin(force_angle)^2 
        Mxy = sin(force_angle) * cos(force_angle) 
        Myy = cos(force_angle)^2

        return fdsg2d.settings.float.([Mxx, Mxy, Myy] .* fdsg2d.time.dt /  fdsg2d.source.rhosrc)

    end
end

function apply_forces!(forces, vx, vy, yid, xid, dx, dy, STF, ti) 

    if length(forces) == 2
        force_x, force_y = forces
        vx[yid, xid] = force_x[ti]
        vy[yid, xid] = force_y[ti]
     
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
