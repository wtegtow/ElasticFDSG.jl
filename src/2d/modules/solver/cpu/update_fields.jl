# =============================================================================
# CPU update fields
# =============================================================================

function update_normal_stresses!(settings::Settings,
                                domain::Domain,
                                elastic::Elastic,
                                fields::Fields,
                                pml::Pml,
                                time::Time)


    Threads.@threads for y in 2+(settings.N-1):domain.ny-(settings.N-1)
    for x in 1+(settings.N-1):domain.nx-settings.N

        # EFFECTIVE ELASTIC PARAMETER
        c11eff = 0.5 * (elastic.c11[y,x+1] + elastic.c11[y,x])
        c13eff = 0.5 * (elastic.c13[y,x+1] + elastic.c13[y,x])
        c33eff = 0.5 * (elastic.c33[y,x+1] + elastic.c33[y,x])

        # SPACE DERIVATIVE
        vx_x = 0
        vy_y = 0

        for i in 1:settings.N
            vx_x += settings.c[i]/domain.dx * (fields.vx[y,x+i]-fields.vx[y,x-(i-1)])
            vy_y += settings.c[i]/domain.dy * (fields.vy[y+(i-1),x]-fields.vy[y-i,x])
        end

        # UPDATE PML REGION
        if domain.pml_points_hash[y,x] != -1 

            pml_idx = domain.pml_points_hash[y,x]


            pml.vx_x_old[pml_idx] =  pml.b_x_odd[x] *  pml.vx_x_old[pml_idx] +  pml.a_x_odd[x] * vx_x
            pml.vy_y_old[pml_idx] =  pml.b_y_evn[y] *  pml.vy_y_old[pml_idx] +  pml.a_y_evn[y] * vy_y

            # UPDATE SPACE DERIVATIVE
            vx_x = vx_x/pml.K_x_odd[x] + pml.vx_x_old[pml_idx]
            vy_y = vy_y/pml.K_y_evn[y] + pml.vy_y_old[pml_idx]

        end 

        # UPDATE STRESSES
        fields.sxx[y,x] += time.dt * (c11eff * vx_x +  c13eff * vy_y)
        fields.syy[y,x] += time.dt * (c13eff * vx_x +  c33eff * vy_y)

    end
    end

end;



function update_sxy!(settings::Settings,
                    domain::Domain,
                    elastic::Elastic,
                    fields::Fields,
                    pml::Pml,
                    time::Time)

    Threads.@threads for y in 1+(settings.N-1):domain.ny-settings.N
    for x in 2+(settings.N-1):domain.nx-(settings.N-1)

        # EFFECTIVE ELASTIC PARAMETER
        c55eff = 0.5 * (elastic.c55[y+1,x] + elastic.c55[y,x])

        # SPACE DERIVATIVE
        vy_x = 0
        vx_y = 0
        for i in 1:settings.N
            vy_x +=  settings.c[i]/domain.dx * (fields.vy[y,x+(i-1)]-fields.vy[y,x-i])
            vx_y +=  settings.c[i]/domain.dy * (fields.vx[y+i,x]-fields.vx[y-(i-1),x])
        end

        # UPDATE PML REGION
        if domain.pml_points_hash[y,x] != -1 

            pml_idx = domain.pml_points_hash[y,x]

            pml.vy_x_old[pml_idx] = pml.b_x_evn[x] * pml.vy_x_old[pml_idx] + pml.a_x_evn[x] * vy_x
            pml.vx_y_old[pml_idx] = pml.b_y_odd[y] * pml.vx_y_old[pml_idx] + pml.a_y_odd[y] * vx_y

            # UPDATE SPACE DERIVATIVE
            vy_x = vy_x/pml.K_x_evn[x] + pml.vy_x_old[pml_idx]
            vx_y = vx_y/pml.K_y_odd[y] + pml.vx_y_old[pml_idx]
        end

        # UPDATE STRESSES
        fields.sxy[y,x] += time.dt * c55eff * (vy_x + vx_y)

    end
    end
end




function update_vx!(settings::Settings,
                    domain::Domain,
                    elastic::Elastic,
                    fields::Fields,
                    pml::Pml,
                    time::Time)

    Threads.@threads for y in 2+(settings.N-1):domain.ny-(settings.N-1)
    for x in 2+(settings.N-1):domain.nx-(settings.N-1)

        # EFFECTIVE ELASTIC PARAMETER
        # rho is defined on vx-grid points

        # SPACE DERIVATIVE
        sxx_x = 0
        sxy_y = 0
        for i in 1:settings.N 
            sxx_x += settings.c[i]/domain.dx * (fields.sxx[y,x+(i-1)]-fields.sxx[y,x-i])
            sxy_y += settings.c[i]/domain.dy * (fields.sxy[y+(i-1),x]-fields.sxy[y-i,x])
        end

        # UPDATE PML REGION
        if domain.pml_points_hash[y,x] != -1 

            pml_idx = domain.pml_points_hash[y,x]

            pml.sxx_x_old[pml_idx] = pml.b_x_evn[x] * pml.sxx_x_old[pml_idx] + pml.a_x_evn[x] * sxx_x
            pml.sxy_y_old[pml_idx] = pml.b_y_evn[y] * pml.sxy_y_old[pml_idx] + pml.a_y_evn[y] * sxy_y

            # UPDATE SPACE DERIVATIVE
            sxx_x = sxx_x/pml.K_x_evn[x] + pml.sxx_x_old[pml_idx]
            sxy_y = sxy_y/pml.K_y_evn[y] + pml.sxy_y_old[pml_idx]
        end

        # UPDATE VELOCITY
        fields.vx[y,x] += time.dt/elastic.rho[y,x] * (sxx_x + sxy_y)

    end
    end
end



function update_vy!(settings::Settings,
    domain::Domain,
    elastic::Elastic,
    fields::Fields,
    pml::Pml,
    time::Time)

    Threads.@threads for y in 1+(settings.N-1):domain.ny-settings.N
    for x in 1+(settings.N-1):domain.nx-settings.N

        # EFFECTIVE ELASTIC PARAMETER
        rhoeff = 0.25 * (elastic.rho[y+1,x] + elastic.rho[y,x] + elastic.rho[y,x+1] + elastic.rho[y+1,x+1])

        # SPACE DERIVATIVE
        sxy_x = 0
        syy_y = 0
        for i in 1:settings.N 
            sxy_x += settings.c[i]/domain.dx * (fields.sxy[y,x+i] - fields.sxy[y,x-(i-1)])
            syy_y += settings.c[i]/domain.dy * (fields.syy[y+i,x] - fields.syy[y-(i-1),x])
        end

        # UPDATE PML REGION
        if domain.pml_points_hash[y,x] != -1 

            pml_idx = domain.pml_points_hash[y,x]

            pml.sxy_x_old[pml_idx] = pml.b_x_odd[x] * pml.sxy_x_old[pml_idx] + pml.a_x_odd[x] * sxy_x
            pml.syy_y_old[pml_idx] = pml.b_y_odd[y] * pml.syy_y_old[pml_idx] + pml.a_y_odd[y] * syy_y

            # UPDATE SPACE DERIVATIVE
            sxy_x = sxy_x/pml.K_x_odd[x] + pml.sxy_x_old[pml_idx]
            syy_y = syy_y/pml.K_y_odd[y] + pml.syy_y_old[pml_idx]

        end

        # UPDATE VELOCITY
        fields.vy[y,x] += time.dt/rhoeff * (sxy_x + syy_y)

    end
    end

end;