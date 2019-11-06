macro force_switch_coefA(α::Int)
    esc(:(α * ($(α+4)*rc - $(α+1)*r1) / (rc^$(α+2) * (rc-r1)^2)))
end

macro force_switch_coefB(α::Int)
    esc(:(α * ($(α+3)*rc - $(α+1)*r1) / (rc^$(α+2) * (rc-r1)^3)))
end

macro nonbonded(args...)
    options = get_options(args)

ex = quote
function(atoms::Vector{Atom}, state::Common.State, do_forces=true)

    $(if get(options, :mode, :none) == :force_switch && get(options, :include_lj, true)
        quote
            fsw = cutoff - switch_r1   # force switch width
            fsA6  = @force_switch_coefA 6
            fsB6  = @force_switch_coefB 6
            fsA12 = @force_switch_coefA 12
            fsB12 = @force_switch_coefB 12
            fsC6  = 1.0/(cutoff^12) - fsA6 *(fsw^3)/3.0 - fsB6 *(fsw^4)/4.0
            fsC12 = 1.0/(cutoff^12) - fsA12*(fsw^3)/3.0 - fsB12*(fsw^4)/4.0
        end
    end)

    $(if get(options, :mode, :none) == :force_switch && get(options, :include_coulomb, true)
        quote
            fsw = cutoff - switch_r1   # force switch width
            fsA2 = @force_switch_coefA 2
            fsB2 = @force_switch_coefB 2
            fsC2 = 1.0/(cutoff^12) - fsA2*(fsw^3)/3.0 - fsB2*(fsw^4)/4.0
        end
    end)

    for i=1:(n_atoms-1)
        atomi = atoms[i]
        ptr = state.nbptr[i]
        while state.nblist[ptr] > 0
            @inbounds j = state.nblist[ptr]
            @inbounds atomj = atoms[j]
            ptr += 1

            dijSq = .0
            @inbounds for k=1:3
                deltaij = coords[j, k] - coords[i, k]
                dijSq  += deltaij * deltaij
                vij[k]  = deltaij
            end

            if dij_sq > cutoff_sq
                continue
            end

            # =============== START CONFIGURABLE SECTION ================
            #                 LENNARD-JONES
            $(if get(options, :include_lj, true)
                if get(options, :mode, :cut) == :cut
                    quote
                        sij = atomi.σ + atomj.σ
                        eij = atomi.ϵ * atomj.ϵ
                        lj6 = (sij * sij/dijSq) ^ 3
                        eLJ += eij * (lj6 * lj6 - lj6)
                    end
                elseif get(options, :mode) == :force_switch
                    quote
                        if dij_sq < switch_cut1_sq
                            sij = atomi.σ + atomj.σ
                            eij = atomi.ϵ * atomj.ϵ
                            lj6 = (sij * sij/dijSq) ^ 3
                            eLJ += eij * (lj6 * lj6 - lj6)
                        else
                            dij = sqrt(dijSq)
                            dr = dij - r1

                            1.0/(dij^6)  - fsA6*( dr^3)/3.0 - fsB6*( dr^4)/4.0 - fsC6
                            1.0/(dij^12) - fsA12*(dr^3)/3.0 - fsB12*(dr^4)/4.0 - fsC12

                        end

                        

                        eLJ += 1.0/
                    end
                end

            end)

            #                 COULOMB
            $(if get(options, :include_coulomb, true)
                quote
                    ecoul = atomi.q * atomj.q / sqrt(dijSq)
                    eCoulomb += ecoul
                end
            end)
            # ================ END CONFIGURABLE SECTION =================
            
            # =============== START CONFIGURABLE SECTION ================
            $(

            )
            # ================ END CONFIGURABLE SECTION =================

        end

    end

end # end fcn declaration
end # end quote
    ex
end # end macro
