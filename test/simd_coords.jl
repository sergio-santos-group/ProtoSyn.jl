
module LATTICE
    @enum TYPE begin
        primitive       = 1
        body_centered   = 2
        face_centered   = 3
    end
end

function generate_template(hf::Vector{Float64}, type::LATTICE.TYPE)::Array{Float64, 2}
    atoms = [0.0 0.0 0.0]
    if type == LATTICE.body_centered
        atoms = vcat(atoms, [hf[1] hf[2] hf[3]])
    elseif type == LATTICE.face_centered
        atoms = vcat(atoms, [hf[1] hf[2] 0.0; hf[1] 0.0 hf[3]; 0.0 hf[2] hf[3]])
    end
    return atoms
end


function calculate_n_atoms(rep::Vector{Int64}, type::LATTICE.TYPE,
    closed::Bool = true)::Int64

    if type == LATTICE.primitive
        n_atoms = (rep[1] + closed) * (rep[2] + closed) * (rep[3] + closed)
    elseif type == LATTICE.body_centered
        n_atoms = (rep[1] + closed) * (rep[2] + closed) * (rep[3] + closed)
        n_atoms += rep[1]*rep[2]*rep[3]
    elseif type == LATTICE.face_centered
        n_atoms = (rep[1] + closed) * (rep[2] + closed) * (rep[3] + closed)
        n_atoms += rep[1]*rep[2]*rep[3] * 3
        if closed
            n_atoms += rep[1]*rep[2] + rep[2] * rep[3] + rep[3] * rep[1]
        end
    end
    return n_atoms
end


function generate_cubic_lattice(size::Float64, rep::Int64,
    type::LATTICE.TYPE, closed::Bool = true)
    return generate_lattice([size, size, size], [rep, rep, rep], type, closed)
end


function generate_lattice(side_length::Vector{Float64}, rep::Vector{Int64},
    type::LATTICE.TYPE, closed::Bool = true)
    # Side length in nm
    # Generate the template
    f = (side_length./rep)
    template = generate_template(f./2, type)
    # Create empty list of atoms
    n_atoms  = calculate_n_atoms(rep, type, closed)
    xyz = zeros(4, n_atoms)
    # Fill State with repetitions of template
    c = 1
    for i in 0:rep[1]
        for j in 0:rep[2]
            for k in 0:rep[3]
                for index in 1:size(template, 1)
                    new_pos = template[index, :] .+ [i*f[1], j*f[2], k*f[3]]
                    if closed && all(new_pos .<= side_length) || all(new_pos .< side_length)
                        # state.xyz[c, :] = new_pos'
                        xyz[1:3, c] .= new_pos
                        c +=1
                    end
                end
            end
        end
    end
    return xyz
end


# @time x = generate_cubic_lattice(10.0, 31, LATTICE.face_centered, false)
# println(size(x))
# println(x[:,1:10])

# open("tmp.xyz", "w") do fout
#     n = size(x,2)
#     x .*= 10.0
#     println(fout, "$n")
#     println(fout, "treta")
#     for i=1:n
#         println(fout, "X  $(x[1,i])  $(x[2,i])  $(x[3,i])")
#     end
# end