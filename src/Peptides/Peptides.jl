module Peptides

using ..Common

#region ENUMS
baremodule DihedralTypes
using Base: @enum
@enum Type begin
    undef = 0
    phi   = 1
    psi   = 2
    omega = 3
    chi1  = 4
    chi2  = 5
    chi3  = 6
    chi4  = 7
    chi5  = 8
end
end



baremodule SecondaryStructureTypes
using Base: @enum
@enum Type begin
    undef = 0
    helix = 1
    sheet = 2
    coil  = 3
end
end
#endregion


mutable struct SecondaryStructureMetadata <: Common.AbstractMetadata
    index::Int
    name::String
    start::Common.Residue
    stop::Common.Residue
    type::Enum{<:Integer}
end


mutable struct BlockMetadata <: Common.AbstractMetadata
    index::Int
    name::String
    pivot::Int
    movable::Vector{Int}
    blocks::Vector{Common.Residue}
end


mutable struct CrankshaftMetadata <: Common.AbstractMetadata
    index
    name
    movable::Vector{Int}
    movable2
end

function treta(residues::Vector{Common.ResidueMetadata}, str::String)
    sec_md = SecondaryStructureMetadata[]
    re = r"(?P<helix>H+)|(?P<sheet>E+)"i
    for (index,m) in enumerate(eachmatch(re, str))
        if m[:helix] != nothing
            name = string("H", index)
            type = SecondaryStructureTypes.helix
        else
            name = string("S", index)
            type = SecondaryStructureTypes.sheet
        end
        r1 = residues[m.offset]
        r2 = residues[m.offset+length(m.match)-1]
        smd = SecondaryStructureMetadata(index, name, r1, r2, type)
        push!(sec_md, smd)
    end
    return sec_md
end


function finddihedrals(metadata::Common.Metadata)
    # list of all dihedrals
    dihedrals = Vector{Common.Dihedral}()

    # start by grouping all atoms by residues
    #   r1 => [at1,at2, ...]
    #   r2 => [atm,atn, ...]
    #   rn => [...]
    atsbyres = Dict{Common.Residue, Vector{Common.AtomMetadata}}()
    for atom in metadata.atoms
        push!(get!(atsbyres, atom.residue, []), atom)
    end

    # generate regex for identifying required atom names
    queries = [
        Regex("(.*):(?P<name>[A-GI-Z]$q[1]?):(.*)")
        for q in "BGDEZH"
    ]
    
    # list al residues, sorted by their index
    residues = sort(collect(keys(atsbyres)), by=r->r.index)
    for residue in residues
        
        # get atoms for this residue
        atoms = atsbyres[residue]

        # get a map from atom name to atom index
        name2index = Dict(at.name => at for at in atoms)
        
        #region BACKBONE-DIHEDRALS
        N  = name2index["N"]
        CA = name2index["CA"]
        C  = name2index["C"]
        
        # find ϕ (phi)
        if N.connects != nothing
            prev_C = findfirst(i->metadata.atoms[i].name=="C", N.connects)
            if prev_C != nothing
                type = DihedralTypes.phi
                mov = metadata.stree.indices[1:Common.getindex(metadata.stree, CA.index)]
                dihd = Common.Dihedral(N.connects[prev_C], N.index, CA.index, C.index, mov, type)
                push!(dihedrals, dihd)
            end
        end
        # find ψ (psi)
        if C.connects != nothing
            next_N = findfirst(i->metadata.atoms[i].name=="N", C.connects)
            if next_N != nothing
                type = DihedralTypes.psi
                mov = metadata.stree.indices[1:Common.getindex(metadata.stree, C.index)]
                dihd = Common.Dihedral(N.index, CA.index, C.index, C.connects[next_N], mov, type)
                push!(dihedrals, dihd)
            end
        end
        #endregion

        #region SIDECHAIN-DIHEDRALS
        if residue.name == "PRO"
            # no sidechains are considered for Prolines
            continue
        end

        # identify the number of possible chis in this residue
        n_chis = residue.name in ("HIS", "PHE", "TRP", "TYR") ? 3 : 6
        
        # concatenate all atom names in this residue
        # such that the regex can be applied to identify
        # the target atom name 
        atnames = string(":", join(map(atom -> atom.name, atoms), ":"), ":")

        path = ["N", "CA"]
        for χ = 1:n_chis
            m = match(queries[χ], atnames)
            if m == nothing
                continue
            end
            
            push!(path, m[:name])
            if length(path) > 3 
                a1, a2, a3, a4 = map(name -> name2index[name].index, path[end-3:end])
                type = DihedralTypes(χ-1)
                mov = metadata.stree.indices[1:Common.getindex(metadata.stree, a3)]
                dihd = Common.Dihedral(a1,a2,a3,a4,mov,type)
                push!(dihedrals, dihd)
                # deleteat!(path, 1)
            end
        end

        #endregion

    end
    dihedrals
end



end

# f(mutators::Vector{X}) where {X<:AbstractMutator}


# struct Sampler{F<:Function}
#     apply!::F
#     mutators::T...
# end

# Sampler(f::F, mutators::T...) where {F<:Function, T<:AbstractMutator} = begin
#     # function(state::Common.State)
#     #     f(state, mutators)
#     # end
#     Sampler(mutators, (state::Common.State)->f(state, mutators))
# end

# Sampler(mutators::T...) where {F<:Function, T<:AbstractMutator} = begin
#     Sampler(mutators, (state::Common.State)-> begin
#         for mutator in mutators
#             apply!(state, mutator)
#         end)
# end

# sampler = Sampler(f, mut1, mut2)
# sampler.apply!(state)

# my_apply = Mutators.@apply function (state, mut_confs)

#     if apply!(state, mutators_config[1]) == 0
#         apply!(state, mutators_config[2])
#     end
# end

# my_tuner = Mutators.@tuner function (mut_confs, driver_state)
#     if driver_state.step % 100 == 0
#         mutators_config[1].step_size = 1.0
#         mutators_config[2].step_size = 1.0
#     end
# end


# my_sampler = Sampler(my_apply, mut1_conf, mut2_conf, ...) = Sampler(my_apply, nothing, mut1_conf, mut2_conf, ...)
# my_sampler = Sampler(my_apply, my_tuner, mut1_conf, mut2_conf, ...)

# Sampler(function(state::Common.state, mutators_config::Vector{X}) where {X<:AbstractMutatorConfig}
#     if apply!(state, mutators_config[1]) == 0
#         apply!(state, mutators_config[2])
#     end
# end,
# mut1_conf, mut2_conf, [...])





# struct Evaluator where {T<AbstractForcefieldContainer}
#     evaluate!::Function
#     components::Union{T, Vector{T}}
# end
# Evaluator(f, args...)
# Evaluator(f, [arg1, arg2, ...])


# evaluator.evaluate!(state, evaluator.components)
# my_evaluator = Evaluator(Forcefield.@evaluate(), bonds, angles)
# my_evaluator = Evaluator(f, bonds)

# abstract type AbstractForcefieldComponent end
# const AbstractForcefieldContainer = Vector{AbstractForcefieldComponent}

# container = AbstractForcefieldContainer()

# bonds,
# angles
# top
# my_evaluator = Forcefield.@evaluator function(state, components::Union{T,Vector{T}}) where {T<:AbstractForcefieldContainer}

# # AS CHAMADAS

# sampler.apply!(state, sampler.mutators)
# if sampler.tune! != nothing
#     sampler.tune!(sampler.mutators, driver_state)
# end


# evaluator.evaluate!(state, evaluator.components, do_forces)
