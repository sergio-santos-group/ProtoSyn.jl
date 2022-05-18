using ProtoSyn

# Using the ProtoSyn.eval(:(...)), PolarSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    export ProteinSelection

    """
    # TODO: Documentation
    """
    mutable struct ProteinSelection{M, T} <: AbstractSelection
        ProteinSelection() = begin
            new{ProtoSyn.Stateless, Residue}()
        end
    end

    # --- Select -------------------------------------------------------------------
    function select(::ProteinSelection, container::ProtoSyn.AbstractContainer)

        n_residues = count_residues(container)
        mask = Mask{Residue}(n_residues)

        known_aminoacids = [x[1] for x in collect(keys(Peptides.grammar.variables)) if x !== "B"]
        for residue in eachresidue(container)
            if (residue.name in keys(ProtoSyn.three_2_one)) && (ProtoSyn.three_2_one[residue.name] in known_aminoacids)
                mask[residue.index] = true
            end
        end
        return mask
    end

    state_mode_type(::ProteinSelection{M, T}) where {M, T} = M
    selection_type(::ProteinSelection{M, T}) where {M, T} = T

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, ps::ProteinSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        println(io, lead*"ProteinSelection › ($(selection_type(ps)))")
    end
end))