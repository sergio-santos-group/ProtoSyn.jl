module Abstract

# Abstract types definition
abstract type DriverConfig end
abstract type DriverState  end
abstract type MutatorConfig end
abstract type Sampler end
abstract type Evaluator end
abstract type Test end
abstract type ForcefieldComponent <: Test end
const ForcefieldComponentContainer = Vector{T} where {T <: ForcefieldComponent}

end # module