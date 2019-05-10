module Abstract

# Abstract types definition
abstract type DriverConfig end
struct NullDriverConfig <: DriverConfig end
abstract type DriverState  end
abstract type MutatorConfig end
abstract type Sampler end
abstract type Evaluator end
abstract type ForcefieldComponent end
const ForcefieldComponentContainer = Vector{T} where {T <: ForcefieldComponent}

end # module