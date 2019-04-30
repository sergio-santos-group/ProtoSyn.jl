# ProtoSyn Syntax Guidelines

## Variable definition
- Spaces should be used in attribution operations: `a = 1` instead of `a=1`;
- Spaces should be used after commas (in 2D indexing and multiple atrributions): `a[i, k]` instead of `a[i,k]`
- Multi-word **variables** should use the underscore notation: `cut_off = 1.0` instead of `cutOff = 1.0` (except for the name of energy components, who should use camelcase);
- Multi-word **constructors** should use camelcase notation: `MonteCarloConfig` instead of `monte_carlo_config`;
- Constructors should have a default constructor and an overloaded `Base.show` method for printing the struct;
- Dictionary entries should use camelcase notation: `energy.comp["eBonds"]` instead of `energy.comp["e_bonds"]`
- Variables defining 
    - **distances**  should start with a "**d**": `d12`;
    - **vectors**  should start with a "**v**": `v12`;
    - **constants**  should start with a "**c**": `c12`;
    - **squared constants**  should end with "**Sq**": `c12Sq`;
- New arrays should follow the `Array{Float64, 1}()` notation (and not `Float64[]`, for example);
- 1D arrays should be named "Vectors": `Vector{Atom}()`
- Variables initialized as zero should follow the `0.0` notation (instead of only `.0`, for example)
- Whenever feasible, indenting should follow the begin-end block structure, Python-like;
- Functions that alter *in-place* one of the inputs should follow Julia convention and end with a "`!`";

## Operations
- Squaring a number should be done using the `a ^ 2` notation (and not `a * a`, for example);
- Operations should be separated by a space: `a * b` instead of `a*b`;
- Functions should explictly define the return value, even though not necessary: `return a`
- Number incrementations should use the `a += 1` notation (and not `a = a + 1`, for example)

## Other
- Energy components use `e` as a prefix: `eAngles` instead of `Angles`