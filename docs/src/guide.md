# Package guide

## Installation

Protosyn is still not a registered package. Please download or clone the [source](https://github.com/sergio-santos-group/ProtoSyn.jl/tree/6233638b9947e5e697c4f2b871b1ea301a6acde5) code from GitHub.

## Usage

What follows is a step-by-step guide on how to use the ProtoSyn library.

### General overview and workflow

ProtoSyn is a stuctural sampling library designed to explore the conformational space of molecules. Specifically, ProtoSyn was developed to study proteins and how they fold in 3D space. With ProtoSyn, the user is able to easily integrate different modules in order to perform distinct actions that change and evaluate the system state.

The **ProtoSyn Flow** is as follows:
1. Load the system initial [`Common.State`](@ref) from a structural file, such as a .pdb or a .gro;
2. Load all the necessary topologies, describing the system bonded and non-bonded interactions;
3. Choose a [Driver](@ref Drivers). This will change the system conformation step after step;
4. For certain Drivers, a sampler is required. This function is responsible for mutating the system in a certain way. Several [Mutators](@ref) are available.
5. Define the evaluating function. Regardless of the Driver chosen, this function will evaluate the system *fitness* and determine the next step.
6. Define the callback function. At this point the user's program requires a way to comunicate with the outside. Callback functions are employed to output the produced information by the Driver.
7. Run the program and do science!

### Examples

As an example, a Steepest Descent Algorithm written in Julia will be explained in detail. As explained in the **ProtoSyn Flow** we have the following steps:

- **Load the system initial [`Common.State`](@ref).**

```julia-repl
julia> state = Common.load_from_pdb("protein.pdb")
Common.State(size=56, energy=Null, xyz=[1.1 2.1 1.2; 2.2 1.2 5.2, ...], forces=[0.0 0.0 0.0; 0.0 0.0 0.0, ...], atnames=["C", "O", ...])
julia> state.energy = Forcefield.Energy()
```

- **Load all the necessary topologies.**

    The topology can be loaded by several ways. ProtoSyn by default expects a JSON file depicting the necessary information. For the Steepest Descent Algorithm, bonds, angles, dihedrals and non-bonded parameters should be depicted in the input JSON file.
```julia-repl
julia> topology = Forcefield.load_from_json("topology.json")
Forcefield.Topology(
atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name="N", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],
bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],
angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],
dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])
```

- **Choose a [Driver](@ref Drivers).**

    This step includes loading the necessary runtime parameters. These are specific to the chosen Driver.
```julia-repl
julia> params = Drivers.SteepestDescent.ConfigParameters(10000, 100, 1e-3, 0.1)
Drivers.SteepestDescent.ConfigParameters(n_steps=10000, log_freq=100, f_tol=1e-3, max_step=0.1)
```

- **Define the sampler.**

    As stated in the [documentation](@ref Drivers), the Steepest Descent Driver does not require a sampler, as the Driver itself is responsible for choosing the next step without the aid of a [Mutator](@ref Mutators).

- **Define the evaluating function.**

    When defining the necessary functions in ProtoSyn, careful care needs to be taken to match the expected function signature, described in the [documentation](@ref Drivers).
```julia-repl
julia> function my_evaluator!(st::Common.State, do_forces::Bool)
        return Forcefield.evalenergy!(topology, st, cut_off=1.2, do_forces=do_forces)
    end
my_evaluator! (generic function with 1 method)
```
!!! tip
    Even though the `my_evaluator!` function does not directly receive the topology as an argument, it is still able to access it as it was defined in the main body of our program.

- **Define the callback function.**

    Altough optional, this function allows the user to easily retrieve the information being produced by the Driver. ProtoSyn includes some functions for this, such as to [Print](@ref) the current state to a structural file.
```julia-repl
julia> output = open("output.xyz", "w")
IOStream(<output.xyz>)
julia> function my_callback(st::Common.State, step::Int)
            Print.as_xyz(st, ostream = output, title = "Step $step")
        end
my_callback (generic function with 1 method)
```

- **Run the program and do science!**

    All the necessary variables and functions have been defined in the main body of our program and it is ready to be deployed.
```julia-repl
julia> Drivers.SteepestDescent.run!(state, my_evaluator!, params, callback = my_callback)
```

For more detailed information, please reference to the Manual.