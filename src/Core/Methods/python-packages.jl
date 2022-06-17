using PyCall

function print_installation_error(package_name::String, installation_instructions::String, current_module::String)
    println()
    @warn """
    📍 ProtoSyn was not able to identify `$package_name` in this system.
    PyCall is currently configured to use the Python version at $(PyCall.current_python()).
    In order to use the SeqDes energy function component, make sure:
    - To set ENV["PYTHON"] to the path of python executable you wish to use, run Pkg.build("PyCall") and re-launch Julia and ProtoSyn.
    - That `$package_name` is installed in the machine trying to load ProtoSyn.

    In order to install `$package_name`, follow the following instructions:
    $installation_instructions

    ProtoSyn will continue loading, but the `$current_module` module will be unavailable.
    To surpress further warnings for unavailable energy function components, set
    the JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC environment flag and re-launch
    Julia and ProtoSyn. 
    \$ export JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC=false
    Optionally, add the above line to ~/.bashrc to persistently supress warnings
    in further sessions.

    """
end # function

function is_python_package_installed(package_name::String, installation_instructions::String, current_module::String)

    try
        pyimport(package_name)
    catch LoadError
        if !("JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC" in keys(ENV))
            ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] = true
        end
        if ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] === "true"
            print_installation_error(
                package_name,
                installation_instructions, 
                current_module)
        end # if

        return false
    end #try

    return true
end # function