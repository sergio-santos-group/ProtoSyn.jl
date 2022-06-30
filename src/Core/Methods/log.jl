using Logging

function set_logger_to_debug()
    logger = ConsoleLogger(stdout, Logging.Debug)
    global_logger(logger)
    global_logger(logger) # Must be done twice...
    printstyled("ProtoSyn will now display debug snippets, informations, warnings and error messages. ", color = :green)
    println("Debug: ✓ | Info: ✓ | Warnings: ✓ | Errors: ✓")
end

function set_logger_to_info()
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)
    global_logger(logger) # Must be done twice...
    printstyled("ProtoSyn will now display informations, warnings and error messages. ", color = :blue)
    println("Debug: ⨯ | Info: ✓ | Warnings: ✓ | Errors: ✓")
end

function set_logger_to_warn()
    logger = ConsoleLogger(stdout, Logging.Warn)
    global_logger(logger)
    global_logger(logger) # Must be done twice...
    printstyled("ProtoSyn will now display warnings and error messages. ", color = :yellow)
    println("Debug: ⨯ | Info: ⨯ | Warnings: ✓ | Errors: ✓")
end

function set_logger_to_error()
    logger = ConsoleLogger(stdout, Logging.Error)
    global_logger(logger)
    global_logger(logger) # Must be done twice...
    printstyled("ProtoSyn will now only display error messages. ", color = :red)
    println("Debug: ⨯ | Info: ⨯ | Warnings: ⨯ | Errors: ✓")
end

function print_loading(msg::String; color::Symbol = :blue)
    printstyled("[ Loading: ", color = :blue, bold = true)
    printstyled(msg*"\n", color = color)
end

set_logger_to_error()