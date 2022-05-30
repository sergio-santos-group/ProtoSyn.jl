using Logging

function set_logger_to_debug()
    logger = ConsoleLogger(stderr, Logging.Debug)
    global_logger(logger)
    printstyled("ProtoSyn will now display debug snippets, informations, warnings and error messages.\n", color = :green)
    println("Debug: ✓ | Info: ✓ | Warnings: ✓ | Errors: ✓")
end

function set_logger_to_info()
    logger = ConsoleLogger(stderr, Logging.Info)
    global_logger(logger)
    printstyled("ProtoSyn will now display informations, warnings and error messages.\n", color = :blue)
    println("Debug: ⨯ | Info: ✓ | Warnings: ✓ | Errors: ✓")
end

function set_logger_to_warn()
    logger = ConsoleLogger(stderr, Logging.Warn)
    global_logger(logger)
    printstyled("ProtoSyn will now display warnings and error messages.\n", color = :yellow)
    println("Debug: ⨯ | Info: ⨯ | Warnings: ✓ | Errors: ✓")
end

function set_logger_to_error()
    logger = ConsoleLogger(stderr, Logging.Error)
    global_logger(logger)
    printstyled("ProtoSyn will now only display error messages.\n", color = :red)
    println("Debug: ⨯ | Info: ⨯ | Warnings: ⨯ | Errors: ✓")
end

set_logger_to_info()