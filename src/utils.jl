# logger
const LOG_LEVELS = Dict(
    :debug => 1,
    :info  => 2,
    :warn  => 3,
    :error => 4,
)
const LOG_LEVEL = Ref(:info)

function set_log_level!(level::Symbol)
    haskey(LOG_LEVELS, level) || throw(ArgumentError("Invalid log level: $level"))
    LOG_LEVEL[] = level
end

macro logger(level, msg)
    quote
        local level = $(esc(level))
        if LOG_LEVELS[level] >= LOG_LEVELS[LOG_LEVEL[]]
            local timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
            local file = basename($(String(__source__.file)))
            local line = $(__source__.line)
            if LOG_LEVELS[level] >= 3
                println("[$(uppercase(String(level)))] [$timestamp] $($(esc(msg)))")
            else 
                println("[$(uppercase(String(level)))] [$timestamp] $($(esc(msg))) ($file:$line)")
            end
        end
    end
end