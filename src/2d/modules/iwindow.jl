mutable struct MessageLog
    messages::Vector{String}
end

function add_message!(Log::MessageLog, text::String)
    push!(Log.messages, text)
end

function iwindow(log::MessageLog, fdsg2d::FDSG2D)
    if fdsg2d.settings.showinfo
        # Clear console
        #print("\e[H\e[2J")

        # Info-Window
        @printf("\n")
        @printf("╔══════════════════════════════════\n")
        @printf("║ FDSG 2D - Summary             \n")
        @printf("╠══════════════════════════════════\n")
        @printf("║ x-coords: %-1s: %-1s: %-1s        \n", fdsg2d.domain.x0,round(fdsg2d.domain.dx,digits=2),round(fdsg2d.domain.xend,digits=2))
        @printf("║ y-coords: %-1s: %-1s: %-1s        \n", fdsg2d.domain.y0,round(fdsg2d.domain.dy,digits=2),round(fdsg2d.domain.yend,digits=2))
        @printf("║ #DoFs: %-1d               \n", fdsg2d.domain.number_of_nodes)
        @printf("║ Appr. Memory: %-1.2f GB              \n", fdsg2d.domain.memory)
        @printf("║ time: %-1s: %-1s: %-1s             \n", fdsg2d.time.t0,round(fdsg2d.time.dt,digits=5),fdsg2d.time.tend)
        @printf("║ nt: %-1d               \n", fdsg2d.time.nt)
        @printf("║ device: %-1s               \n", fdsg2d.settings.device_name)
        @printf("║ source x: %-1.2f             \n", fdsg2d.source.x)
        @printf("║ source y: %-1.2f             \n", fdsg2d.source.y)
        @printf("╚══════════════════════════════════\n")

        # Message logs
        for message in log.messages
            println(message)
        end
    end
end