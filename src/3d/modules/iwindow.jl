mutable struct MessageLog
    messages::Vector{String}
end

function add_message!(Log::MessageLog, text::String)
    push!(Log.messages, text)
end

function iwindow(log::MessageLog, fdsg3d::FDSG3D)
    if fdsg3d.settings.showinfo

        # Info-Window
        @printf("\n")
        @printf("╔══════════════════════════════════\n")
        @printf("║ FDSG 3D - Summary                \n")
        @printf("╠══════════════════════════════════\n")
        @printf("║ x-coords: %-1s: %-1s: %-1s       \n", fdsg3d.domain.x0,round(fdsg3d.domain.dx,digits=2),round(fdsg3d.domain.xend,digits=2))
        @printf("║ y-coords: %-1s: %-1s: %-1s       \n", fdsg3d.domain.y0,round(fdsg3d.domain.dy,digits=2),round(fdsg3d.domain.yend,digits=2))
        @printf("║ z-coords: %-1s: %-1s: %-1s       \n", fdsg3d.domain.z0,round(fdsg3d.domain.dz,digits=2),round(fdsg3d.domain.zend,digits=2))
        @printf("║ #DoFs: %-1d                      \n", fdsg3d.domain.number_of_nodes)
        @printf("║ Appr. Memory: %-1.2f GB          \n", fdsg3d.domain.memory)
        @printf("║ time: %-1s: %-1s: %-1s           \n", fdsg3d.time.t0,round(fdsg3d.time.dt,digits=5),fdsg3d.time.tend)
        @printf("║ nt: %-1d                         \n", fdsg3d.time.nt)
        @printf("║ device: %-1s                     \n", fdsg3d.settings.device_name)
        @printf("║ source x: %-1.2f                 \n", fdsg3d.source.x)
        @printf("║ source y: %-1.2f                 \n", fdsg3d.source.y)
        @printf("║ source z: %-1.2f                 \n", fdsg3d.source.z)
        @printf("╚══════════════════════════════════\n")

        # Message logs
        for message in log.messages
            println(message)
        end
    end
end