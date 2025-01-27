function check_device(device)

    device_list = ["cpu", "gpu-metal", "gpu-cuda"]


    # cpu
    if device == "cpu"
        return "cpu"
    
    # nvidia gpu
    elseif device == "gpu-cuda"

        device_ids = CUDA.devices()
        if isempty(device_ids)
            error("Cuda device not found.")
        else
            return CUDA.name(CUDA.device!(0))
        end

    # apple gpu
    elseif device == "gpu-metal"

        device_list = Metal.MTL.devices()

        if isempty(device_list)
            error("Metal device not found.")
        else 
            return device_list[1].name
        end

    # exit 
    else 
        error("Device: $device not found in $device_list")
    end

end;