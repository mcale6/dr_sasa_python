#pragma once

#include <string>

class GPUConfig {
public:
    enum class ComputeMode {
        CPU,
        OpenCL,
        CUDA
    };

    struct DeviceInfo {
        std::string name;
        std::string vendor;
        bool available;
    };

    static ComputeMode get_preferred_mode() {
        #ifdef __APPLE__
            return ComputeMode::CPU;  // Always use CPU on Apple
        #else
            // For other platforms, could add OpenCL/CUDA detection here
            return ComputeMode::CPU;
        #endif
    }

    static DeviceInfo get_device_info(ComputeMode mode) {
        DeviceInfo info;
        info.name = "CPU";
        info.vendor = "Generic";
        info.available = true;
        return info;
    }
};