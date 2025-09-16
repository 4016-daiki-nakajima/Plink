#include <iostream>
#include <random>
#include "miniaudio.h"

// Callback function to generate white noise
void data_callback(ma_device* device, void* output, const void* input, ma_uint32 frameCount)
{
    float* out = static_cast<float*>(output);
    static std::default_random_engine generator;
    static std::uniform_real_distribution<float> distribution(-1.0f, 1.0f);

    for (ma_uint32 i = 0; i < frameCount * device->playback.channels; ++i) {
        out[i] = distribution(generator);
    }
}

int main()
{
    // Initialize miniaudio context
    ma_context context;
    if (ma_context_init(NULL, 0, NULL, &context) != MA_SUCCESS) {
        std::cerr << "Failed to initialize context." << std::endl;
        return -1;
    }

    // Configure device for playback
    ma_device_config deviceConfig = ma_device_config_init(ma_device_type_playback);
    deviceConfig.playback.format = ma_format_f32;
    deviceConfig.playback.channels = 2; // Stereo
    deviceConfig.sampleRate = 44100;
    deviceConfig.dataCallback = data_callback;

    // Initialize playback device
    ma_device device;
    if (ma_device_init(&context, &deviceConfig, &device) != MA_SUCCESS) {
        std::cerr << "Failed to initialize playback device." << std::endl;
        ma_context_uninit(&context);
        return -1;
    }

    // Start playback
    if (ma_device_start(&device) != MA_SUCCESS) {
        std::cerr << "Failed to start playback device." << std::endl;
        ma_device_uninit(&device);
        ma_context_uninit(&context);
        return -1;
    }

    std::cout << "Playing white noise... Press Enter to stop." << std::endl;
    std::cin.get();

    // Clean up
    ma_device_uninit(&device);
    ma_context_uninit(&context);

    return 0;
}
