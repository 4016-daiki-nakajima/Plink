#include <iostream>
#include <vector>
#include "miniaudio.h"
#include <cmath>
// Sample rate and buffer size
constexpr ma_uint32 sampleRate = 44100;
constexpr ma_uint32 bufferSize = sampleRate; // 1 second of audio

// Global variables
std::vector<float> myAudioData(bufferSize); // Your precomputed audio data
ma_uint32 playbackPosition = 0;

// Callback function to stream audio data
void data_callback(ma_device* device, void* output, const void* input, ma_uint32 frameCount)
{
    float* out = static_cast<float*>(output);
    ma_uint32 channels = device->playback.channels;

    for (ma_uint32 frame = 0; frame < frameCount; ++frame) {
        for (ma_uint32 channel = 0; channel < channels; ++channel) {
            if (playbackPosition < bufferSize) {
                out[frame * channels + channel] = myAudioData[playbackPosition];
            } else {
                out[frame * channels + channel] = 0.0f; // Silence after buffer ends
            }
        }
        ++playbackPosition;
    }
}

int main()
{
    // Initialize your audio data (example: fill with a sine wave)
    float frequency = 440.0f; // A4 note
    for (ma_uint32 i = 0; i < bufferSize; ++i) {
        myAudioData[i] = 0.5f * float( bufferSize - i ) / bufferSize * sinf((2.0f * M_PI * frequency * i) / sampleRate);
    }

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
    deviceConfig.sampleRate = sampleRate;
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

    std::cout << "Playing audio buffer... Press Enter to stop." << std::endl;
    std::cin.get();

    // Clean up
    ma_device_uninit(&device);
    ma_context_uninit(&context);

    return 0;
}
