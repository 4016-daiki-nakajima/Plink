
#define MINIAUDIO_IMPLEMENTATION
#include "miniaudio.h"
#include <cmath>

// Constants for the sine wave
const double PI = 3.14159265358979323846;
const double FREQUENCY = 440.0; // A4 note
const double AMPLITUDE = 0.25;  // Volume
double phase = 0.0;

// Data callback function
void data_callback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount)
{
    float* out = (float*)pOutput;
    for (ma_uint32 i = 0; i < frameCount; ++i) {
        out[i] = static_cast<float>(AMPLITUDE * sin(phase));
        phase += 2.0 * PI * FREQUENCY / pDevice->sampleRate;
        if (phase >= 2.0 * PI) {
            phase -= 2.0 * PI;
        }
    }
}

int main()
{
    // Device configuration
    ma_device_config config = ma_device_config_init(ma_device_type_playback);
    config.playback.format   = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate        = 48000;
    config.dataCallback      = data_callback;

    // Initialize the device
    ma_device device;
    if (ma_device_init(NULL, &config, &device) != MA_SUCCESS) {
        return -1;
    }

    // Start playback
    if (ma_device_start(&device) != MA_SUCCESS) {
        ma_device_uninit(&device);
        return -1;
    }

    // Keep the application running to allow continuous playback
    printf("Press Enter to quit...\n");
    getchar();

    // Clean up
    ma_device_uninit(&device);
    return 0;
}
