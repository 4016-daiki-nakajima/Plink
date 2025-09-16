#pragma once

#include <iostream>
#include <vector>
#include "miniaudio.h"
#include "common.h"

namespace SimpleModal
{
    class AudioManager
    {
    public:
        AudioManager(ma_uint32 sampleRate = 44100, ma_uint32 bufferSize = 44100);
        ~AudioManager();

        // Initialize the audio system
        bool initialize();

        // Start/Stop playback
        bool startPlayback();
        bool stopPlayback();

        // Set audio data for playback
        void setAudioData(const std::vector<float> &data);

        // Reset playback position
        void resetPlayback();

        // Static callback for miniaudio
        static void data_callback(ma_device *device, void *output, const void *input, ma_uint32 frameCount);

    private:
        // Apply fade out to the end of the buffer to prevent clicks
        void applyFadeOut(std::vector<float> &buffer, int fadeLength);

        ma_context context;
        ma_device device;
        ma_device_config deviceConfig;

        std::vector<float> audioBuffer;
        ma_uint32 sampleRate;
        ma_uint32 bufferSize;
        ma_uint32 playbackPosition;

        bool isInitialized;
        bool isPlaying;

        // Number of samples for fade out (e.g., 1000 samples = ~23ms at 44.1kHz)
        static constexpr int FADE_OUT_SAMPLES = 1000;
    };
} // namespace SimpleModal