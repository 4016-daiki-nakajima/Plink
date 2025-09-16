#include "Audio/AudioManager.h"
#include "Utils/utils.h"
#include <cmath>

namespace SimpleModal
{

    AudioManager::AudioManager(ma_uint32 sampleRate, ma_uint32 bufferSize)
        : sampleRate(sampleRate), bufferSize(bufferSize), playbackPosition(0), isInitialized(false), isPlaying(false)
    {
        audioBuffer.resize(bufferSize, 0.0f);
    }

    AudioManager::~AudioManager()
    {
        if (isInitialized)
        {
            if (isPlaying)
            {
                stopPlayback();
            }
            ma_device_uninit(&device);
            ma_context_uninit(&context);
        }
    }

    bool AudioManager::initialize()
    {
        if (isInitialized)
        {
            return true;
        }

        // Initialize miniaudio context
        if (ma_context_init(NULL, 0, NULL, &context) != MA_SUCCESS)
        {
            std::cerr << "Failed to initialize audio context." << std::endl;
            return false;
        }

        // Configure device for playback
        deviceConfig = ma_device_config_init(ma_device_type_playback);
        deviceConfig.playback.format = ma_format_f32;
        deviceConfig.playback.channels = 2; // Stereo
        deviceConfig.sampleRate = sampleRate;
        deviceConfig.dataCallback = data_callback;
        deviceConfig.pUserData = this; // Pass this instance to callback

        // Initialize playback device
        if (ma_device_init(&context, &deviceConfig, &device) != MA_SUCCESS)
        {
            std::cerr << "Failed to initialize playback device." << std::endl;
            ma_context_uninit(&context);
            return false;
        }

        Utils::printSuccess("[AudioManager] Successfully initialized playback device: " + std::string(device.playback.name));

        isInitialized = true;
        return true;
    }

    bool AudioManager::startPlayback()
    {
        if (!isInitialized)
        {
            std::cerr << "Audio system not initialized." << std::endl;
            return false;
        }

        if (isPlaying)
        {
            return true;
        }

        if (ma_device_start(&device) != MA_SUCCESS)
        {
            std::cerr << "Failed to start playback device." << std::endl;
            return false;
        }

        isPlaying = true;
        return true;
    }

    bool AudioManager::stopPlayback()
    {
        if (!isPlaying)
        {
            return true;
        }

        if (ma_device_stop(&device) != MA_SUCCESS)
        {
            std::cerr << "Failed to stop playback device." << std::endl;
            return false;
        }

        isPlaying = false;
        return true;
    }

    void AudioManager::applyFadeOut(std::vector<float> &buffer, int fadeLength)
    {
        if (buffer.size() < fadeLength)
        {
            return; // Buffer too small for fade out
        }

        // Apply cosine fade out (smoother than linear fade)
        for (int i = 0; i < fadeLength; ++i)
        {
            float factor = 0.5f * (1.0f + cos(M_PI * i / fadeLength)); // 1 -> 0 smooth transition
            int idx = buffer.size() - fadeLength + i;
            buffer[idx] *= factor;
        }
    }

    void AudioManager::setAudioData(const std::vector<float> &data)
    {
        if (data.empty())
        {
            audioBuffer.clear();
            audioBuffer.resize(bufferSize, 0.0f);
            return;
        }

        // Find the last non-zero sample (or very small value)
        size_t lastNonZeroIdx = data.size() - 1;
        const float threshold = 1e-4f; // Threshold for "silence"
        while (lastNonZeroIdx > 0 && std::abs(data[lastNonZeroIdx]) < threshold)
        {
            --lastNonZeroIdx;
        }

        // Add some padding after the last non-zero sample for fade out
        size_t paddedSize = std::min(lastNonZeroIdx + FADE_OUT_SAMPLES + 1, data.size());

        if (paddedSize > bufferSize)
        {
            // Truncate if input data is larger than buffer
            audioBuffer.assign(data.begin(), data.begin() + bufferSize);
        }
        else
        {
            // Copy data and pad with zeros if smaller
            audioBuffer.assign(data.begin(), data.begin() + paddedSize);
            audioBuffer.resize(bufferSize, 0.0f);
        }

        // Apply fade out
        applyFadeOut(audioBuffer, FADE_OUT_SAMPLES);

        resetPlayback();
    }

    void AudioManager::resetPlayback()
    {
        playbackPosition = 0;
    }

    void AudioManager::data_callback(ma_device *device, void *output, const void *input, ma_uint32 frameCount)
    {
        AudioManager *audioManager = static_cast<AudioManager *>(device->pUserData);
        float *out = static_cast<float *>(output);
        ma_uint32 channels = device->playback.channels;

        for (ma_uint32 frame = 0; frame < frameCount; ++frame)
        {
            for (ma_uint32 channel = 0; channel < channels; ++channel)
            {
                if (audioManager->playbackPosition < audioManager->bufferSize)
                {
                    out[frame * channels + channel] = audioManager->audioBuffer[audioManager->playbackPosition];
                }
                else
                {
                    out[frame * channels + channel] = 0.0f; // Silence after buffer ends
                }
            }
            ++audioManager->playbackPosition;
        }
    }

} // namespace SimpleModal