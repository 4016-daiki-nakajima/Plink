#include "utils.h"

const double speed_of_sound = 343.0;
const double air_density = 1.2041;

std::vector<float> Utils::setSpectrumColor(float val, float minVal, float maxVal)
{
    if (val > maxVal)
        val = maxVal;
    if (val < minVal)
        val = minVal;

    float t = 2 * ((val - minVal) / (maxVal - minVal)) - 1;

    std::vector<float> blue = {0.0f, 0.2f, 1.0f};
    std::vector<float> black = {0.0f, 0.0f, 0.0f};
    std::vector<float> orange = {1.0f, 0.1f, 0.0f};

    float r, g, b;
    if (t < 0.0f)
    {
        float s = -t;
        r = black[0] + (blue[0] - black[0]) * s;
        g = black[1] + (blue[1] - black[1]) * s;
        b = black[2] + (blue[2] - black[2]) * s;
    }
    else
    {
        float s = t;
        r = black[0] + (orange[0] - black[0]) * s;
        g = black[1] + (orange[1] - black[1]) * s;
        b = black[2] + (orange[2] - black[2]) * s;
    }

    return {r, g, b};
}
