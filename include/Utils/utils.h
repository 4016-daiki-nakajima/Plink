#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <vector>

extern const double speed_of_sound;
extern const double air_density;

#define TIC(name) auto start_##name = std::chrono::high_resolution_clock::now(); 
#define TOC(name) \
    auto end_##name = std::chrono::high_resolution_clock::now(); \
    std::cout << Utils::CYAN << "[Time] " << #name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end_##name- start_##name).count() << "ms" << Utils::RESET << std::endl; \

namespace Utils
{

    // ANSI escape codes for colors
    const std::string RED = "\033[31m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string BLUE = "\033[34m";
    const std::string MAGENTA = "\033[35m";
    const std::string CYAN = "\033[36m";
    const std::string RESET = "\033[0m";

    // Print colored text
    inline void printRed(const std::string &text)
    {
        std::cout << RED << text << RESET << std::endl;
    }

    inline void printError(const std::string &text)
    {
        printRed(text);
    }
    // ------------------------------------------------------------

    inline void printGreen(const std::string &text)
    {
        std::cout << GREEN << text << RESET << std::endl;
    }
    inline void printSuccess(const std::string &text)
    {
        printGreen(text);
    }
    // ------------------------------------------------------------

    inline void printYellow(const std::string &text)
    {
        std::cout << YELLOW << text << RESET << std::endl;
    }

    inline void printWarning(const std::string &text)
    {
        printYellow(text);
    }
    // ------------------------------------------------------------

    inline void printBlue(const std::string &text)
    {
        std::cout << BLUE << text << RESET << std::endl;
    }

    inline void printMagenta(const std::string &text)
    {
        std::cout << MAGENTA << text << RESET << std::endl;
    }

    // ------------------------------------------------------------
    inline void printCyan(const std::string &text)
    {
        std::cout << CYAN << text << RESET << std::endl;
    }
    inline void printInfo(const std::string &text)
    {
        printCyan(text);
    }
    // ------------------------------------------------------------

    // Print colored text without newline
    inline void printColoredText(const std::string &text, const std::string &color)
    {
        std::cout << color << text << RESET;
    }

    std::vector<float> setSpectrumColor(float val, float minVal=-1.0f, float maxVal=1.0f);

} // namespace Utils

#endif // UTILS_H
