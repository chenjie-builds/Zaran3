#pragma once
#include <random>
namespace zaran
{
inline double RandDouble(double min, double max)
{
    static std::random_device rd;                     // Obtain a random number from hardware
    static std::mt19937 eng(rd());                    // Seed the generator
    std::uniform_real_distribution<> distr(min, max); // Define the range
    return distr(eng);
}
inline int RandInt(int min, int max)
{
    static std::random_device rd;                    // Obtain a random number from hardware
    static std::mt19937 eng(rd());                   // Seed the generator
    std::uniform_int_distribution<> distr(min, max); // Define the range
    return distr(eng);
}
} // namespace zaran