#ifndef GSNUNF_TIMER_H
#define GSNUNF_TIMER_H

#include <chrono>
#include <iostream>
#include <string>

namespace gsnunf {

class Timer
{
public:
    Timer( std::string delimiter = "," ) : m_delimiter(delimiter) {
        m_startTime = std::chrono::high_resolution_clock::now();
    }
    ~Timer() {
        stop();
    }
    void stop() {
        auto endTime = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTime).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
        auto duration = end - start;

        std::cout << duration << m_delimiter;
    }
private:
    std::chrono::time_point< std::chrono::high_resolution_clock > m_startTime;
    std::string m_delimiter;
};

} // namespace gsnunf

#endif // GSNUNF_TIMER_H

