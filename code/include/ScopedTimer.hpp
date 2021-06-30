
#pragma once

#include <iostream>
#include <chrono>

class ScopedTimer {

    std::chrono::steady_clock::time_point m_begin;
    std::string m_prefix;

public:
    ScopedTimer() : m_begin(std::chrono::steady_clock::now()) { }

    explicit ScopedTimer(const std::string& s) : ScopedTimer() {
        m_prefix = s;
    }

    long long elapsed() {
        const auto t = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(t - m_begin).count();
    }

    ~ScopedTimer() {
        if (!m_prefix.empty())
            std::cout << m_prefix << "\t... done in " << elapsed() << " ms" << std::endl;
    }
};
