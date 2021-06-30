
#pragma once

#include <iostream>

class ProgressIndicator {
    size_t m_n;
    bool m_active;

public:
    explicit ProgressIndicator(size_t n, bool active = true)
    : m_n(n)
    , m_active(active)
    {
        if(m_active)
            std::cout << "progress 0%\r" << std::flush;
    }

    void tick(int i) {
        if(!m_active) return;
        if (100ll * (i + 1) / m_n != 100ll * i / m_n)
            std::cout << "progress " << 100ll * (i + 1) / m_n << "%\r" << std::flush;
    }

};
