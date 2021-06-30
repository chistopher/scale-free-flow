
#pragma once

#include <array>
#include <tuple>

#include <utils.h>

namespace {

    auto testGraph = readGraph("testdata.txt");
    static constexpr std::array pairs{ // some pairs with flow values
            std::tuple(1, 2, 3),
            std::tuple(13, 378, 5),
            std::tuple(100, 200, 4),
            std::tuple(200, 100, 4),
    };

}
