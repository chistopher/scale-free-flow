
#pragma once

#include <vector>
#include <array>
#include <utility>

namespace MinimalDinics {

    struct edge {
        int from, to;
        int flow, cap;
        std::size_t rev;
    };

    long long dinics(int s, int t, std::vector<std::vector<edge>>& adj);

    std::vector<std::vector<edge>> buildNetwork(const std::vector<std::pair<int,int>>& edges, int n);
    std::vector<std::vector<edge>> buildNetwork(const std::vector<std::array<int,3>>& edges, int n);

    std::vector<int> cutSide(int s, int t, const std::vector<std::vector<edge>>& adj, bool sideS);
}
