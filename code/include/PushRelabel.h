
#pragma once

#include <vector>
#include <utility>
#include <memory>

#include <utils.h>

namespace GoldbergTarjan {

    class InternalBoostAlg; // FWD

    struct Edge {
        double capacity;
        double residual_capacity;
        size_t from;
        size_t to;
        size_t reverse_edge;
    };

    struct Node {
        size_t edge_begin;
        size_t edge_end;
    };

    struct LinearizedNetwork {
        size_t n;
        std::vector<Node> adj; // adjacency list
        std::vector<Edge> data; // consecutive edge data sorted by source end
        explicit LinearizedNetwork(const EdgeList& graph);
    };


    class PushRelabel {
        LinearizedNetwork m_graph;
        std::unique_ptr<InternalBoostAlg> m_internal_alg;
        enum { EMPTY, INIT, PREFLOW, FLOW } m_state = EMPTY;
    public:
        explicit PushRelabel(const EdgeList& graph);
        ~PushRelabel();
        void init(int s, int t); // must be called once before each flow
        double maxFlow(int s, int t);
        void convertPreflow(); // can be called once after each flow
        void deinit(); // can be called anytime to make the next init faster
        const auto& graph() const { return m_graph; };
    };

} // end namespace PR

