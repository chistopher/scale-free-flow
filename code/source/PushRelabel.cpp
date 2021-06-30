
#include <PushRelabel.h>

#include <cassert>
#include <iostream>

#include <utils.h>

#include <boost/graph/push_relabel_max_flow.hpp>


// some stuff needed for boost
struct IntIterator {
    size_t val = 0;
    IntIterator() = default;
    IntIterator(size_t start) : val(start) {}
    auto operator!=(const IntIterator& other) { return val != other.val; }
    auto operator==(const IntIterator& other) { return val == other.val; }
    auto operator*() { return val; }
    auto operator++() { ++val; }
    auto operator--() { --val; }
};

template<>
struct boost::graph_traits<GoldbergTarjan::LinearizedNetwork> {
    using vertex_descriptor = size_t;
    using edge_descriptor = size_t;
    using vertex_iterator = IntIterator;
    using out_edge_iterator = IntIterator;
    using vertices_size_type = size_t;
    using edges_size_type = size_t;
};

namespace GoldbergTarjan {
    struct EdgeCapMap { std::vector<Edge>& data; };
    struct EdgeResCapMap { std::vector<Edge>& data; };
    struct EdgeRevMap { std::vector<Edge>& data; };
    auto get(EdgeCapMap& m, size_t e) { return m.data[e].capacity; }
    //auto get(EdgeResCapMap& m, size_t e) { return m.data[e].residual_capacity < 1e-6 ? 0.0 : m.data[e].residual_capacity; } // very small == nothing
    double get(EdgeResCapMap& m, size_t e) { return m.data[e].residual_capacity; } // very small == nothing
    auto get(EdgeRevMap& m, size_t e) { return m.data[e].reverse_edge; }
    // auto put(EdgeCapMap& m, size_t e, double val) { m.data[e].capacity = val; }
    auto put(EdgeResCapMap& m, size_t e, double val) { m.data[e].residual_capacity = val; }
    // auto put(EdgeRevMap& m, size_t e, size_t val) { m.data[e].reverse_edge = val; }

    auto num_vertices(LinearizedNetwork& g) { return g.n; }
    auto num_edges(LinearizedNetwork& g) { return g.data.size(); }
    auto vertices(LinearizedNetwork& g) { return std::pair(IntIterator(0),IntIterator(g.n)); }
    auto out_edges(size_t v, LinearizedNetwork& g) { return std::pair(IntIterator(g.adj[v].edge_begin), IntIterator(g.adj[v].edge_end)); }
    auto source(size_t e, LinearizedNetwork& g) { return g.data[e].from; }
    auto target(size_t e, LinearizedNetwork& g) { return g.data[e].to; }
}

// pimpl pattern hidden data
namespace GoldbergTarjan {
    class InternalBoostAlg {
    public:
        boost::detail::push_relabel<LinearizedNetwork, EdgeCapMap, EdgeResCapMap, EdgeRevMap, boost::identity_property_map, double> m_alg;
        InternalBoostAlg(LinearizedNetwork& graph, size_t s, size_t t)
            : m_alg(graph, {graph.data}, {graph.data}, {graph.data}, static_cast<size_t>(s), static_cast<size_t>(t), boost::identity_property_map()) {
        }
    };
}
// end boost stuff

using namespace std;


GoldbergTarjan::LinearizedNetwork::LinearizedNetwork(const EdgeList &graph) : n(graph.n) {
    if(!graph.isDirected) printf("lin net supports only directed for now"), exit(0);

    adj.resize(n);
    data.resize(2*std::size(graph.edges));

    vector<size_t> deg(n,0); // first stores degree, then index for next edge
    for(auto [u,v] : graph.edges) deg[u]++, deg[v]++;
    adj[0].edge_begin = 0;
    for(int i=0; i<n; ++i) {
        adj[i].edge_end = adj[i].edge_begin + deg[i];
        deg[i] = adj[i].edge_begin;
        if(i+1<n) adj[i+1].edge_begin = adj[i].edge_end;
    }
    assert(adj[n-1].edge_end == 2*size(graph.edges));

    for(int i=0; i<size(graph.edges); ++i) {
        auto [u,v] = graph.edges[i];
        auto w = graph.hasWeights ? graph.weights[i] : 1.0;
        auto e1 = deg[u];
        data[e1] = Edge{w,w,(size_t)u,(size_t)v,0};
        deg[u]++;
        auto e2 = deg[v];
        data[e2] = Edge{0,0,(size_t)v,(size_t)u,0};
        deg[v]++;
        data[e1].reverse_edge = e2;
        data[e2].reverse_edge = e1;
    }

#ifndef NDEBUG
    assert(adj.front().edge_begin == size_t(0));
    assert(adj.back().edge_end == size(data));
    for(int i=0; i<n; ++i) {
        for(auto [be,en] = adj[i]; be != en; ++be) {
            assert(be<=en && en<=size(data));
            assert(data[be].from==i);
            assert(data[be].reverse_edge < size(data));
            assert(data[data[be].reverse_edge].from==data[be].to);
            assert(data[data[be].reverse_edge].to==i);
        }
    }
#endif
}

GoldbergTarjan::PushRelabel::PushRelabel(const EdgeList &graph)
: m_graph(graph.isDirected ? graph : undirectedToDirected(graph))
{ }

GoldbergTarjan::PushRelabel::~PushRelabel() {};

void GoldbergTarjan::PushRelabel::init(int s, int t) {
    if(m_state != EMPTY) cout << "WARNING: PR state was not empty before init\n";
    m_internal_alg = make_unique<InternalBoostAlg>(m_graph, static_cast<size_t>(s), static_cast<size_t>(t));
    m_state = INIT;
}

double GoldbergTarjan::PushRelabel::maxFlow(int s, int t) {

    if(m_state != INIT || !m_internal_alg || m_internal_alg->m_alg.src != s || m_internal_alg->m_alg.sink != t)
        cout << "ERROR: PR was not or incorrect initialized\n", exit(0);

    auto flow = m_internal_alg->m_alg.maximum_preflow();
    assert(m_internal_alg->m_alg.is_optimal());
    m_state = PREFLOW;

    return flow;
}

void GoldbergTarjan::PushRelabel::convertPreflow() {
    if(m_state != PREFLOW)
        cout << "ERROR: must compute preflow before flow\n", exit(0);
    m_internal_alg->m_alg.convert_preflow_to_flow();
    assert(m_internal_alg->m_alg.is_flow()); // may break when using non-integer weights
    m_state = FLOW;
}

void GoldbergTarjan::PushRelabel::deinit() {
    if(m_state==EMPTY) cout << "WARNING: PR deinit  called while empty\n";
    m_internal_alg.reset();
    m_state = EMPTY;
}


