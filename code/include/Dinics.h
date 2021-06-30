
#pragma once

#include <vector>
#include <array>

#include <utils.h>
#include <DinicsStats.h>

namespace Dinics {

struct Edge {
    int to;
    double flow, cap;
    Edge* rev;
    bool isSaturated() const { return cap-flow < 1e-6; }
};

struct Neighbors {
    Edge* m_begin;
    Edge* m_end;
    auto begin() { return m_begin; }
    auto end() { return m_end; }
    const auto begin() const { return m_begin; }
    const auto end() const { return m_end; }
    Edge& operator[](size_t idx) { return m_begin[idx]; }
    const Edge& operator[](size_t idx) const { return m_begin[idx]; }
    size_t size() const { return m_end - m_begin; }
};
inline size_t size(const Neighbors& nei) { return nei.m_end - nei.m_begin; }

// a node type with conditional member mixins
template<bool COND> struct HasDistT {};
template<bool COND> struct HasStamp {};
template<> struct HasDistT<true> { int distT = 0; };
template<> struct HasStamp<true> { long long stamp = -1; };
template<bool BI, bool STAMPS>
struct Node : public HasDistT<BI>, public HasStamp<STAMPS> {
    Neighbors edges; // outgoing edges
    int distS = 0;
    int nextEdge = 0;
};


template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
class Dinics {
public:
    using node_type = Node<BI,STAMPS>;
    constexpr static bool isBI = BI;
    constexpr static bool isSTAMPS = STAMPS;
    constexpr static bool isRESET = RESET;
    constexpr static bool isSKIP = SKIP;
    constexpr static bool isADJREV = ADJREV;

    explicit Dinics(const EdgeList& graph);

    void resetFlow();
    double maxFlow(int s, int t, bool collectStats = false);

    const auto& graph() const { return m_graph; };
    const auto& stats() const { return m_stats; };

protected:

    // internal methods
    void clearNodeData();
    double dfs(int v, int t, double aug);
    bool bfs(int s, int t);
    bool bfsImplUni(int s, int t);
    bool bfsImplBi(int s, int t);

    // gathers stats based on last BFS
    void gatherStatsBFS(long long previousFlow, bool reachable);
    void gatherStatsDFS();

    // graph data
    std::vector<node_type> m_graph;
    std::vector<Edge> m_edges;
    size_t m_n;

    // bfs data
    std::vector<int> m_qs; // queues
    std::vector<int> m_qt; // queues
    int m_maxLayerS;
    int m_maxLayerT;

    std::vector<double*> m_changed; // all changed flow values to reset for next flow
    long long m_time; // current timestamp

    // statistics about last maxFlow computation
    std::vector<RoundStatistics> m_stats;
};

// representatives of different modifications ranked by gain
using Dinics0Vanilla= Dinics<0,0,0,0,0>;
using Dinics1Bi     = Dinics<1,0,0,0,0>;
using Dinics2Reset  = Dinics<1,0,1,0,0>;
using Dinics3Stamps = Dinics<1,1,1,0,0>;
using Dinics4Skip   = Dinics<1,1,1,1,0>;
using Dinics5OPT    = Dinics<1,1,1,1,1>;

// prevent compiler to generate templates
extern template class Dinics<0,0,0,0,0>;
extern template class Dinics<1,0,0,0,0>;
extern template class Dinics<0,1,0,0,0>;
extern template class Dinics<1,1,0,0,0>;
extern template class Dinics<0,0,1,0,0>;
extern template class Dinics<1,0,1,0,0>;
extern template class Dinics<0,1,1,0,0>;
extern template class Dinics<1,1,1,0,0>;
extern template class Dinics<0,0,0,1,0>;
extern template class Dinics<1,0,0,1,0>;
extern template class Dinics<0,1,0,1,0>;
extern template class Dinics<1,1,0,1,0>;
extern template class Dinics<0,0,1,1,0>;
extern template class Dinics<1,0,1,1,0>;
extern template class Dinics<0,1,1,1,0>;
extern template class Dinics<1,1,1,1,0>;
extern template class Dinics<0,0,0,0,1>;
extern template class Dinics<1,0,0,0,1>;
extern template class Dinics<0,1,0,0,1>;
extern template class Dinics<1,1,0,0,1>;
extern template class Dinics<0,0,1,0,1>;
extern template class Dinics<1,0,1,0,1>;
extern template class Dinics<0,1,1,0,1>;
extern template class Dinics<1,1,1,0,1>;
extern template class Dinics<0,0,0,1,1>;
extern template class Dinics<1,0,0,1,1>;
extern template class Dinics<0,1,0,1,1>;
extern template class Dinics<1,1,0,1,1>;
extern template class Dinics<0,0,1,1,1>;
extern template class Dinics<1,0,1,1,1>;
extern template class Dinics<0,1,1,1,1>;
extern template class Dinics<1,1,1,1,1>;

} // end namespace Dinics
