
#include <Dinics.h>

#include <cassert>
#include <algorithm>

const int INF = 1e9;

namespace Dinics {

// helper that uses given space to simulate a queue
// only works if less pushes than size of vec
template<typename T>
struct QueueWrapper {
    explicit QueueWrapper(std::vector<T>& mem) : vec(mem) {};
    std::vector<T>& vec;
    size_t beg = 0;
    size_t end = 0;
    void push(T node) {
        vec[end++] = node;
    }
    T pop() {
        return vec[beg++];
    }
    auto size() {
        return end-beg;
    }
};

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::Dinics(const EdgeList &graph)
: m_graph(graph.n)
, m_edges(2*size(graph.edges))
, m_n(graph.n)
, m_qs(graph.n)
, m_qt(graph.n)
, m_maxLayerS(0)
, m_maxLayerT(0)
, m_changed()
, m_time(0)
, m_stats()
{
    if(graph.isDirected && (isADJREV && isBI)) std::cout << "cannot apply ADJREV improvement to directed graph", exit(0);

    std::vector<size_t> deg(m_n,0); // first stores degree, then index for next edge
    for(auto [u,v] : graph.edges) deg[u]++, deg[v]++;
    m_graph[0].edges.m_begin = m_edges.data();
    for(int i=0; i<m_n; ++i) {
        m_graph[i].edges.m_end = m_graph[i].edges.m_begin + deg[i];
        deg[i] = m_graph[i].edges.m_begin - m_edges.data();
        if(i+1<m_n) m_graph[i+1].edges.m_begin = m_graph[i].edges.m_end;
    }
    assert(m_graph[m_n-1].edges.m_end - m_edges.data() == 2*size(graph.edges));

    for(int i=0; i<size(graph.edges); ++i) {
        auto [a,b] = graph.edges[i];
        auto w1 = graph.hasWeights ? graph.weights[i] : 1.0;
        auto w2 = graph.isDirected ? 0.0 : w1;
        auto e1 = &(m_edges[deg[a]++] = Edge{b,0,w1,nullptr});
        auto e2 = &(m_edges[deg[b]++] = Edge{a,0,w2,nullptr});
        e1->rev = e2;
        e2->rev = e1;
    }

#ifndef NDEBUG
    assert(m_graph.front().edges.m_begin == m_edges.data());
    assert(m_graph.back().edges.m_end == m_edges.data() + m_edges.size());
    for(int i=0; i<m_n; ++i) {
        for(auto [be,en] = m_graph[i].edges; be != en; ++be) {
            assert(m_edges.data() <= be && be<=en && en<=m_edges.data()+m_edges.size());
            assert(m_edges.data() <= be->rev && be->rev < m_edges.data() + m_edges.size());
            assert(be->rev->to==i);
        }
    }
#endif
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
void Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::resetFlow() {
    m_maxLayerS = 0;
    m_maxLayerT = 0;
    if constexpr (RESET) {
        for(auto* f_ptr : m_changed) *f_ptr = 0;
        m_changed.clear();
    } else {
        for(auto& e : m_edges)
            e.flow = 0;
    }
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
double Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::maxFlow(int s, int t, bool collectStats) {
    if(m_maxLayerS || m_maxLayerT) std::cout << "ERROR: Dincis was not reset between flows\n", exit(0);
    m_stats.clear();
    auto flow = 0.0;
    while(true) {
        clearNodeData();
        auto found = bfs(s, t);
        if(collectStats) gatherStatsBFS(static_cast<long long>(flow),found);
        if(!found) break;
        double aug;
        while((aug = dfs(s, t, INF))>0) flow += aug;
        if(collectStats) gatherStatsDFS();
    }
    return flow;
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
void Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::clearNodeData() {
    if constexpr (!STAMPS) {
        for(node_type & node : m_graph) {
            node.distS = INF;
            node.nextEdge = 0;
            if constexpr (BI) node.distT = INF;
        }
    } else {
        m_time++;
    }
}


template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
double Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::dfs(int v, int t, double aug) {
    if (aug < 1e-6) return 0;
    if (v == t)
        return aug;

    auto& node = m_graph[v];
    if constexpr (STAMPS) assert(node.stamp==m_time);

    const auto deg = std::size(node.edges);
    for (int& i = node.nextEdge; i<deg; ++i) {
        auto& e = node.edges[i];
        auto& nei = m_graph[e.to];

        // #### long block of code to determine if we can use this edge ####
        if (e.isSaturated()) continue;
        if constexpr (STAMPS) if(nei.stamp != m_time) continue; // old timestamp is not relevant for sure
        if constexpr (BI) { // BIDIRECTIONAL CASE
            if constexpr (SKIP) {
                if(node.distS+1 > m_maxLayerS && node.distT-1 != nei.distT) continue; // if beyond S space, only visit with decreasing t-dist
                if(node.distS+1 ==m_maxLayerS && nei.distT == INF) continue; // do not visit into last unless in intersection
                if(node.distS+1 < m_maxLayerS && node.distS+1 != nei.distS) continue; // if in S space, only visit with increasing s-dist
            } else {
                if(node.distS+1 != nei.distS && node.distT-1 != nei.distT) continue; // proceed if either S dist increases or T dist decreases
            }
        } else { // UNIDIRECTIONAL CASE
            if constexpr (SKIP) if(nei.distS == m_maxLayerS && e.to != t) continue; // skip nodes in last layer that are not t
            if(node.distS+1 != nei.distS) continue; // only visit with increasing s-dist
        }
        // #### end of long block ####

        auto pushed = dfs(e.to, t, std::min(aug, e.cap - e.flow));
        if (pushed <= 0) continue;

        if constexpr (RESET) {
            m_changed.push_back(&e.flow);
            m_changed.push_back(&e.rev->flow);
        }

        e.flow += pushed;
        e.rev->flow -= pushed;
        return pushed;
    }
    return 0;
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
bool Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::bfs(int s, int t) {
    if constexpr (!BI) {
        return bfsImplUni(s,t);
    } else {
        return bfsImplBi(s,t);
    }
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
bool Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::bfsImplUni(int s, int t) {
    m_maxLayerS = 0;
    m_graph[s].distS = 0;
    if constexpr (STAMPS) {
        m_graph[s].stamp = m_time;
        m_graph[s].nextEdge = 0;
    }
    QueueWrapper<int> q(m_qs); q.push(s);
    bool found = false;
    while(q.size()) {
        m_maxLayerS++;
        auto layerSize = q.size();
        while (layerSize--) {
            auto& node = m_graph[q.pop()];
            for(auto& e : node.edges) {
                node_type& nei = m_graph[e.to];
                if(e.isSaturated()) continue;
                if constexpr (STAMPS) {
                    if(nei.stamp == m_time) continue;
                    nei.stamp = m_time;
                    nei.nextEdge = 0;
                } else {
                    if(nei.distS != INF) continue;
                }
                nei.distS = node.distS+1;
                q.push(e.to);
                found = found || e.to==t;
                if constexpr (ADJREV) if(found) return true;
            }
        }
        if constexpr (!ADJREV) if(found) return true; // stop bfs only when done with t's layer
    }
    return false;
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
bool Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::bfsImplBi(int s, int t) {
    if constexpr (BI) { // trick to prevent compile errors when this is instantiated for !BI

    m_maxLayerS = 0;
    m_maxLayerT = 0;

    // sets distances and updates stamp
    auto visitNode = [this](node_type &node, int distS, int distT) {
        node.distS = distS;
        node.distT = distT;
        if constexpr (STAMPS) {
            node.stamp = m_time;
            node.nextEdge = 0;
        }
    };

    visitNode(m_graph[s], 0, INF);
    visitNode(m_graph[t], INF, 0);

    // prepare queues
    QueueWrapper qs(m_qs);
    qs.push(s);
    QueueWrapper qt(m_qt);
    qt.push(t);

    auto qs_edges = std::size(m_graph[s].edges); // TODO measure this or make this optional
    auto qt_edges = std::size(m_graph[t].edges);

    bool found = false;
    while (!found && qs.size() && qt.size()) {

        if (qs_edges < qt_edges) {
            // do one layer forward
            m_maxLayerS++;
            qs_edges = 0;
            auto layerSize = qs.size();
            while (layerSize--) {
                for (const auto &e : m_graph[qs.pop()].edges) {
                    node_type &nei = m_graph[e.to];
                    if (e.isSaturated()) continue;
                    bool unseen;
                    if constexpr (STAMPS) unseen = nei.stamp != m_time;
                    else unseen = nei.distS == INF && nei.distT == INF;
                    if (unseen) { // unseen node
                        visitNode(nei, m_maxLayerS, INF);
                        qs.push(e.to);
                        qs_edges += size(nei.edges);
                    } else if (nei.distT != INF) { // node already seen by backwards search
                        found = true;
                        nei.distS = m_maxLayerS;
                    }
                }
            }
        } else {
            // do one layer backward
            m_maxLayerT++;
            qt_edges = 0;
            auto layerSize = qt.size();
            while (layerSize--) {
                for (const auto &e : m_graph[qt.pop()].edges) {
                    node_type &nei = m_graph[e.to];
                    if constexpr (ADJREV) {
                        if (e.cap - -e.flow < 1e-6) continue; // cool trick (assumes cap of twin edges is same)
                    } else {
                        if (e.rev->isSaturated())
                            continue; // look at residual capacity on twin edge
                    }
                    bool unseen;
                    if constexpr (STAMPS) unseen = nei.stamp != m_time;
                    else unseen = nei.distS == INF && nei.distT == INF;
                    if (unseen) { // new node
                        visitNode(nei, INF, m_maxLayerT);
                        qt.push(e.to);
                        qt_edges += size(nei.edges);
                    } else if (nei.distS != INF) { // node already seen by forward search
                        found = true;
                        nei.distT = m_maxLayerT;
                    }
                }
            }
        }

    } // end search
    return found;

    } else return s+t; // end if constexpr (BI) // used s,t to disable unused warning
}


template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
void Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::gatherStatsBFS(long long previousFlow, bool reachable) {

    // accessors that handle stamps
    auto distS = [this](const node_type& node) {
        if constexpr (STAMPS) if(node.stamp!=m_time) return INF;
        return node.distS;
    };
    auto distT = [this](const node_type& node) {
        if constexpr (STAMPS) if(node.stamp!=m_time) return INF;
        if constexpr (BI) return node.distT;
        return INF;
    };

    RoundStatistics round{};
    round.flowBefore = previousFlow;
    round.distOfSink = reachable ? m_maxLayerS + m_maxLayerT : 0;
    round.fromS.lastLayer = m_maxLayerS;
    round.fromT.lastLayer = m_maxLayerT;

    assert(!BI || !reachable || any_of(begin(m_graph), end(m_graph), [&](auto& node) {
        return distS(node)!=INF && distT(node)!=INF;
    }));
    for(auto& node : m_graph) {
        assert(distS(node)==INF || distS(node) <= m_maxLayerS);
        assert(distT(node)==INF || distT(node) <= m_maxLayerT);

        // forward search
        if(distS(node)<m_maxLayerS) {
            assert(distT(node)==INF);
            round.fromS.nodesBeforeLastLayer++;
            round.fromS.edgesBeforeLastLayer += std::size(node.edges);
        } else if(distS(node)==m_maxLayerS) {
            round.fromS.nodesInLastLayer++;
            round.fromS.edgesInLastLayer += std::size(node.edges);
        }

        // backward search
        if(distT(node)<m_maxLayerT) {
            assert(distS(node)==INF);
            round.fromT.nodesBeforeLastLayer++;
            round.fromT.edgesBeforeLastLayer += std::size(node.edges);
        } else if(distT(node)==m_maxLayerT) {
            round.fromT.nodesInLastLayer++;
            round.fromT.edgesInLastLayer += std::size(node.edges);
        }

        if(distS(node)==m_maxLayerS && distT(node)==m_maxLayerT) {
            round.nodesInInter++;
            round.edgesInInter += std::size(node.edges);
        }
    }

    m_stats.push_back(round);
}

template<bool BI, bool STAMPS, bool RESET, bool SKIP, bool ADJREV>
void Dinics<BI, STAMPS, RESET, SKIP, ADJREV>::gatherStatsDFS() {

    // accessors that handle stamps
    auto distS = [this](const node_type& node) {
        if constexpr (STAMPS) if(node.stamp!=m_time) return INF;
        return node.distS;
    };
    auto distT = [this](const node_type& node) {
        if constexpr (STAMPS) if(node.stamp!=m_time) return INF;
        if constexpr (BI) return node.distT;
        return INF;
    };

    RoundStatistics& round = m_stats.back();
    for(auto& node : m_graph) {
        // forward search
        if(distS(node)<m_maxLayerS) {
            round.dfsSpace.nodesInS += (node.nextEdge>0);
            round.dfsSpace.edgesFromS += node.nextEdge;
        } else if(distS(node)==m_maxLayerS) {
            round.dfsSpace.nodesInSlast += (node.nextEdge>0);
            round.dfsSpace.edgesFromSlast += node.nextEdge;
        }

        // backward search
        if(distT(node)<m_maxLayerT) {
            round.dfsSpace.nodesInT += (node.nextEdge>0);
            round.dfsSpace.edgesFromT += node.nextEdge;
        } else if(distT(node)==m_maxLayerT) {
            // dfs can never reach last layer of backwards search (except intersection)
        }

        if(distS(node)==m_maxLayerS && distT(node)==m_maxLayerT) {
            round.dfsSpace.nodesInInter += (node.nextEdge>0);
            round.dfsSpace.edgesFromInter += node.nextEdge;
        }
    }
}

// explicit instantiations
template class Dinics<0,0,0,0,0>;
template class Dinics<1,0,0,0,0>;
template class Dinics<0,1,0,0,0>;
template class Dinics<1,1,0,0,0>;
template class Dinics<0,0,1,0,0>;
template class Dinics<1,0,1,0,0>;
template class Dinics<0,1,1,0,0>;
template class Dinics<1,1,1,0,0>;
template class Dinics<0,0,0,1,0>;
template class Dinics<1,0,0,1,0>;
template class Dinics<0,1,0,1,0>;
template class Dinics<1,1,0,1,0>;
template class Dinics<0,0,1,1,0>;
template class Dinics<1,0,1,1,0>;
template class Dinics<0,1,1,1,0>;
template class Dinics<1,1,1,1,0>;
template class Dinics<0,0,0,0,1>;
template class Dinics<1,0,0,0,1>;
template class Dinics<0,1,0,0,1>;
template class Dinics<1,1,0,0,1>;
template class Dinics<0,0,1,0,1>;
template class Dinics<1,0,1,0,1>;
template class Dinics<0,1,1,0,1>;
template class Dinics<1,1,1,0,1>;
template class Dinics<0,0,0,1,1>;
template class Dinics<1,0,0,1,1>;
template class Dinics<0,1,0,1,1>;
template class Dinics<1,1,0,1,1>;
template class Dinics<0,0,1,1,1>;
template class Dinics<1,0,1,1,1>;
template class Dinics<0,1,1,1,1>;
template class Dinics<1,1,1,1,1>;


} // end namespace Dinics
