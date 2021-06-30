
#include <iostream>

namespace Dinics {

    struct RoundStatistics {
        static constexpr auto CSVHEADER = "distOfSink,flowBefore," // general
                "lastLayerS,nodesBeforeLastLayerS,nodesInLastLayerS,edgesBeforeLastLayerS,edgesInLastLayerS," // forward BFS
                "lastLayerT,nodesBeforeLastLayerT,nodesInLastLayerT,edgesBeforeLastLayerT,edgesInLastLayerT," // backward BFS
                "nodesInInter,edgesInInter," // BFS intersection
                "dfsNodesBeforeLastLayerS,dfsNodesInLastLayerS,dfsNodesInInter,dfsNodesBeforeLastLayerT," // DFS nodes
                "dfsEdgesBeforeLastLayerS,dfsEdgesInLastLayerS,dfsEdgesInInter,dfsEdgesBeforeLastLayerT"; // DFS edges

        // general stuff
        long long distOfSink = 0; // dist from s to t
        long long flowBefore = 0;

        // BFS search space
        struct SearchSpaceBFS {
            long long lastLayer = 0;
            long long nodesBeforeLastLayer = 0; // nodes with fully explored adjacent edges
            long long nodesInLastLayer = 0; // nodes for which we do not check all edges in BFS
            long long edgesBeforeLastLayer = 0; // edges adjacent to nodes before last BFS layer
            long long edgesInLastLayer = 0; // edges adjacent to nodes in last BFS layer
        };
        SearchSpaceBFS fromS;
        SearchSpaceBFS fromT;
        long long nodesInInter= 0; // intersection of forward and backward search
        long long edgesInInter = 0;

        // DFS search space
        struct SearchSpaceDFS { // excludes the actual augmented paths
            long long nodesInS = 0;
            long long nodesInSlast = 0;
            long long nodesInInter = 0;
            long long nodesInT = 0;

            long long edgesFromS = 0;
            long long edgesFromSlast = 0;
            long long edgesFromInter = 0;
            long long edgesFromT = 0;
        };
        SearchSpaceDFS dfsSpace;
    };

    // writes all values comma separated to os according to RoundStatistics::CSVHEADER
    std::ostream& operator<<(std::ostream& os, const RoundStatistics& s);

} // end namespace Dinics
