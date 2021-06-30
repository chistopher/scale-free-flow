
#pragma once

#include <utils.h>

template <typename, typename, typename> class Graph; // FWD

namespace BoykovKolmogorov {

    class BKAlgorithm {
        Graph<double,double,double>* m_graph;
        bool first = true;
    public:
        BKAlgorithm(const EdgeList& graph);
        ~BKAlgorithm();
        double maxFlow(int s, int t);
    };

} // end namespace BoykovKolmogorov

