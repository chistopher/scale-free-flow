
#include <DinicsStats.h>

namespace Dinics {

std::ostream &operator<<(std::ostream& os, const RoundStatistics &s) {
    os << s.distOfSink << ',' << s.flowBefore << ',';
    for(auto& search : {s.fromS, s.fromT}) {
        os << search.lastLayer            << ','
           << search.nodesBeforeLastLayer << ','
           << search.nodesInLastLayer     << ','
           << search.edgesBeforeLastLayer << ','
           << search.edgesInLastLayer     << ',';
    }
    os << s.nodesInInter << ',' << s.edgesInInter << ',';
    os << s.dfsSpace.nodesInS   << ',' << s.dfsSpace.nodesInSlast   << ',' << s.dfsSpace.nodesInInter   << ',' << s.dfsSpace.nodesInT << ',';
    os << s.dfsSpace.edgesFromS << ',' << s.dfsSpace.edgesFromSlast << ',' << s.dfsSpace.edgesFromInter << ',' << s.dfsSpace.edgesFromT;
    return os;
}

} // end namespace Dinics
