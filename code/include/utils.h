
#pragma once

#include <vector>
#include <utility>
#include <string>
#include <array>
#include <map>
#include <random>

struct EdgeList {
    std::vector<std::pair<int,int>> edges;
    std::vector<double> weights;
    bool isDirected = false;
    bool hasWeights = false;
    size_t n = 0;
};

// reads a graph from a file
EdgeList readGraph(const std::string& file);

// returns vector of degrees
std::vector<unsigned> degrees(const EdgeList& edges);

// find largest connected component and returns it as an edge list (ids are remapped) along with its size
EdgeList largestCC(const EdgeList& graph, std::map<int,int>* out = nullptr);

// for each undirected edge, creates two directed ones
EdgeList undirectedToDirected(const EdgeList& graph);

// import and export EdgeList to DIMACS; a -1 for s or t means no s-t pair is given
void exportDIMACS(const EdgeList& g, const std::string& file, int s=-1, int t=-1);
std::tuple<EdgeList,int,int> importDIMACS(const std::string& file);

// print graph stats
void printGraphStats(const EdgeList& graph, const std::string& file = "");

// choose random terminal pairs in degree interval
std::vector<std::pair<int,int>> generatePairs(const EdgeList& graph, int num, int mndeg, int mxdeg, std::mt19937& gen);
