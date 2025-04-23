#pragma once

#include <robin-map/robin_map.h>

#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "dup_metrics.h"

using std::priority_queue;
using std::unordered_map;
using std::unordered_set;
using std::vector;

/*  */
class BamWrap;
class ReadEnds;
class ReadNameParser;

/*
 * optical duplicationgraph
 */
template <class Node>
struct Graph {           // set？
    vector<Node> nodes;  // 
    unordered_map<Node, int> nodeIdxMap;
    unordered_map<int, unordered_set<int>> neighbors;  // 

    int addNode(const Node &singleton) {
        int idx = -1;
        if (nodeIdxMap.find(singleton) == nodeIdxMap.end()) {
            nodes.push_back(singleton);
            idx = nodes.size() - 1;
            nodeIdxMap[singleton] = idx;
            neighbors[idx].clear();
        } else
            idx = nodeIdxMap[singleton];

        return idx;
    }

    /* bidirectional and public */
    void addEdge(Node &left, Node &right) {
        int leftIndex = addNode(left);
        if (left == right)
            return;
        int rightIndex = addNode(right);
        addNeighbor(leftIndex, rightIndex);
        addNeighbor(rightIndex, leftIndex);
    }

    void addNeighbor(int fromNode, int toNode) {
        unordered_set<int> &fromNodesNeighbors = neighbors[fromNode];
        if (fromNodesNeighbors.find(toNode) == fromNodesNeighbors.end())
            fromNodesNeighbors.insert(toNode);
    }

    /**
     * returns the cluster map of connected components
     *
     * @return Nodes that point to the same integer are in the same cluster.
     */
    void cluster(unordered_map<Node, int> *pClusterMap) {
        auto &clusterMap = *pClusterMap;
        vector<int> vCluster(nodes.size(), 0);
        for (int i = 0; i < vCluster.size(); ++i) vCluster[i] = i;
        for (int i = 0; i < neighbors.size(); ++i) {
            for (auto &j : neighbors[i]) joinNodes(j, i, &vCluster);
        }
        for (auto &node : nodes) {
            clusterMap[node] = findRepNode(vCluster, nodeIdxMap[node]);
        }
    }

    // Part of Union-Find with Path Compression that joins to nodes to be part of the same cluster.
    static void joinNodes(int nodeId1, int nodeId2, vector<int> *pGrouping) {
        auto &grouping = *pGrouping;
        int repNode1 = findRepNode(grouping, nodeId1);
        int repNode2 = findRepNode(grouping, nodeId2);
        if (repNode1 == repNode2)
            return;
        grouping[repNode1] = repNode2;
    }

    // Part of Union-Find with Path Compression to determine the duplicate set a particular UMI belongs to.
    static int findRepNode(vector<int> &grouping, int nodeId) {
        int representativeUmi = nodeId;  // All UMIs of a duplicate set will have the same reprsentativeUmi.
        while (representativeUmi != grouping[representativeUmi]) representativeUmi = grouping[representativeUmi];
        while (nodeId != representativeUmi) {
            int newUmiId = grouping[nodeId];
            grouping[nodeId] = representativeUmi;
            nodeId = newUmiId;
        }
        return representativeUmi;
    }
};

/*
 * read
 */
int16_t computeDuplicateScore(BamWrap &bw);

/*
 * Builds a read ends object that represents a single read.
 * read
 */
void buildReadEnds(BamWrap &bw, int64_t index, ReadNameParser &rnParser, ReadEnds *pKey);

/*
 * pairend read end
 */
void modifyPairedEnds(const ReadEnds &fragEnd, ReadEnds *pPairedEnds);

/**
 * Looks through the set of reads and identifies how many of the duplicates are
 * in fact optical duplicates, and stores the data in the instance level histogram.
 * Additionally sets the transient isOpticalDuplicate flag on each read end that is
 * identified as an optical duplicate.
 * 
 */
void trackOpticalDuplicates(vector<const ReadEnds *> &readEndsArr, const ReadEnds *pBestRe);

/*
 * Estimates the size of a library based on the number of paired end molecules observed
 * and the number of unique pairs observed.
 */
int64_t estimateLibrarySize(int64_t readPairs, int64_t uniqueReadPairs);

/**
 * Estimates the ROI (return on investment) that one would see if a library was sequenced to
 * x higher coverage than the observed coverage.
 **/
double estimateRoi(long estimatedLibrarySize, double x, long pairs, long uniquePairs);

/**
 * Calculates a histogram using the estimateRoi method to estimate the effective yield
 * doing x sequencing for x=1..10.
 */
void calculateRoiHistogram(DuplicationMetrics &metrics);