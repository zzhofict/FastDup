#include "md_funcs.h"

#include <math.h>
#include <stdint.h>
#include <util/bam_buf.h>
#include <util/murmur3.h>
#include <util/profiling.h>

#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "dup_metrics.h"
#include "md_args.h"
#include "read_ends.h"
#include "read_name_parser.h"

using std::cerr;
using std::endl;
using std::map;
using std::set;
using std::unordered_map;
using std::vector;

namespace nsgv {
extern MarkDupsArg gMdArg;
extern DuplicationMetrics gMetrics;
};

/*
 * 计算read的分数
 */
int16_t computeDuplicateScore(BamWrap &bw) {
    int16_t score = 0;
    switch (nsgv::gMdArg.DUPLICATE_SCORING_STRATEGY) {
    case nsmd::SUM_OF_BASE_QUALITIES:
        // two (very) long reads worth of high-quality bases can go over
        // Short.MAX_VALUE/2 and risk overflow.
        score += (int16_t)min(bw.GetSumOfBaseQualities(), INT16_MAX / 2);
        break;
    case nsmd::TOTAL_MAPPED_REFERENCE_LENGTH:
        if (!bw.GetReadUnmappedFlag())
            // no need to remember the score since this scoring mechanism is
            // symmetric
            score = (int16_t)min(bw.GetReferenceLength(), INT16_MAX / 2);
        break;
    case nsmd::RANDOM:
        // The RANDOM score gives the same score to both reads so that they get
        // filtered together. it's not critical do use the readName since the
        // scores from both ends get added, but it seem to be clearer this way.
        score += (short)(Murmur3::Instance().HashUnencodedChars(bw.query_name()) & 0b11111111111111);
        // subtract Short.MIN_VALUE/4 from it to end up with a number between
        // 0 and Short.MAX_VALUE/2. This number can be then discounted in case
        // the read is not passing filters. We need to stay far from overflow so
        // that when we add the two scores from the two read mates we do not
        // overflow since that could cause us to chose a failing read-pair
        // instead of a passing one.
        score -= INT16_MIN / 4;
    default:
        break;
    }
    // make sure that filter-failing records are heavily discounted. (the
    // discount can happen twice, once for each mate, so need to make sure we do
    // not subtract more than Short.MIN_VALUE overall.)
    score += bw.GetReadFailsVendorQualityCheckFlag() ? (int16_t)(INT16_MIN / 2) : 0;

    return score;
}

/*
 * Builds a read ends object that represents a single read.
 * 用来表示一个read的特征结构
 */
void buildReadEnds(BamWrap &bw, int64_t index, ReadNameParser &rnParser, ReadEnds *pKey) {
    auto &k = *pKey;
    auto &bc = bw.b->core;
    k.read1FirstOfPair = bw.GetFirstOfPairFlag();
    k.read1ReferenceIndex = bc.tid;
    k.read1Coordinate = (bc.flag & BAM_FREVERSE) ? bw.GetUnclippedEnd() : bw.GetUnclippedStart();
    k.orientation = (bc.flag & BAM_FREVERSE) ? ReadEnds::R : ReadEnds::F;

    k.read1IndexInFile = index;
    k.score = computeDuplicateScore(bw);
    // Doing this lets the ends object know that it's part of a pair
    if (bw.GetReadPairedFlag() && !bw.GetMateUnmappedFlag()) {
        k.read2ReferenceIndex = bc.mtid;
    }
    // Fill in the location information for optical duplicates
    if (!ReadNameParser::sWrongNameFormat)
        rnParser.AddLocationInformation(bw.query_name(), pKey);
    else
        nsgv::gMdArg.CHECK_OPTICAL_DUP = false;
    // cout << k.tile << ' ' << k.x << ' ' << k.y << endl;
    // 计算位置key
    k.posKey = BamWrap::bam_global_pos(k.read1ReferenceIndex, k.read1Coordinate);  // << 1 | k.orientation;
}

/* 对找到的pairend read end添加一些信息 */
void modifyPairedEnds(const ReadEnds &fragEnd, ReadEnds *pPairedEnds) {
    auto &pairedEnds = *pPairedEnds;

    int64_t bamIdx = fragEnd.read1IndexInFile;
    const int matesRefIndex = fragEnd.read1ReferenceIndex;
    const int matesCoordinate = fragEnd.read1Coordinate;
    // Set orientationForOpticalDuplicates, which always goes by the first then
    // the second end for the strands.  NB: must do this before updating the
    // orientation later.
    if (fragEnd.read1FirstOfPair) {
        pairedEnds.orientationForOpticalDuplicates =
            ReadEnds::GetOrientationByte(fragEnd.IsNegativeStrand(), pairedEnds.orientation == ReadEnds::R);
    } else {
        pairedEnds.orientationForOpticalDuplicates =
            ReadEnds::GetOrientationByte(pairedEnds.orientation == ReadEnds::R, fragEnd.IsNegativeStrand());
    }
    // If the other read is actually later, simply add the other read's data as
    // read2, else flip the reads
    if (matesRefIndex > pairedEnds.read1ReferenceIndex ||
        (matesRefIndex == pairedEnds.read1ReferenceIndex && matesCoordinate >= pairedEnds.read1Coordinate)) {
        pairedEnds.read2ReferenceIndex = matesRefIndex;
        pairedEnds.read2Coordinate = matesCoordinate;
        pairedEnds.read2IndexInFile = bamIdx;
        pairedEnds.orientation =
            ReadEnds::GetOrientationByte(pairedEnds.orientation == ReadEnds::R, fragEnd.IsNegativeStrand());

        // if the two read ends are in the same position, pointing in opposite
        // directions, the orientation is undefined and the procedure above will
        // depend on the order of the reads in the file. To avoid this, we set
        // it explicitly (to FR):
        if (pairedEnds.read2ReferenceIndex == pairedEnds.read1ReferenceIndex &&
            pairedEnds.read2Coordinate == pairedEnds.read1Coordinate && pairedEnds.orientation == ReadEnds::RF) {
            pairedEnds.orientation = ReadEnds::FR;
        }
    } else {
        pairedEnds.read2ReferenceIndex = pairedEnds.read1ReferenceIndex;
        pairedEnds.read2Coordinate = pairedEnds.read1Coordinate;
        pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
        pairedEnds.read1ReferenceIndex = matesRefIndex;
        pairedEnds.read1Coordinate = matesCoordinate;
        pairedEnds.read1IndexInFile = bamIdx;
        pairedEnds.orientation =
            ReadEnds::GetOrientationByte(fragEnd.IsNegativeStrand(), pairedEnds.orientation == ReadEnds::R);
        pairedEnds.posKey = fragEnd.posKey;
    }
    pairedEnds.score += fragEnd.score;
}

static inline bool closeEnough(const ReadEnds *lhs, const ReadEnds *rhs, const int distance) {
    return lhs != rhs &&  // no comparing an object to itself (checked using object identity)!
           (lhs->tile != ReadEnds::NO_VALUE) &&
           (rhs->tile != ReadEnds::NO_VALUE) &&  // no comparing objects without locations
           lhs->tile == rhs->tile &&             // and the same tile
           abs(lhs->x - rhs->x) <= distance && abs(lhs->y - rhs->y) <= distance;
}

static inline bool closeEnoughShort(const ReadEnds *lhs, const ReadEnds *rhs, const int distance) {
    return lhs != rhs && abs(lhs->x - rhs->x) <= distance && abs(lhs->y - rhs->y) <= distance;
}

/**
 * Finds which reads within the list of duplicates that are likely to be optical/co-localized duplicates of
 * one another. Within each cluster of optical duplicates that is found, one read remains un-flagged for
 * optical duplication and the rest are flagged as optical duplicates.  The set of reads that are considered
 * optical duplicates are indicated by returning "true" at the same index in the resulting boolean[] as the
 * read appeared in the input list of physical locations.
 *
 * @param list   a list of reads that are determined to be duplicates of one another
 * @param keeper a single PhysicalLocation that is the one being kept as non-duplicate, and thus should never be
 *               annotated as an optical duplicate. May in some cases be null, or a PhysicalLocation not
 *               contained within the list! (always not be null!)
 * @return a boolean[] of the same length as the incoming list marking which reads are optical duplicates
 */
static void findOpticalDuplicates(vector<const ReadEnds *> &readEndsArr, const ReadEnds *pBestRe,
                                  vector<bool> *pOpticalDuplicateFlags) {
    const int DEFAULT_OPTICAL_DUPLICATE_DISTANCE = 100;
    const int DEFAULT_MAX_DUPLICATE_SET_SIZE = 300000;

    vector<bool> &opticalDuplicateFlags = *pOpticalDuplicateFlags;
    // opticalDuplicateFlags.push_back(true); // for test
    int len = readEndsArr.size();
    // If there is only one or zero reads passed in (so there are obviously no optical duplicates),
    // or if there are too many reads (so we don't want to try to run this expensive n^2 algorithm),
    // then just return an array of all false
    if (len < 2 || len > DEFAULT_MAX_DUPLICATE_SET_SIZE) {
        return;
    }

    if (len >= 4) {
        /**
         * Compute the optical duplicates correctly in the case where the duplicate group could end up with transitive
         * optical duplicates
         * getOpticalDuplicatesFlagWithGraph
         */
        // Make a graph where the edges are reads that lie within the optical duplicate pixel distance from each other,
        // we will then use the union-find algorithm to cluster the graph and find optical duplicate groups
        Graph<int> opticalDistanceRelationGraph;
        unordered_map<int, vector<int>> tileRGmap;
        int keeperIndex = -1;
        for (int i = 0; i < readEndsArr.size(); ++i) {
            const ReadEnds *currentLoc = readEndsArr[i];
            if (currentLoc == pBestRe)
                keeperIndex = i;
            if (currentLoc->tile != ReadEnds::NO_VALUE) {
                int key = currentLoc->tile;  // 只处理一个样本，所以只有一个read group
                tileRGmap[key].push_back(i);
            }
            opticalDistanceRelationGraph.addNode(i);
        }
        // Since because finding adjacent optical duplicates is an O(n^2) operation, we can subdivide the input into its
        // readgroups in order to reduce the amount of redundant checks across readgroups between reads.
        // fillGraphFromAGroup
        for (auto &entry : tileRGmap) {
            auto &groupList = entry.second;
            if (groupList.size() > 1) {
                for (int i = 0; i < groupList.size(); ++i) {
                    int iIndex = groupList[i];
                    const ReadEnds *currentLoc = readEndsArr[iIndex];
                    for (int j = i + 1; j < groupList.size(); ++j) {
                        int jIndex = groupList[j];
                        const ReadEnds *other = readEndsArr[jIndex];
                        if (closeEnoughShort(currentLoc, other, DEFAULT_OPTICAL_DUPLICATE_DISTANCE))
                            opticalDistanceRelationGraph.addEdge(iIndex, jIndex);
                    }
                }
            }
        }
        // Keep a map of the reads and their cluster assignments
        unordered_map<int, int> opticalDuplicateClusterMap;
        opticalDistanceRelationGraph.cluster(&opticalDuplicateClusterMap);
        unordered_map<int, int> clusterToRepresentativeRead;
        int keeperCluster = -1;
        // Specially mark the keeper as specifically not a duplicate if it exists
        if (keeperIndex >= 0) {
            keeperCluster = opticalDuplicateClusterMap[keeperIndex];
            clusterToRepresentativeRead[keeperCluster] = keeperIndex;
        }

        for (auto &entry : opticalDuplicateClusterMap) {
            int recordIndex = entry.first;
            int recordAssignedCluster = entry.second;
            // If its not the first read we've seen for this cluster, mark it as an optical duplicate
            auto repReadItr = clusterToRepresentativeRead.find(recordAssignedCluster);
            if (repReadItr != clusterToRepresentativeRead.end() && recordIndex != keeperIndex) {
                const ReadEnds *representativeLoc = readEndsArr[repReadItr->second];
                const ReadEnds *currentRecordLoc = readEndsArr[recordIndex];
                // If not in the keeper cluster, then keep the minX -> minY valued duplicate (note the tile must be
                // equal for reads to cluster together)
                if (!(keeperIndex >= 0 && recordAssignedCluster == keeperCluster) &&
                    (currentRecordLoc->x < representativeLoc->x ||
                     (currentRecordLoc->x == representativeLoc->x && currentRecordLoc->y < representativeLoc->y))) {
                    // Mark the old min as an optical duplicate, and save the new min
                    opticalDuplicateFlags[repReadItr->second] = true;
                    clusterToRepresentativeRead[recordAssignedCluster] = recordIndex;
                } else {
                    // If a smaller read has already been visited, mark the test read as an optical duplicate
                    opticalDuplicateFlags[recordIndex] = true;
                }
            } else {
                clusterToRepresentativeRead[recordAssignedCluster] = recordIndex;
            }
        }
    } else {
        /**
         * Compute optical duplicates quickly in the standard case where we know that there won't be any transitive
         * distances to worry about. Note, this is guaranteed to be correct when there are at most 2x reads from a
         * readgroup or 3x with the keeper present
         * getOpticalDuplicatesFlagFast
         */
        // First go through and compare all the reads to the keeper
        for (int i = 0; i < len; ++i) {
            const ReadEnds *other = readEndsArr[i];
            opticalDuplicateFlags[i] = closeEnough(pBestRe, other, DEFAULT_OPTICAL_DUPLICATE_DISTANCE);
        }
        // Now go through and do each pairwise comparison not involving the actualKeeper
        for (int i = 0; i < len; ++i) {
            const ReadEnds *lhs = readEndsArr[i];
            if (lhs == pBestRe)  // no comparisons to actualKeeper since those are all handled above
                continue;
            for (int j = i + 1; j < len; ++j) {
                const ReadEnds *rhs = readEndsArr[j];
                if (rhs == pBestRe)  // no comparisons to actualKeeper since those are all handled above
                    continue;
                if (opticalDuplicateFlags[i] && opticalDuplicateFlags[j])
                    continue;  // both already marked, no need to check
                if (closeEnough(lhs, rhs, DEFAULT_OPTICAL_DUPLICATE_DISTANCE)) {
                    // At this point we want to mark either lhs or rhs as duplicate. Either could have been marked
                    // as a duplicate of the keeper (but not both - that's checked above), so be careful about which
                    // one to now mark as a duplicate.
                    int index = opticalDuplicateFlags[j] ? i : j;
                    opticalDuplicateFlags[index] = true;
                }
            }
        }
    }
}

/**
 * Looks through the set of reads and identifies how many of the duplicates are
 * in fact optical duplicates, and stores the data in the instance level histogram.
 *
 * We expect only reads with FR or RF orientations, not a mixture of both.
 *
 * In PCR duplicate detection, a duplicates can be a have FR and RF when fixing the orientation order to the first end
 * of the mate.  In optical duplicate detection, we do not consider them duplicates if one read as FR and the other RF
 * when we order orientation by the first mate sequenced (read #1 of the pair).
 */
static int checkOpticalDuplicates(vector<const ReadEnds *> &readEndsArr, const ReadEnds *pBestRe) {
    vector<bool> opticalDuplicateFlags(readEndsArr.size(), false);
    // find OpticalDuplicates
    findOpticalDuplicates(readEndsArr, pBestRe, &opticalDuplicateFlags);
    int opticalDuplicates = 0;
    for (int i = 0; i < opticalDuplicateFlags.size(); ++i) {
        ReadEnds *pRe = const_cast<ReadEnds *>(readEndsArr[i]);
        if (opticalDuplicateFlags[i]) {
            ++opticalDuplicates;
            pRe->isOpticalDuplicate = true;
        } else {
            pRe->isOpticalDuplicate = false;
        }
    }
    return opticalDuplicates;
}

/**
 * 记录光学原因造成的冗余
 */
void trackOpticalDuplicates(vector<const ReadEnds *> &readEndsArr, const ReadEnds *pBestRe) {
    bool hasFR = false, hasRF = false;
    // Check to see if we have a mixture of FR/RF
    for (auto pRe : readEndsArr) {
        if (ReadEnds::FR == pRe->orientationForOpticalDuplicates)
            hasFR = true;
        else if (ReadEnds::RF == pRe->orientationForOpticalDuplicates)
            hasRF = true;
    }

    // Check if we need to partition since the orientations could have changed
    int nOpticalDup;
    if (hasFR && hasRF) {  // need to track them independently
        vector<const ReadEnds *> trackOpticalDuplicatesF;
        vector<const ReadEnds *> trackOpticalDuplicatesR;
        // Split into two lists: first of pairs and second of pairs,
        // since they must have orientation and same starting end
        for (auto pRe : readEndsArr) {
            if (ReadEnds::FR == pRe->orientationForOpticalDuplicates)
                trackOpticalDuplicatesF.push_back(pRe);
            else if (ReadEnds::RF == pRe->orientationForOpticalDuplicates)
                trackOpticalDuplicatesR.push_back(pRe);
            else
                cerr << "Found an unexpected orientation: " << pRe->orientation << endl;
        }
        // track the duplicates
        int nOpticalDupF = checkOpticalDuplicates(trackOpticalDuplicatesF, pBestRe);
        int nOpticalDupR = checkOpticalDuplicates(trackOpticalDuplicatesR, pBestRe);
        nOpticalDup = nOpticalDupF + nOpticalDupR;
    } else {  // No need to partition
        nOpticalDup = checkOpticalDuplicates(readEndsArr, pBestRe);
    }

    // 统计信息，trackDuplicateCounts
    ++nsgv::gMetrics.DuplicateCountHist[readEndsArr.size()];
    if (readEndsArr.size() > nOpticalDup)
        ++nsgv::gMetrics.NonOpticalDuplicateCountHist[readEndsArr.size() - nOpticalDup];
    if (nOpticalDup)
        ++nsgv::gMetrics.OpticalDuplicatesCountHist[nOpticalDup + 1];
}

/**
 * Estimates the size of a library based on the number of paired end molecules observed
 * and the number of unique pairs observed.
 * <p>
 * Based on the Lander-Waterman equation that states:
 * C/X = 1 - exp( -N/X )
 * where
 * X = number of distinct molecules in library
 * N = number of read pairs
 * C = number of distinct fragments observed in read pairs
 */
int64_t estimateLibrarySize(int64_t readPairs, int64_t uniqueReadPairs) {
    int64_t librarySize = 0;
    int64_t readPairDuplicates = readPairs - uniqueReadPairs;
    auto f = [](double x, double c, double n) { return c / x - 1 + exp(-n / x); };
    if (readPairs > 0 && readPairDuplicates > 0) {
        double m = 1.0;
        double M = 100.0;
        if (uniqueReadPairs >= readPairs || f(m * uniqueReadPairs, uniqueReadPairs, readPairs) < 0) {
            cerr << "Invalid values for pairs and unique pairs: " << readPairs << ", " << uniqueReadPairs << endl;
            return 0;
        }
        // find value of M, large enough to act as other side for bisection method
        while (f(M * uniqueReadPairs, uniqueReadPairs, readPairs) > 0) {
            M *= 10.0;
        }
        // use bisection method (no more than 40 times) to find solution
        for (int i = 0; i < 40; ++i) {
            double r = (m + M) / 2.0;
            double u = f(r * uniqueReadPairs, uniqueReadPairs, readPairs);
            if (u == 0)
                break;
            else if (u > 0)
                m = r;
            else if (u < 0)
                M = r;
        }
        return uniqueReadPairs * (m + M) / 2.0;
    }
    return librarySize;
}

/**
 * Estimates the ROI (return on investment) that one would see if a library was sequenced to
 * x higher coverage than the observed coverage.
 *
 * @param estimatedLibrarySize the estimated number of molecules in the library
 * @param x                    the multiple of sequencing to be simulated (i.e. how many X sequencing)
 * @param pairs                the number of pairs observed in the actual sequencing
 * @param uniquePairs          the number of unique pairs observed in the actual sequencing
 * @return a number z <= x that estimates if you had pairs*x as your sequencing then you
 * would observe uniquePairs*z unique pairs.
 */
double estimateRoi(long estimatedLibrarySize, double x, long pairs, long uniquePairs) {
    return estimatedLibrarySize * (1 - exp(-(x * pairs) / estimatedLibrarySize)) / uniquePairs;
}

/**
 * Calculates a histogram using the estimateRoi method to estimate the effective yield
 * doing x sequencing for x=1..10.
 */
void calculateRoiHistogram(DuplicationMetrics &metrics) {
    int64_t uniquePairs = metrics.READ_PAIRS_EXAMINED - metrics.READ_PAIR_DUPLICATES;
    metrics.CoverageMult.resize(101);
    for (int x = 1; x <= 100; ++x) {
        metrics.CoverageMult[x] =
            estimateRoi(metrics.ESTIMATED_LIBRARY_SIZE, x, metrics.READ_PAIRS_EXAMINED, uniquePairs);
    }
}