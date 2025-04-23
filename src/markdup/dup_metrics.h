#pragma once

#include <stdint.h>

#include <string>
#include <vector>

#include "md_types.h"

using std::string;
using std::vector;

/*
 * 
 */
struct DuplicationMetrics {
    /**
     * The library on which the duplicate marking was performed.
     */
    string LIBRARY = "";
    /**
     * The number of mapped reads examined which did not have a mapped mate pair,
     * either because the read is unpaired, or the read is paired to an unmapped mate.
     */
    uint64_t UNPAIRED_READS_EXAMINED = 0;
    /**
     * The number of mapped read pairs examined. (Primary, non-supplemental)
     */
    uint64_t READ_PAIRS_EXAMINED = 0;
    /**
     * The number of reads that were either secondary or supplementary
     */
    uint64_t SECONDARY_OR_SUPPLEMENTARY_RDS = 0;
    /**
     * The total number of unmapped reads examined. (Primary, non-supplemental)
     */
    uint64_t UNMAPPED_READS = 0;
    /**
     * The number of fragments that were marked as duplicates.
     */
    uint64_t UNPAIRED_READ_DUPLICATES = 0;
    /**
     * The number of read pairs that were marked as duplicates.
     */
    uint64_t READ_PAIR_DUPLICATES = 0;
    /**
     * The number of read pairs duplicates that were caused by optical duplication.
     * Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
     */
    uint64_t READ_PAIR_OPTICAL_DUPLICATES = 0;
    /**
     * The fraction of mapped sequence that is marked as duplicate.
     */
    double PERCENT_DUPLICATION = 0.0;
    /**
     * The estimated number of unique molecules in the library based on PE duplication.
     */
    uint64_t ESTIMATED_LIBRARY_SIZE = 0;

    // 
    vector<double> CoverageMult;

    // ，4，4
    MDMap DuplicateCountHist;
    MDMap NonOpticalDuplicateCountHist;
    MDMap OpticalDuplicatesCountHist;

    //  for test
    MDSet<int64_t> singletonReads;
    MDSet<int64_t> dupReads;  // for test
    MDSet<int64_t> bestReads;
};