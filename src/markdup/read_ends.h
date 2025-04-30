/*
Description: read
ends，

Copyright : All right reserved by ICT

Author : Zhang Zhonghai
Date : 2023/11/3
*/
#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <util/bam_wrap.h>

#include <algorithm>

/**
 * Small interface that provides access to the physical location information
 * about a cluster. All values should be defaulted to -1 if unavailable.
 * ReadGroup and Tile should only allow non-zero positive integers, x and y
 * coordinates may be negative.
 */
struct PhysicalLocation {
    static const int NO_VALUE = -1;
    /**
     * Small class that provides access to the physical location information
     * about a cluster. All values should be defaulted to -1 if unavailable.
     * Tile should only allow non-zero positive integers, x and y coordinates
     * must be non-negative. This is different from PhysicalLocationShort in
     * that the x and y positions are ints, not shorts thus, they do not
     * overflow within a HiSeqX tile.
     */
    int16_t tile = -1;
    // int32_t x = -1;
    // int32_t y = -1;
    // This is a bug in picard Markduplicates, because some tile coordinates exceede the range of int16_t
    int16_t x = -1;
    int16_t y = -1;
};

/* read ends，picard ReadEndsForMarkDuplicates*/
struct ReadEnds : PhysicalLocation {
    static const int8_t F = 0, R = 1, FF = 2, FR = 3, RR = 4, RF = 5;
    /* bam */
    bool read1FirstOfPair = true;
    /* ReadEnds */
    /** Little struct-like class to hold read pair (and fragment) end data for
     * duplicate marking. */
    // int16_t libraryId; // ，
    int8_t orientation = -1;
    int32_t read1ReferenceIndex = -1;
    int32_t read1Coordinate = -1;
    int32_t read2ReferenceIndex = -1;
    // This field is overloaded for flow based processing as the end coordinate of read 1. (paired reads not supported)
    int32_t read2Coordinate = -1;
    /* Additional information used to detect optical dupes */
    // int16_t readGroup = -1; bamread
    // group，normaltumor
    /** For optical duplicate detection the orientation matters regard to 1st or
     * 2nd end of a mate */
    int8_t orientationForOpticalDuplicates = -1;
    /** A *transient* flag marking this read end as being an optical duplicate.
     */
    bool isOpticalDuplicate = false;

    /* ReadEndsForMarkDuplicates */
    /** Little struct-like class to hold read pair (and fragment) end data for
     * MarkDuplicatesWithMateCigar **/
    int16_t score = 0;
    int64_t read1IndexInFile = -1;
    int64_t read2IndexInFile = -1;
    int64_t duplicateSetSize = -1;

    /* ReadEndsForMarkDuplicatesWithBarcodes () */
    // int32_t barcode = 0; // primary barcode for this read (and pair)
    // int32_t readOneBarcode = 0; // read one barcode, 0 if not present
    // int32_t readTwoBarcode = 0; // read two barcode, 0 if not present or not
    // paired

    /* zzh */
    int64_t posKey = -1;  //  return (int64_t)tid <<
                          // MAX_CONTIG_LEN_SHIFT | (int64_t)pos; （clip，map）

    /* ，readends，task */
    int oprateTime = 0;

    /* pairend read， */
    static int8_t GetOrientationByte(bool read1NegativeStrand, bool read2NegativeStrand) {
        if (read1NegativeStrand) {
            if (read2NegativeStrand)
                return RR;
            else
                return RF;
        } else {
            if (read2NegativeStrand)
                return FR;
            else
                return FF;
        }
    }

    /* readends */
    static bool AreComparableForDuplicates(const ReadEnds &lhs, const ReadEnds &rhs, bool compareRead2) {
        bool areComparable = true;
        areComparable = lhs.read1ReferenceIndex == rhs.read1ReferenceIndex &&
                        lhs.read1Coordinate == rhs.read1Coordinate && lhs.orientation == rhs.orientation;
        if (areComparable && compareRead2) {
            areComparable =
                lhs.read2ReferenceIndex == rhs.read2ReferenceIndex && lhs.read2Coordinate == rhs.read2Coordinate;
        }
        return areComparable;
    }

    /*  */
    bool IsPositiveStrand() const { return orientation == F; }

    /* pairend */
    bool IsPaired() const { return read2ReferenceIndex != -1; }

    bool IsNegativeStrand() const { return orientation == R; }

    // AreComparableForDuplicates
    static inline bool ReadLittleThan(const ReadEnds &a, const ReadEnds &b, bool compareRead2 = false) {
        int comp = a.read1ReferenceIndex - b.read1ReferenceIndex;
        if (comp == 0)
            comp = a.read1Coordinate - b.read1Coordinate;
        if (compareRead2) {
            if (comp == 0)
                comp = a.read2ReferenceIndex - b.read2ReferenceIndex;
            if (comp == 0)
                comp = a.read2Coordinate - b.read2Coordinate;
        }
        if (comp == 0)
            comp = a.orientation - b.orientation;

        return comp < 0;
    }

    static bool FragLittleThan(const ReadEnds &a, const ReadEnds &b) {
        int comp = a.read1ReferenceIndex - b.read1ReferenceIndex;
        if (comp == 0)
            comp = a.read1Coordinate - b.read1Coordinate;
        if (comp == 0)  // ，bam，read2ReferenceIndex,
            comp = a.orientation - b.orientation;
        if (comp == 0)
            comp = a.read2ReferenceIndex - b.read2ReferenceIndex;
        if (comp == 0)
            comp = a.read2Coordinate - b.read2Coordinate;
//        if (comp == 0)
//            comp = a.tile - b.tile;
//        if (comp == 0)
//            comp = a.x - b.x; // picardbug，shortx，y，
//        if (comp == 0)
//            comp - a.y - b.y;
        if (comp == 0)
            comp = (int)(a.read1IndexInFile - b.read1IndexInFile);
        if (comp == 0)
            comp = (int)(a.read2IndexInFile - b.read2IndexInFile);
        return comp < 0;
    }

    static bool PairLittleThan(const ReadEnds &a, const ReadEnds &b) {
        int comp = a.read1ReferenceIndex - b.read1ReferenceIndex;
        if (comp == 0)
            comp = a.read1Coordinate - b.read1Coordinate;
        if (comp == 0)
            comp = a.read2ReferenceIndex - b.read2ReferenceIndex;
        if (comp == 0)
            comp = a.read2Coordinate - b.read2Coordinate;
        if (comp == 0)  // ，，
            comp = a.orientation - b.orientation;
//        if (comp == 0)
//            comp = a.tile - b.tile;
//        if (comp == 0)
//            comp = a.x - b.x; // picardbug，shortx，y，
//        if (comp == 0)
//            comp - a.y - b.y;
        if (comp == 0)
            comp = (int)(a.read1IndexInFile - b.read1IndexInFile);
        if (comp == 0)
            comp = (int)(a.read2IndexInFile - b.read2IndexInFile);
        return comp < 0;
    }

    static bool CorLittleThan(const ReadEnds &a, const ReadEnds &b) {
        int comp = a.read1ReferenceIndex - b.read1ReferenceIndex;
        if (comp == 0)
            comp = a.read1Coordinate - b.read1Coordinate;
        if (comp == 0)
            comp = a.read2ReferenceIndex - b.read2ReferenceIndex;
        if (comp == 0)
            comp = a.read2Coordinate - b.read2Coordinate;
        if (comp == 0)  // ，，
            comp = a.orientation - b.orientation;
        return comp < 0;
    }

    // for pairs only
    int64_t Left() { return BamWrap::bam_global_pos(read1ReferenceIndex, read1Coordinate); }
    int64_t Right() { return BamWrap::bam_global_pos(read2ReferenceIndex, read2Coordinate); }
};

struct ReadEndsHash {
    std::size_t operator()(const ReadEnds &o) const { return std::hash<int64_t>()(o.read1IndexInFile); }
};

struct ReadEndsEqual {
    bool operator()(const ReadEnds &o1, const ReadEnds &o2) const { return o1.read1IndexInFile == o2.read1IndexInFile; }
};
