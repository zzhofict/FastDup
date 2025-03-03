/*
Description: read
ends结构体主要用来标记冗余，包含一些序列的测序过程中的物理信息等

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

/* 包含了所有read ends信息，如picard里边的 ReadEndsForMarkDuplicates*/
struct ReadEnds : PhysicalLocation {
    static const int8_t F = 0, R = 1, FF = 2, FR = 3, RR = 4, RF = 5;
    /* 保留一些bam记录中的数据 */
    bool read1FirstOfPair = true;
    /* ReadEnds中的成员变量 */
    /** Little struct-like class to hold read pair (and fragment) end data for
     * duplicate marking. */
    // int16_t libraryId; // 没用，不考虑多样本
    int8_t orientation = -1;
    int32_t read1ReferenceIndex = -1;
    int32_t read1Coordinate = -1;
    int32_t read2ReferenceIndex = -1;
    // This field is overloaded for flow based processing as the end coordinate of read 1. (paired reads not supported)
    int32_t read2Coordinate = -1;
    /* Additional information used to detect optical dupes */
    // int16_t readGroup = -1; 一般经过比对后的bam文件只有一个read
    // group，normal或者tumor
    /** For optical duplicate detection the orientation matters regard to 1st or
     * 2nd end of a mate */
    int8_t orientationForOpticalDuplicates = -1;
    /** A *transient* flag marking this read end as being an optical duplicate.
     */
    bool isOpticalDuplicate = false;

    /* ReadEndsForMarkDuplicates中的成员变量 */
    /** Little struct-like class to hold read pair (and fragment) end data for
     * MarkDuplicatesWithMateCigar **/
    int16_t score = 0;
    int64_t read1IndexInFile = -1;
    int64_t read2IndexInFile = -1;
    int64_t duplicateSetSize = -1;

    /* ReadEndsForMarkDuplicatesWithBarcodes中的成员变量 (好像用不到) */
    // int32_t barcode = 0; // primary barcode for this read (and pair)
    // int32_t readOneBarcode = 0; // read one barcode, 0 if not present
    // int32_t readTwoBarcode = 0; // read two barcode, 0 if not present or not
    // paired

    /* zzh增加的成员变量 */
    int64_t posKey = -1;  // 根据位置信息生成的关键字 return (int64_t)tid <<
                          // MAX_CONTIG_LEN_SHIFT | (int64_t)pos; （包含clip的序列，也就是可能比map结果更偏左）

    /* 用来做一些判断，因为一些readends会做多次操作，比如task之间有重叠等等 */
    int oprateTime = 0;

    /* 根据pairend read的比对方向，来确定整体的比对方向 */
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

    /* 比较两个readends是否一样（有个冗余） */
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

    /* 比对方向是否正向 */
    bool IsPositiveStrand() const { return orientation == F; }

    /* pairend是否合适的比对上了 */
    bool IsPaired() const { return read2ReferenceIndex != -1; }

    bool IsNegativeStrand() const { return orientation == R; }

    // 对于相交的数据进行比对，a是否小于b，根据AreComparableForDuplicates函数得来
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
        if (comp == 0)  // 这个放在坐标比较之前，因为在解析bam的时候，可能有的给read2ReferenceIndex赋值了,有的则没赋值
            comp = a.orientation - b.orientation;
        if (comp == 0)
            comp = a.read2ReferenceIndex - b.read2ReferenceIndex;
        if (comp == 0)
            comp = a.read2Coordinate - b.read2Coordinate;
        if (comp == 0)
            comp = a.tile - b.tile;
        if (comp == 0)
            comp = a.x - b.x; // 由于picard的bug，用short类型来表示x，y，导致其可能为负数
        if (comp == 0)
            comp - a.y - b.y;
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
        if (comp == 0)  // 这个放在坐标比较了之后，把坐标范围的放在之前，这样对分段数据块比较好处理
            comp = a.orientation - b.orientation;
        if (comp == 0)
            comp = a.tile - b.tile;
        if (comp == 0)
            comp = a.x - b.x; // 由于picard的bug，用short类型来表示x，y，导致其可能为负数
        if (comp == 0)
            comp - a.y - b.y;
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
        if (comp == 0)  // 这个放在坐标比较了之后，把坐标范围的放在之前，这样对分段数据块比较好处理
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