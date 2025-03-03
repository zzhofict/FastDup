#pragma once

#include <robin-map/robin_map.h>
#include <robin-map/robin_set.h>

#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "read_ends.h"
#include "util/bam_buf.h"

using std::priority_queue;
using std::set;
using std::string;
using std::unordered_set;
using std::vector;

/* 存放未匹配readend相同位点的所有readend */
struct UnpairedREInfo {
    int64_t taskSeq;  // 对应第几轮计算
    ReadEnds unpairedRE;
};

/* 对于一个pair数据，一个完整的计算点，包含read1的比对位置和read2的比对位置 */
struct CalcKey {
    int8_t orientation = -1;
    int32_t read1ReferenceIndex = -1;
    int32_t read1Coordinate = -1;
    int32_t read2ReferenceIndex = -1;
    int32_t read2Coordinate = -1;
    int64_t left = -1;
    int64_t right = -1;

    CalcKey() {}
    static CalcKey MaxCK() {
        CalcKey ck;
        ck.left = ck.right = INT64_MAX;
        return ck;
    }
    static CalcKey MinCK() {
        CalcKey ck;
        ck.left = ck.right = INT64_MIN;
        return ck;
    }

    CalcKey(const ReadEnds &re) {
        orientation = re.orientation;
        read1ReferenceIndex = re.read1ReferenceIndex;
        read1Coordinate = re.read1Coordinate;
        read2ReferenceIndex = re.read2ReferenceIndex;
        read2Coordinate = re.read2Coordinate;
        left = Read1Pos();
        right = Read2Pos();
    }

    int64_t Read1Pos() const { return BamWrap::bam_global_pos(read1ReferenceIndex, read1Coordinate); }
    int64_t Read2Pos() const { return BamWrap::bam_global_pos(read2ReferenceIndex, read2Coordinate); }
    int64_t Left() const { return left; }
    int64_t Right() const { return right; }

    bool operator<(const CalcKey &o) const {
        int comp = read1ReferenceIndex - o.read1ReferenceIndex;
        if (comp == 0)
            comp = read1Coordinate - o.read1Coordinate;
        if (comp == 0)
            comp = read2ReferenceIndex - o.read2ReferenceIndex;
        if (comp == 0)
            comp = read2Coordinate - o.read2Coordinate;
        // 需要orientation，因为要跟排序的比较方式和顺序一致
        if (comp == 0)
            comp = orientation - o.orientation;
        return comp < 0;
    }
    bool operator <= (const CalcKey &o) const {
        return *this < o || *this == o;
    }
    bool operator==(const CalcKey &o) const {
        return read1ReferenceIndex == o.read1ReferenceIndex && read1Coordinate == o.read1Coordinate &&
               orientation == o.orientation && read2ReferenceIndex == o.read2ReferenceIndex &&
               read2Coordinate == o.read2Coordinate;
    }
    bool operator<(const ReadEnds &o) const {
        int comp = read1ReferenceIndex - o.read1ReferenceIndex;
        if (comp == 0)
            comp = read1Coordinate - o.read1Coordinate;
        if (comp == 0)
            comp = read2ReferenceIndex - o.read2ReferenceIndex;
        if (comp == 0)
            comp = read2Coordinate - o.read2Coordinate;
        // 需要orientation，因为要跟排序的比较方式和顺序一致
        if (comp == 0)
            comp = orientation - o.orientation;
        return comp < 0;
    }
    std::size_t operator()() const {
        size_t h1 = read1ReferenceIndex;
        h1 = (h1 << 40) | (read1Coordinate << 8) | orientation;
        size_t h2 = read2ReferenceIndex;
        h2 = (h2 << 32) | read2Coordinate;
        return std::hash<int64_t>()(h1) ^ std::hash<int64_t>()(h2);
    }
};

struct CalcKeyHash {
    std::size_t operator()(const CalcKey &o) const { return o(); }
};

struct CalcKeyEqual {
    bool operator()(const CalcKey &o1, const CalcKey &o2) const { return o1 == o2; }
};

/* 用来记录冗余索引相关的信息 */
struct DupInfo {
    int16_t dupSet = 0;  // dup set size
    uint16_t repIdxHigh = 0;  // 这一批冗余中的非冗余read 代表的索引
    uint32_t repIdxLow = 0;
    int64_t idx;

    DupInfo() : DupInfo(-1, 0, 0) {}
    DupInfo(int64_t idx_) : DupInfo(idx_, 0, 0) {}
    DupInfo(int64_t idx_, int64_t repIdx_, int dupSet_) : idx(idx_), dupSet(dupSet_) {
        repIdxHigh = repIdx_ >> 32;
        repIdxLow = (uint32_t)repIdx_;
    }
    int64_t GetRepIdx() {
        int64_t repIdx = repIdxHigh;
        repIdx = (repIdx << 32) | repIdxLow;
        return repIdx;
    }
    bool operator<(const DupInfo &o) const { return idx < o.idx; }
    bool operator>(const DupInfo &o) const { return idx > o.idx; }
    operator int64_t() const { return idx; }
};

struct DupInfoHash {
    std::size_t operator()(const DupInfo &o) const { return std::hash<int64_t>()(o.idx); }
};

struct DupInfoEqual {
    bool operator()(const DupInfo &o1, const DupInfo &o2) const { return o1.idx == o2.idx; }
    bool operator()(const DupInfo &o1, const int64_t &o2) const { return o1.idx == o2; }
    bool operator()(const int64_t &o1, const DupInfo &o2) const { return o1 == o2.idx; }
};

using MDMap = tsl::robin_map<int64_t, int64_t>;

template <typename T>
// using MDSet = set<T>;
// using MDSet = unordered_set<T>;
using MDSet = tsl::robin_set<T>;

template <typename T>
// using DPSet = set<T>;
// using DPSet = unordered_set<T, DupInfoHash, DupInfoEqual>;
using DPSet = tsl::robin_set<T, DupInfoHash, DupInfoEqual>;

template <typename T>
using CalcSet = set<T>;
// using CalcSet = tsl::robin_set<T, CalcKeyHash>;

using ReadEndsSet = tsl::robin_set<ReadEnds, ReadEndsHash, ReadEndsEqual>;

/* 当遗留数据在当前任务找到了pair read后，进行冗余计算时候存放结果的数据结构 */
struct TaskSeqDupInfo {
    DPSet<DupInfo> dupIdx;
    MDSet<int64_t> opticalDupIdx;
    DPSet<DupInfo> repIdx;
    MDSet<int64_t> notDupIdx;
    MDSet<int64_t> notOpticalDupIdx;
    MDSet<int64_t> notRepIdx;
};

/* 保存有未匹配pair位点的信息，包括read end数组和有几个未匹配的read end */
struct UnpairedPosInfo {
    int unpairedNum = 0;
    int64_t taskSeq;
    vector<ReadEnds> pairArr;
    MDSet<string> readNameSet;
};
// typedef unordered_map<string, UnpairedREInfo> UnpairedNameMap;
// typedef unordered_map<int64_t, UnpairedPosInfo> UnpairedPositionMap;

typedef tsl::robin_map<string, UnpairedREInfo> UnpairedNameMap;  // 以read name为索引，保存未匹配的pair read
//typedef map<string, UnpairedREInfo> UnpairedNameMap;  // 以read name为索引，保存未匹配的pair read
typedef tsl::robin_map<int64_t, UnpairedPosInfo> UnpairedPositionMap;  // 以位点为索引，保存该位点包含的对应的所有read和该位点包含的剩余未匹配的read的数量
// typedef map<CalcKey, vector<ReadEnds>> CkeyReadEndsMap;  // 以calckey为关键字，保存在相邻数据块之前找到的匹配readEnds
typedef unordered_map<CalcKey, vector<ReadEnds>, CalcKeyHash, CalcKeyEqual>
    CkeyReadEndsMap;  // 以calckey为关键字，保存在相邻数据块之前找到的匹配readEnds
// typedef tsl::robin_map<CalcKey, vector<ReadEnds>, CalcKeyHash, CalcKeyEqual> CkeyReadEndsMap;  //
// 以calckey为关键字，保存在相邻数据块之前找到的匹配readEnds

/* 单线程处理冗余参数结构体 */
struct MarkDupDataArg {
    int64_t taskSeq;                     // 任务序号
    int64_t bamStartIdx;                 // 当前vBam数组中第一个bam记录在整体bam中所处的位置
    vector<BamWrap *> bams;              // 存放待处理的bam read
    vector<ReadEnds> pairs;              // 成对的reads
    vector<ReadEnds> frags;              // 暂未找到配对的reads
    DPSet<DupInfo> pairDupIdx;           // pair的冗余read的索引
    MDSet<int64_t> pairOpticalDupIdx;    // optical冗余read的索引
    DPSet<DupInfo> fragDupIdx;           // frag的冗余read的索引
    DPSet<DupInfo> pairRepIdx;           // pair的dupset代表read的索引
    MDSet<int64_t> pairSingletonIdx;     // 某位置只有一对read的所有read pair个数
    UnpairedNameMap unpairedDic;         // 用来寻找pair end
    UnpairedPositionMap unpairedPosArr;  // 存放未匹配的ReadEnd对应位点的所有ReadEnd，为了避免重复存储
};

/*
 * 优先队列，用最小堆来实现对所有冗余索引的排序
 */
template <typename T>
struct PairArrIdIdx {
    int arrId = 0;
    uint64_t arrIdx = 0;
    T dupIdx = 0;
};

template <typename T>
struct IdxGreaterThan {
    bool operator()(const PairArrIdIdx<T> &a, const PairArrIdIdx<T> &b) { return a.dupIdx > b.dupIdx; }
};

template <typename T>
struct DupIdxQueue {
    // 将冗余索引和他对应的task vector对应起来

    // 由于是多个task来查找冗余，所以每次找到的冗余index都放在一个独立的vector中，vector之间可能有重叠，所以需要用一个最小堆来维护
    vector<vector<T>> *dupIdx2DArr;
    priority_queue<PairArrIdIdx<T>, vector<PairArrIdIdx<T>>, IdxGreaterThan<T>> minHeap;
    uint64_t popNum = 0;

    int Init(vector<vector<T>> *_dupIdx2DArr) {
        dupIdx2DArr = _dupIdx2DArr;
        if (dupIdx2DArr == nullptr) {
            return -1;
        }
        for (int i = 0; i < dupIdx2DArr->size(); ++i) {
            auto &v = (*dupIdx2DArr)[i];
            if (!v.empty()) {
                minHeap.push({i, 1, v[0]});
            }
        }
        return 0;
    }

    T Pop() {
        T ret = -1;
        if (!minHeap.empty()) {
            auto idx = minHeap.top();
            minHeap.pop();
            ++popNum;
            ret = idx.dupIdx;
            auto &v = (*dupIdx2DArr)[idx.arrId];
            if (v.size() > idx.arrIdx) {
                minHeap.push({idx.arrId, idx.arrIdx + 1, v[idx.arrIdx]});
            }
        }
        return ret;
    }

    uint64_t Size() {
        uint64_t len = 0;
        if (dupIdx2DArr != nullptr) {
            for (auto &v : *dupIdx2DArr) {
                len += v.size();
            }
        }
        return len - popNum;
    }

    uint64_t RealSize() {
        uint64_t len = 0;
        auto preTop = minHeap.top();
        DupInfo dupIdx = this->Pop();
        DupInfo nextDup = dupIdx;
        auto topIdx = minHeap.top();

        // ofstream ofs("n.dup"); ofstream ofs1("n-all.dup");

        while (dupIdx != -1) {
            ++len;
            while (true) {
                topIdx = minHeap.top();
                nextDup = this->Pop();
                if (nextDup != dupIdx) {
                    dupIdx = nextDup;
                    break;
                } else {
                    cout << "the same dup: " << dupIdx << '\t' << preTop.arrId << '\t' << preTop.arrIdx << '\t'
                         << preTop.dupIdx << '\t' << topIdx.arrId << '\t' << topIdx.arrIdx << '\t' << topIdx.dupIdx
                         << endl;
                }
            }
            
            // ofs << topIdx.dupIdx << endl; ofs1 << topIdx.arrId << '\t' << topIdx.arrIdx << '\t' << topIdx.dupIdx << endl;
            
            dupIdx = nextDup;
            preTop = topIdx;
        }
        // ofs.close(); ofs1.close();
        cout << "RealSize: " << len << endl;
        return len;
    }
};