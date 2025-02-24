#pragma once

#include <inttypes.h>
#include <spdlog/spdlog.h>
#include <util/yarn.h>

#include "md_types.h"

struct ReadData {
    vector<BamWrap *> bams;   // read step output
    int64_t bamStartIdx = 0;  // 每轮读入的bam数组，起始位置在全局bam中的索引
    int64_t taskSeq = 0;      // 任务序号
};

struct GenREData {
    vector<vector<ReadEnds>> pairsArr;       // 成对的reads
    vector<vector<ReadEnds>> fragsArr;       // 暂未找到配对的reads
    vector<UnpairedNameMap> unpairedDicArr;  // 用来寻找pair end
    void Init(int nThread) {
        for (int i = 0; i <= nThread; ++i) {  // 比线程多一个，主要是pairs多一个，其他没用
            pairsArr.push_back(vector<ReadEnds>());
            fragsArr.push_back(vector<ReadEnds>());
            unpairedDicArr.push_back(UnpairedNameMap());
        }
    }
    UnpairedNameMap unpairedDic;         // 代替sort step中一部分计算
    size_t byteSize() {
        size_t bytes = 0;
        for (auto &v : pairsArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &v : fragsArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &m : unpairedDicArr) bytes += m.size() * 100;
        bytes += unpairedDic.size() * 100;
        return bytes;
    }
};

struct SortMarkData {
    vector<ReadEnds> pairs;              // 成对的reads
    vector<ReadEnds> frags;              // 暂未找到配对的reads
    UnpairedNameMap unpairedDic;         // 用来寻找pair end
    size_t byteSize() {
        size_t bytes = 0;
        for (auto &r : pairs) bytes += sizeof(r);
        for (auto &r : frags) bytes += sizeof(r);
        bytes += unpairedDic.size() * 100;
        return bytes;
    }
};

struct SortData {
    volatile void *dataPtr;  // SortMarkData pointer
};

struct MarkDupData {
    int64_t taskSeq = 0;               // 任务序号
    DPSet<DupInfo> pairDupIdx;         // pair的冗余read的索引
    MDSet<int64_t> pairOpticalDupIdx;  // optical冗余read的索引
    DPSet<DupInfo> fragDupIdx;         // frag的冗余read的索引
    DPSet<DupInfo> pairRepIdx;         // pair的dupset代表read的索引
    CkeyReadEndsMap ckeyReadEndsMap;

    volatile void *dataPtr;  // SortMarkData pointer

    void clear() {
        fragDupIdx.clear();
        pairDupIdx.clear();
        pairOpticalDupIdx.clear();
        pairRepIdx.clear();
        ckeyReadEndsMap.clear();
    }

    size_t byteSize() {
        size_t bytes = 0;
        bytes += pairDupIdx.size() * 100;
        bytes += pairOpticalDupIdx.size() * 100;
        bytes += fragDupIdx.size() * 100;
        bytes += pairRepIdx.size() * 100;
        return bytes;
    }
};

struct DupResult {
    vector<vector<DupInfo>> dupIdxArr;
    vector<vector<int64_t>> opticalDupIdxArr;
    vector<vector<DupInfo>> repIdxArr;
    size_t byteSize() {
        size_t bytes = 0;
        size_t tmp = 0;
        for (auto &v : dupIdxArr)
            for (auto &r : v) tmp += sizeof(r);
        spdlog::info("dupIdxArr size : {} GB", tmp / 1024.0 / 1024 / 1024);
        bytes += tmp;
        tmp = 0;
        for (auto &v : opticalDupIdxArr)
            for (auto &r : v) tmp += sizeof(r);
        spdlog::info("opticalDupIdxArr size : {} GB", tmp / 1024.0 / 1024 / 1024);
        bytes += tmp;
        tmp = 0;
        for (auto &v : repIdxArr)
            for (auto &r : v) tmp += sizeof(r);
        spdlog::info("repIdxArr size : {} GB", tmp / 1024.0 / 1024 / 1024);
        bytes += tmp;
        spdlog::info("result size : {} GB", bytes / 1024.0 / 1024 / 1024);

        return bytes;
    }
};

struct IntersectData {
    UnpairedNameMap unpairedDic;  // 用来寻找pair end
    CkeyReadEndsMap ckeyReadEndsMap;

    // 每个task对应一个vector
    vector<vector<DupInfo>> &dupIdxArr;
    vector<vector<int64_t>> &opticalDupIdxArr;
    vector<vector<DupInfo>> &repIdxArr;

    IntersectData(DupResult *resPtr)
        : dupIdxArr(resPtr->dupIdxArr), opticalDupIdxArr(resPtr->opticalDupIdxArr), repIdxArr(resPtr->repIdxArr) {}

    size_t byteSize() {
        size_t bytes = 0;
        bytes += unpairedDic.size() * 100;
        for (auto &v : dupIdxArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &v : opticalDupIdxArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &v : repIdxArr)
            for (auto &r : v) bytes += sizeof(r);
        spdlog::info("result size : {}", bytes);

        return bytes;
    }
};

// 记录流水线状态，task的序号，以及某阶段是否结束
struct PipelineArg {
    static const int GENBUFNUM = 2;
    static const int SORTBUFNUM = 2;
    static const int MARKBUFNUM = 4;
    uint64_t readOrder = 0;
    uint64_t genREOrder = 0;
    uint64_t sortOrder = 0;
    uint64_t markDupOrder = 0;
    uint64_t intersectOrder = 0;
    int numThread = 0;

    volatile int readFinish = 0;
    volatile int genREFinish = 0;
    volatile int sortFinish = 0;
    volatile int markDupFinish = 0;

    yarn::lock_t *readSig;
    yarn::lock_t *genRESig;
    yarn::lock_t *sortSig;
    yarn::lock_t *markDupSig;

    PipelineArg(DupResult *resPtr) : intersectData(resPtr) {
        readSig = yarn::NEW_LOCK(0);   // 最大值1, 双buffer在bambuf中实现了，对调用透明
        genRESig = yarn::NEW_LOCK(0);  // 最大值2, 双buffer
        sortSig = yarn::NEW_LOCK(0);
        markDupSig = yarn::NEW_LOCK(0);
        for (int i = 0; i < SORTBUFNUM; ++i) {
            sortData[i].dataPtr = &sortMarkData[i];
        }
        for (int i = 0; i < MARKBUFNUM; ++i) {
            markDupData[i].dataPtr = &sortMarkData[i + SORTBUFNUM];
        }
    }

    SortMarkData sortMarkData[SORTBUFNUM + MARKBUFNUM];

    // for step-1 read
    ReadData readData;
    // for step-2 generate readends
    GenREData genREData[GENBUFNUM];
    // for step-3 sort each thread frags and pairs
    SortData sortData[SORTBUFNUM];
    // for step-4 mark duplicate
    MarkDupData markDupData[MARKBUFNUM];
    // for step-5 deal with intersect data
    IntersectData intersectData;

    size_t byteSize() {
        size_t bytes = 0;

        size_t tmp = 0;
        for (int i = 0; i < SORTBUFNUM + MARKBUFNUM; ++i) tmp += sortMarkData[i].byteSize();
        bytes += tmp;
        spdlog::info("sortMarkData size : {}", tmp);
        for (int i = 0; i < GENBUFNUM; ++i) tmp += genREData[i].byteSize();
        bytes += tmp;
        spdlog::info("genREData size : {}", tmp);
        for (int i = 0; i < MARKBUFNUM; ++i) tmp += markDupData[i].byteSize();
        bytes += tmp;
        spdlog::info("markDupData size : {}", tmp);
        tmp += intersectData.byteSize();
        bytes += tmp;
        spdlog::info("intersectData size : {}", tmp);

        return bytes;
    }
};

struct REArrIdIdx {
    int arrId = 0;        // 数组索引
    uint64_t arrIdx = 0;  // 数组内部当前索引
    const ReadEnds *pE = nullptr;
};

struct REGreaterThan {
    bool operator()(const REArrIdIdx &a, const REArrIdIdx &b) { return *b.pE < *a.pE; }
};

struct ReadEndsHeap {
    // 将冗余索引和他对应的task vector对应起来
    vector<vector<ReadEnds>> *arr2d;
    priority_queue<REArrIdIdx, vector<REArrIdIdx>, REGreaterThan> minHeap;
    uint64_t popNum = 0;

    int Init(vector<vector<ReadEnds>> *_arr2d) {
        arr2d = _arr2d;
        if (arr2d == nullptr) {
            return -1;
        }
        for (int i = 0; i < arr2d->size(); ++i) {
            auto &v = (*arr2d)[i];
            if (!v.empty()) {
                minHeap.push({i, 1, &v[0]});
            }
        }
        return 0;
    }

    const ReadEnds *Pop() {
        const ReadEnds *ret = nullptr;
        if (!minHeap.empty()) {
            auto minVal = minHeap.top();
            minHeap.pop();
            ++popNum;
            ret = minVal.pE;
            auto &v = (*arr2d)[minVal.arrId];
            if (v.size() > minVal.arrIdx) {
                minHeap.push({minVal.arrId, minVal.arrIdx + 1, &v[minVal.arrIdx]});
            }
        }
        return ret;
    }

    uint64_t Size() {
        uint64_t len = 0;
        if (arr2d != nullptr) {
            for (auto &v : *arr2d) {
                len += v.size();
            }
        }
        return len - popNum;
    }
};

// 并行运行mark duplicate
void NewPipeMarkDups();