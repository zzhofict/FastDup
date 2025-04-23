#pragma once

#include <inttypes.h>
#include <spdlog/spdlog.h>
#include <util/yarn.h>

#include "md_types.h"

struct ReadData {
    vector<BamWrap *> bams;   // read step output
    int64_t bamStartIdx = 0;  // bam，bam
    int64_t taskSeq = 0;      // 
};

struct GenREData {
    vector<vector<ReadEnds>> pairsArr;       // reads
    vector<vector<ReadEnds>> fragsArr;       // reads
    vector<UnpairedNameMap> unpairedDicArr;  // pair end
    void Init(int nThread) {
        for (int i = 0; i <= nThread; ++i) {  // ，pairs，
            pairsArr.push_back(vector<ReadEnds>());
            fragsArr.push_back(vector<ReadEnds>());
            unpairedDicArr.push_back(UnpairedNameMap());
        }
    }
    UnpairedNameMap unpairedDic;         // sort step
    size_t byteSize() {
        size_t bytes = 0;
        for (auto &v : pairsArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &v : fragsArr)
            for (auto &r : v) bytes += sizeof(r);
        spdlog::info("genRE readends size : {} GB", bytes / 1024.0 / 1024 / 1024);
        for (auto &m : unpairedDicArr) bytes += m.size() * sizeof(*m.begin());
        bytes += unpairedDic.size() * sizeof(*unpairedDic.begin());
        return bytes;
    }
};

struct SortMarkData {
    vector<ReadEnds> pairs;              // reads
    vector<ReadEnds> frags;              // reads
    UnpairedNameMap unpairedDic;         // pair end
    size_t byteSize() {
        size_t bytes = 0;
        for (auto &r : pairs) bytes += sizeof(r);
        for (auto &r : frags) bytes += sizeof(r);
        spdlog::info("sortmark readends size : {} GB", bytes / 1024.0 / 1024 / 1024);
        bytes += unpairedDic.bucket_count() * sizeof(*unpairedDic.begin());
        return bytes;
    }
};

struct SortData {
    volatile void *dataPtr;  // SortMarkData pointer
};

struct MarkDupData {
    int64_t taskSeq = 0;               // 
    DPSet<DupInfo> pairDupIdx;         // pairread
    MDSet<int64_t> pairOpticalDupIdx;  // opticalread
    DPSet<DupInfo> fragDupIdx;         // fragread
    DPSet<DupInfo> pairRepIdx;         // pairdupsetread
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
    UnpairedNameMap unpairedDic;  // pair end
    CkeyReadEndsMap ckeyReadEndsMap;
    int64_t rightPos = 0;

    // taskvector
    vector<vector<DupInfo>> &dupIdxArr;
    vector<vector<int64_t>> &opticalDupIdxArr;
    vector<vector<DupInfo>> &repIdxArr;

    IntersectData(DupResult *resPtr)
        : dupIdxArr(resPtr->dupIdxArr), opticalDupIdxArr(resPtr->opticalDupIdxArr), repIdxArr(resPtr->repIdxArr) {}

    size_t byteSize() {
        size_t bytes = 0;
        bytes += unpairedDic.size() * sizeof(*unpairedDic.begin());
        for (auto &v : dupIdxArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &v : opticalDupIdxArr)
            for (auto &r : v) bytes += sizeof(r);
        for (auto &v : repIdxArr)
            for (auto &r : v) bytes += sizeof(r);
        spdlog::info("result size : {} GB", bytes / 1024.0 / 1024 / 1024);

        return bytes;
    }
};

// ，task，
struct PipelineArg {
    static const int GENBUFNUM = 2;
    static const int SORTBUFNUM = 2;
    static const int MARKBUFNUM = 3;
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
        readSig = yarn::NEW_LOCK(0);   // 1, bufferbambuf，
        genRESig = yarn::NEW_LOCK(0);  // 2, buffer
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
        spdlog::info("sortMarkData size : {} GB", tmp / 1024.0 / 1024 / 1024);
        for (int i = 0; i < GENBUFNUM; ++i) tmp += genREData[i].byteSize();
        bytes += tmp;
        spdlog::info("genREData size : {} GB", tmp / 1024.0 / 1024 / 1024);
        for (int i = 0; i < MARKBUFNUM; ++i) tmp += markDupData[i].byteSize();
        bytes += tmp;
        spdlog::info("markDupData size : {} GB", tmp / 1024.0 / 1024 / 1024);
        tmp += intersectData.byteSize();
        bytes += tmp;
        spdlog::info("intersectData size : {} GB", tmp / 1024.0 / 1024 / 1024);

        return bytes;
    }
};

struct REArrIdIdx {
    int arrId = 0;        // 
    uint64_t arrIdx = 0;  // 
    const ReadEnds *pE = nullptr;
};

struct REFragGreaterThan {
    bool operator()(const REArrIdIdx &a, const REArrIdIdx &b) { return ReadEnds::FragLittleThan(*b.pE, *a.pE); }
};

struct REPairGreaterThan {
    bool operator()(const REArrIdIdx &a, const REArrIdIdx &b) { return ReadEnds::PairLittleThan(*b.pE, *a.pE);
    }
};

template <typename CMP>
struct ReadEndsHeap {
    // task vector
    vector<vector<ReadEnds>> *arr2d;
    priority_queue<REArrIdIdx, vector<REArrIdIdx>, CMP> minHeap;
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

// mark duplicate
void PipelineMarkDups();