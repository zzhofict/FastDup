#include "pipeline_md.h"

#include <klib/kthread.h>
#include <pthread.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "dup_metrics.h"
#include "md_args.h"
#include "md_funcs.h"
#include "read_ends.h"
#include "read_name_parser.h"
#include "util/bam_buf.h"
#include "util/profiling.h"
#include "util/yarn.h"

using std::cout;
using std::vector;

namespace nsgv {

extern MarkDupsArg gMdArg;                   // 用来测试性能
extern samFile *gInBamFp;                    // 输入的bam文件
extern sam_hdr_t *gInBamHeader;              // 输入的bam文件头信息
extern DuplicationMetrics gMetrics;          // 统计信息
extern vector<ReadNameParser> gNameParsers;
extern PipelineArg gPipe;

};  // namespace nsgv

/* 排序 */
static inline void sortReadEndsArr(vector<ReadEnds> &arr) {
    size_t blockSize = 64 * 1024;
    if (arr.size() < blockSize) {
        std::sort(arr.begin(), arr.end());
        return;
    }
    size_t blockNum = (arr.size() + blockSize - 1) / blockSize;
    size_t crossNum = 1024;
    size_t start, end, i, left, right;
    std::sort(arr.begin(), arr.begin() + blockSize);
    for (i = 1; i < blockNum; ++i) {
        start = i * blockSize;
        end = min(start + blockSize, arr.size());
        std::sort(arr.begin() + start, arr.begin() + end);
        left = crossNum;
        while (!(arr[start - left] < arr[start])) {
            left = left << 1;
            if (left >= blockSize) {
                std::sort(arr.begin(), arr.end());  // 退化到普通排序
                return;
            }
        }
        right = min(crossNum, end - start - 1);

        while (!(arr[start - 1] < arr[start + right])) {
            right = min(right << 1, end - start - 1);
            if (right == end - start - 1)
                break;
        }
        std::sort(arr.begin() + start - left, arr.begin() + start + right);
    }
}

/* 处理一组pairend的readends，标记冗余, 这个函数需要串行运行，因为需要做一些统计*/
static void markDupsForPairs(vector<const ReadEnds *> &vpRe,
                             DPSet<DupInfo> *dupIdx,
                             MDSet<int64_t> *opticalDupIdx,
                             DPSet<DupInfo> *repIdx,
                             MDSet<int64_t> *notDupIdx = nullptr,
                             MDSet<int64_t> *notOpticalDupIdx = nullptr,
                             MDSet<int64_t> *notRepIdx = nullptr) {
    if (vpRe.size() < 2) {
        return;
    }

    int maxScore = 0;
    const ReadEnds *pBest = nullptr;
    /** All read ends should have orientation FF, FR, RF, or RR **/
    for (auto pe : vpRe) {  // 找分数最高的readend
        if (pe->score > maxScore || pBest == nullptr) {
            maxScore = pe->score;
            pBest = pe;
        }
    }

    if (notDupIdx != nullptr) {
        notDupIdx->insert(pBest->read1IndexInFile);
        notDupIdx->insert(pBest->read2IndexInFile);
    }

    if (nsgv::gMdArg.CHECK_OPTICAL_DUP) {  // 检查光学冗余
        // trackOpticalDuplicates
        vector<const ReadEnds *> prevOpticalRe;
        if (notOpticalDupIdx != nullptr) {
            for (auto pe : vpRe) {
                if (pe->isOpticalDuplicate) {
                    prevOpticalRe.push_back(pe);
                }
            }
        }
        trackOpticalDuplicates(vpRe, pBest);
        // 由于重叠问题，可能会更新数据
        if (notOpticalDupIdx != nullptr) {
            for (auto pe : prevOpticalRe) {
                if (!pe->isOpticalDuplicate) {
                    notOpticalDupIdx->insert(pe->read1IndexInFile);
                    notOpticalDupIdx->insert(pe->read2IndexInFile);
                }
            }
        }
    }
    for (auto pe : vpRe) {  // 对非best read标记冗余
        if (pe != pBest) {  // 非best
            dupIdx->insert(DupInfo(pe->read1IndexInFile, pBest->read1IndexInFile, (int16_t)vpRe.size()));  // 添加read1
            if (pe->read2IndexInFile != pe->read1IndexInFile)
                dupIdx->insert(DupInfo(pe->read2IndexInFile, pBest->read1IndexInFile, (int16_t)vpRe.size()));  // read2,
            //  注意这里代表是read1的索引
            //  检查是否optical dup
            if (pe->isOpticalDuplicate && opticalDupIdx != nullptr) {
                opticalDupIdx->insert(pe->read1IndexInFile);
                if (pe->read2IndexInFile != pe->read1IndexInFile)
                    opticalDupIdx->insert(pe->read2IndexInFile);
            }
        }
    }
    // 在输出的bam文件中添加tag, best作为dupset的代表
    if (nsgv::gMdArg.TAG_DUPLICATE_SET_MEMBERS) {
        repIdx->insert(DupInfo(pBest->read1IndexInFile, pBest->read1IndexInFile, (int16_t)vpRe.size()));
        repIdx->insert(DupInfo(pBest->read2IndexInFile, pBest->read1IndexInFile, (int16_t)vpRe.size()));
        if (notRepIdx != nullptr) {
            for (auto pe : vpRe) {
                if (pe != pBest) {
                    notRepIdx->insert(pe->read1IndexInFile);
                    if (pe->read2IndexInFile != pe->read1IndexInFile)
                        notRepIdx->insert(pe->read2IndexInFile);
                }
            }
        }
    }
}

/* 处理一组非paired的readends，标记冗余 */
static void markDupsForFrags(vector<const ReadEnds *> &vpRe, bool containsPairs, DPSet<DupInfo> *dupIdx,
                             MDSet<int64_t> *notDupIdx = nullptr) {
    if (containsPairs) {
        for (auto pe : vpRe) {
            if (!pe->IsPaired()) {
                dupIdx->insert(pe->read1IndexInFile);
            }
        }
    } else {
        int maxScore = 0;
        const ReadEnds *pBest = nullptr;
        for (auto pe : vpRe) {
            if (pe->score > maxScore || pBest == nullptr) {
                maxScore = pe->score;
                pBest = pe;
            }
        }
        if (notDupIdx != nullptr) {
            notDupIdx->insert(pBest->read1IndexInFile);
        }
        for (auto pe : vpRe) {
            if (pe != pBest) {
                dupIdx->insert(pe->read1IndexInFile);
            }
        }
    }
}

/* 找到与readend pos相等的所有readend */
static void getEqualRE(const ReadEnds &re, vector<ReadEnds> &src, vector<ReadEnds> *dst) {
    auto range = std::equal_range(src.begin(), src.end(), re, ReadEnds::PairsLittleThan);  // 只比对位点
    dst->insert(dst->end(), range.first, range.second);
}

#define LOWER_BOUND(tid, nthread, ndata) ((tid) * (ndata) / (nthread))
#define UPPER_BOUND(tid, nthread, ndata) ((tid + 1) * (ndata) / (nthread))

/* 处理pairs */
static void processPairs(vector<ReadEnds> &readEnds, DPSet<DupInfo> *dupIdx, MDSet<int64_t> *opticalDupIdx,
                         DPSet<DupInfo> *repIdx, MDSet<int64_t> *notDupIdx = nullptr,
                         MDSet<int64_t> *notOpticalDupIdx = nullptr, MDSet<int64_t> *notRepIdx = nullptr) {
    // return;
    vector<const ReadEnds *> vpCache;  // 有可能是冗余的reads
    const ReadEnds *pReadEnd = nullptr;
    for (auto &re : readEnds) {
        if (pReadEnd != nullptr && ReadEnds::AreComparableForDuplicates(*pReadEnd, re, true))  // 跟前一个一样
            vpCache.push_back(&re);                                                            // 处理一个潜在的冗余组
        else {
            markDupsForPairs(vpCache, dupIdx, opticalDupIdx, repIdx, notDupIdx, notOpticalDupIdx,
                             notRepIdx);  // 不一样
            vpCache.clear();
            vpCache.push_back(&re);
            pReadEnd = &re;
        }
    }
    markDupsForPairs(vpCache, dupIdx, opticalDupIdx, repIdx, notDupIdx, notOpticalDupIdx, notRepIdx);
}

/* 处理frags */
static void processFrags(vector<ReadEnds> &readEnds, DPSet<DupInfo> *dupIdx, MDSet<int64_t> *notDupIdx = nullptr) {
    bool containsPairs = false;
    bool containsFrags = false;
    vector<const ReadEnds *> vpCache;  // 有可能是冗余的reads
    const ReadEnds *pReadEnd = nullptr;
    for (auto &re : readEnds) {
        if (pReadEnd != nullptr && ReadEnds::AreComparableForDuplicates(*pReadEnd, re, false)) {
            vpCache.push_back(&re);
            containsPairs = containsPairs || re.IsPaired();
            containsFrags = containsFrags || !re.IsPaired();
        } else {
            if (vpCache.size() > 1 && containsFrags) {
                markDupsForFrags(vpCache, containsPairs, dupIdx, notDupIdx);
            }
            vpCache.clear();
            vpCache.push_back(&re);
            pReadEnd = &re;
            containsPairs = re.IsPaired();
            containsFrags = !re.IsPaired();
        }
    }
    if (vpCache.size() > 1 && containsFrags) {
        markDupsForFrags(vpCache, containsPairs, dupIdx, notDupIdx);
    }
}

/* 单线程markdup (第二步)*/
static void markdups(MarkDupDataArg *arg) {
    auto &p = *arg;
    p.fragDupIdx.clear();
    p.pairDupIdx.clear();
    p.pairOpticalDupIdx.clear();
    p.pairRepIdx.clear();
    p.pairSingletonIdx.clear();
    /* generateDuplicateIndexes，计算冗余read在所有read中的位置索引 */
    // 先处理 pair
    processPairs(p.pairs, &p.pairDupIdx, &p.pairOpticalDupIdx, &p.pairRepIdx, &p.pairSingletonIdx);
    // 再处理frag
    processFrags(p.frags, &p.fragDupIdx);
}

/* 获取交叉部分的数据 */
static inline void getIntersectData(vector<ReadEnds> &leftArr, vector<ReadEnds> &rightArr, vector<ReadEnds> *dst,
                                    bool isPairCmp = false) {
    if (leftArr.empty() || rightArr.empty()) {
        // cout << "bad size: " << leftArr.size() << '\t' << rightArr.size() << endl;
        return;
    }
    const size_t leftEndIdx = leftArr.size() - 1;
    const size_t rightStartIdx = 0;
    size_t leftSpan = 0;
    size_t rightSpan = 0;

    while (!ReadEnds::ReadLittleThan(leftArr[leftEndIdx - leftSpan], rightArr[rightStartIdx], isPairCmp)) {
        leftSpan += 1;
        if (leftSpan > leftEndIdx) {  // 上一个的范围被下一个全部包围了（可能会有bug，上上个也与下一个有交集呢？）
            leftSpan = leftArr.size() - 1;
            break;
        }
    }

    while (!ReadEnds::ReadLittleThan(leftArr[leftEndIdx], rightArr[rightSpan], isPairCmp)) {
        rightSpan += 1;
        if (rightSpan == rightArr.size() - 1)  // 同上，可能会有bug
            break;
    }
    dst->insert(dst->end(), leftArr.end() - leftSpan, leftArr.end());
    dst->insert(dst->end(), rightArr.begin(), rightArr.begin() + rightSpan);
    std::sort(dst->begin(), dst->end());
}

/* 将frags重叠部分的dup idx放进数据中 */
static inline void refreshFragDupIdx(DPSet<DupInfo> &dupIdx, MDSet<int64_t> &notDupIdx, MarkDupDataArg *lastArg,
                                     MarkDupDataArg *curArg) {
    auto &lp = *lastArg;
    auto &p = *curArg;
    for (auto idx : dupIdx) {
        lp.fragDupIdx.insert(idx);
        p.fragDupIdx.erase(idx);
    }
    for (auto idx : notDupIdx) {
        lp.fragDupIdx.erase(idx);
        p.fragDupIdx.erase(idx);
    }
}

// 用来分别处理dup和optical dup
static void refeshTaskDupInfo(DPSet<DupInfo> &dupIdx, MDSet<int64_t> &opticalDupIdx, DPSet<DupInfo> &repIdx,
                              MDSet<int64_t> &notDupIdx, MDSet<int64_t> &notOpticalDupIdx, MDSet<int64_t> &notRepIdx,
                              DPSet<DupInfo> &latterDupIdx, MDSet<int64_t> &latterOpticalDupIdx,
                              DPSet<DupInfo> &latterRepIdx, MDSet<int64_t> &latterNotDupIdx,
                              MDSet<int64_t> &latterNotOpticalDupIdx, MDSet<int64_t> &latterNotRepIdx) {
    for (auto idx : dupIdx) latterDupIdx.insert(idx);
    for (auto idx : opticalDupIdx) latterOpticalDupIdx.insert(idx);
    for (auto idx : repIdx) latterRepIdx.insert(idx);
    for (auto idx : notDupIdx) latterNotDupIdx.insert(idx);
    for (auto idx : notOpticalDupIdx) latterNotOpticalDupIdx.insert(idx);
    for (auto idx : notRepIdx) latterNotRepIdx.insert(idx);
}

/* 最后合并数据并排序 */
template <typename DupContainer, typename T>
static void refeshFinalTaskDupInfo(DupContainer &dupIdx, MDSet<int64_t> &notDupIdx, vector<T> &dupArr,
                                   vector<T> &cacheDupIdx, vector<T> &midArr) {
    midArr.resize(0);
    cacheDupIdx.resize(0);
    cacheDupIdx.insert(cacheDupIdx.end(), dupIdx.begin(), dupIdx.end());
    std::sort(cacheDupIdx.begin(), cacheDupIdx.end());

    auto ai = dupArr.begin();
    auto ae = dupArr.end();
    auto bi = cacheDupIdx.begin();
    auto be = cacheDupIdx.end();

    T val = 0;
    while (ai != ae && bi != be) {
        if (*ai < *bi) {
            val = *ai++;
        } else if (*bi < *ai) {
            val = *bi++;
        } else {
            val = *bi++;  // 相等的时候取后面的作为结果
            ai++;
        }
        if (notDupIdx.find(val) == notDupIdx.end()) {
            midArr.push_back(val);
        }
    }
    while (ai != ae) {
        val = *ai++;
        if (notDupIdx.find(val) == notDupIdx.end()) {
            midArr.push_back(val);
        }
    }
    while (bi != be) {
        val = *bi++;
        if (notDupIdx.find(val) == notDupIdx.end()) {
            midArr.push_back(val);
        }
    }
    dupArr = midArr;
}

// for step 2 generate read ends

// multi-thread generate read ends
static void mtGenerateReadEnds(void *data, long idx, int tid) {
    auto &p = *(PipelineArg *)data;
    auto &rnParser = nsgv::gNameParsers[idx];
    int nThread = p.numThread;
    auto &bams = p.readData.bams;
    int64_t bamStartIdx = p.readData.bamStartIdx;
    int64_t taskSeq = p.readData.taskSeq;
    GenREData &genREData = p.genREData[p.genREOrder % p.GENBUFNUM];
    auto &pairs = genREData.pairsArr[idx];
    auto &frags = genREData.fragsArr[idx];
    auto &unpairedDic = genREData.unpairedDicArr[idx];

    pairs.clear();
    frags.clear();
    unpairedDic.clear();

    size_t start_id = LOWER_BOUND(idx, nThread, bams.size());
    size_t end_id = UPPER_BOUND(idx, nThread, bams.size());
    for (size_t i = start_id; i < end_id; ++i) {  // 循环处理每个read
        BamWrap *bw = bams[i];
        const int64_t bamIdx = bamStartIdx + i;
        if (bw->GetReadUnmappedFlag()) {
            if (bw->b->core.tid == -1)
                // When we hit the unmapped reads with no coordinate, no reason to continue (only in coordinate sort).
                break;
        } else if (!bw->IsSecondaryOrSupplementary()) {  // 是主要比对
            ReadEnds fragEnd;
            buildReadEnds(*bw, bamIdx, rnParser, &fragEnd);
            frags.push_back(fragEnd);                                     // 添加进frag集合
            if (bw->GetReadPairedFlag() && !bw->GetMateUnmappedFlag()) {  // 是pairend而且互补的read也比对上了
                string key = bw->query_name();
                if (unpairedDic.find(key) == unpairedDic.end()) {
                    unpairedDic[key] = {taskSeq, fragEnd};
                } else {  // 找到了pairend
                    auto &pairedEnds = unpairedDic.at(key).unpairedRE;
                    modifyPairedEnds(fragEnd, &pairedEnds);
                    pairs.push_back(pairedEnds);
                    unpairedDic.erase(key);  // 删除找到的pairend
                }
            }
        }
    }
    // sortReadEndsArr(frags);
    sort(frags.begin(), frags.end());
    sort(pairs.begin(), pairs.end());
}

static void doGenRE(PipelineArg &pipeArg) {
    // return;
    GenREData &genREData = pipeArg.genREData[pipeArg.genREOrder % pipeArg.GENBUFNUM];
    ReadData &readData = pipeArg.readData;

    // generate read ends
    const int numThread = pipeArg.numThread;

    kt_for(numThread, mtGenerateReadEnds, &pipeArg, numThread);
    // 用来在genRE阶段计算一些Sort阶段应该计算的数据，保持每个step用时更平衡
    // 轮询每个线程中未找到匹配的read，找到匹配的
    genREData.unpairedDic.clear();
    genREData.unpairedPosArr.clear();
    vector<ReadEnds> &pairs = genREData.pairsArr[numThread];
    pairs.clear();
    for (int i = 0; i < numThread; ++i) {
        for (auto &val : genREData.unpairedDicArr[i]) {
            const string &key = val.first;
            const ReadEnds &fragEnd = val.second.unpairedRE;
            if (genREData.unpairedDic.find(key) == genREData.unpairedDic.end()) {
                genREData.unpairedDic[key] = {readData.taskSeq, fragEnd};
            } else {  // 找到了pairend
                auto &pairedEnds = genREData.unpairedDic.at(key).unpairedRE;
                modifyPairedEnds(fragEnd, &pairedEnds);
                pairs.push_back(pairedEnds);
                genREData.unpairedDic.erase(key);  // 删除找到的pairend
            }
        }
    }
    sort(pairs.begin(), pairs.end());

    // create unpaired info
    for (auto &e : genREData.unpairedDic) {
        auto posKey = e.second.unpairedRE.posKey;
        auto &unpairArrInfo = genREData.unpairedPosArr[posKey];
        unpairArrInfo.unpairedNum++;
        unpairArrInfo.taskSeq = readData.taskSeq;
        unpairArrInfo.readNameSet.insert(e.first);
    }
}
// end for step 2 generate read ends

// for step-3 sort
static void doSort(PipelineArg &pipeArg) {
    // return;
    GenREData &genREData = pipeArg.genREData[pipeArg.sortOrder % pipeArg.GENBUFNUM];
    SortData &sortData = pipeArg.sortData[pipeArg.sortOrder % pipeArg.SORTBUFNUM];
    SortMarkData &smd = *(SortMarkData *)sortData.dataPtr;

    smd.unpairedDic = genREData.unpairedDic;
    smd.unpairedPosArr = genREData.unpairedPosArr;

    smd.pairs.clear();
    smd.frags.clear();
    const ReadEnds *pRE;
    ReadEndsHeap pairsHeap, fragsHeap;
    pairsHeap.Init(&genREData.pairsArr);
    while ((pRE = pairsHeap.Pop()) != nullptr) {
        smd.pairs.push_back(*pRE);
    }
    fragsHeap.Init(&genREData.fragsArr);
    while ((pRE = fragsHeap.Pop()) != nullptr) {
        smd.frags.push_back(*pRE);
    }
}
// for step-4 sort
static void doMarkDup(PipelineArg &pipeArg) {
    // return;
    MarkDupData &mdData = pipeArg.markDupData[pipeArg.markDupOrder % pipeArg.MARKBUFNUM];
    SortData &sortData = pipeArg.sortData[pipeArg.markDupOrder % pipeArg.SORTBUFNUM];
    mdData.taskSeq = pipeArg.markDupOrder;
    mdData.fragDupIdx.clear();
    mdData.pairDupIdx.clear();
    mdData.pairOpticalDupIdx.clear();
    mdData.pairRepIdx.clear();

    auto tmpPtr = mdData.dataPtr;
    mdData.dataPtr = sortData.dataPtr;
    sortData.dataPtr = tmpPtr;

    SortMarkData &smd = *(SortMarkData *)mdData.dataPtr;
    //  先处理 pair
    processPairs(smd.pairs, &mdData.pairDupIdx, &mdData.pairOpticalDupIdx, &mdData.pairRepIdx);
    // 再处理frag
    processFrags(smd.frags, &mdData.fragDupIdx);
}

template <typename T>
static void refreshDupIdx(T &srcArr, T &insArr, T &delArr) {
    for (auto dup : srcArr) {
        insArr.insert(dup);
        delArr.erase(dup);
    }
}
template <typename T1, typename T2>
static void refreshNotDupIdx(T1 &srcArr, T2 &delArr1, T2 &delArr2) {
    for (auto dup : srcArr) {
        delArr1.erase(dup);
        delArr2.erase(dup);
    }
}

static void refreshMarkDupData(DPSet<DupInfo> &dupIdx, MDSet<int64_t> &opticalDupIdx, DPSet<DupInfo> &repIdx,
                               MDSet<int64_t> &notDupIdx, MDSet<int64_t> &notOpticalDupIdx, MDSet<int64_t> &notRepIdx,
                               MarkDupData &lp, MarkDupData &p) {
    refreshDupIdx(dupIdx, lp.pairDupIdx, p.pairDupIdx);
    refreshDupIdx(opticalDupIdx, lp.pairOpticalDupIdx, p.pairOpticalDupIdx);
    refreshDupIdx(repIdx, lp.pairRepIdx, p.pairRepIdx);

    refreshNotDupIdx(notDupIdx, lp.pairDupIdx, p.pairDupIdx);
    refreshNotDupIdx(notOpticalDupIdx, lp.pairOpticalDupIdx, p.pairOpticalDupIdx);
    refreshNotDupIdx(notRepIdx, lp.pairRepIdx, p.pairRepIdx);
}

// for step-5 sort
static void doIntersect(PipelineArg &pipeArg) {
    // return;
    IntersectData &g = pipeArg.intersectData;
    MarkDupData &lp = pipeArg.markDupData[(pipeArg.intersectOrder - 1) % pipeArg.MARKBUFNUM];
    MarkDupData &p = pipeArg.markDupData[pipeArg.intersectOrder % pipeArg.MARKBUFNUM];
    SortMarkData &lpSM = *(SortMarkData *)lp.dataPtr;
    SortMarkData &pSM = *(SortMarkData *)p.dataPtr;

    vector<ReadEnds> reArr;
    DPSet<DupInfo> dupIdx;
    MDSet<int64_t> opticalDupIdx;
    DPSet<DupInfo> repIdx;
    MDSet<int64_t> notOpticalDupIdx;
    MDSet<int64_t> notDupIdx;
    MDSet<int64_t> notRepIdx;

    // 先处理重叠的frags
    getIntersectData(lpSM.frags, pSM.frags, &reArr);
    processFrags(reArr, &dupIdx, &notDupIdx);
    refreshDupIdx(dupIdx, lp.fragDupIdx, p.fragDupIdx);
    refreshNotDupIdx(notDupIdx, lp.fragDupIdx, p.fragDupIdx);

    // 再处理重叠的pairs
    reArr.clear();
    dupIdx.clear();
    notDupIdx.clear();
    getIntersectData(lpSM.pairs, pSM.pairs, &reArr, true);

    processPairs(reArr, &dupIdx, &opticalDupIdx, &repIdx, &notDupIdx, &notOpticalDupIdx, &notRepIdx);
    refreshMarkDupData(dupIdx, opticalDupIdx, repIdx, notDupIdx, notOpticalDupIdx, notRepIdx, lp, p);

    //  处理之前未匹配的部分
    map<CalcKey, int64_t> recalcPos;
    CalcSet<CalcKey> alreadyAdd;  // 与该位点相同的pair都添加到数组里了
    MDSet<int64_t> addToGlobal;
    int64_t prevLastPos = 0, nextFirstPos = 0;
    if (lpSM.frags.size() > 0)
        prevLastPos = lpSM.frags.back().posKey;
    if (pSM.frags.size() > 0)
        nextFirstPos = pSM.frags[0].posKey;

    for (auto &prevUnpair : lpSM.unpairedDic) {  // 遍历上一个任务中的每个未匹配的read
        auto &readName = prevUnpair.first;
        auto &prevPosInfo = prevUnpair.second;
        auto prevFragEnd = prevPosInfo.unpairedRE;                      // 未匹配的read end
        if (pSM.unpairedDic.find(readName) != pSM.unpairedDic.end()) {  // 在当前这个任务里找到了这个未匹配的read
            auto &nextPosInfo = pSM.unpairedDic[readName];
            auto &nextFragEnd = nextPosInfo.unpairedRE;
            int64_t prevPosKey = prevFragEnd.posKey;
            modifyPairedEnds(nextFragEnd, &prevFragEnd);  // 在某些clip情况下，poskey可能是后面的read
            int64_t nextPosKey = max(prevPosKey, nextFragEnd.posKey);
            CalcKey ck = {prevPosKey, nextPosKey};
            UnpairedPosInfo *prevUnpairInfoP = nullptr;
            UnpairedPosInfo *nextUnpairInfoP = nullptr;

            if (lpSM.unpairedPosArr.find(prevPosKey) != lpSM.unpairedPosArr.end())
                prevUnpairInfoP = &lpSM.unpairedPosArr[prevPosKey];
            if (pSM.unpairedPosArr.find(prevPosKey) != pSM.unpairedPosArr.end())
                nextUnpairInfoP = &pSM.unpairedPosArr[prevPosKey];
            // pos分为两种情况，根据poskey(pair中两个read分别的pos)的位置确定
            //  1.
            //  prevpos在交叉部分之前，nextpos在交叉部分之后，这种情况不需要获取pairarr中的数据;
            //  2.
            //  prevpos在交叉部分之前，nextpos在交叉部分，需要获取lp中的相等read pair进行重新计算
            //  复杂情况1. g中包含prevPosKey对应的unpair，p中有对应的pair，此时应该把这些pair考虑进去
            //  3.
            //  prevpos在交叉部分，nextpos在交叉部分之后，需要获取p中的相等read pair进行重新计算
            //  复杂情况2. p中是否包含prevPosKey对应的unpair
            //  4.
            //  prevpos在交叉部分，nextpos在交叉部分，需要获取lp和p中的相等read pair进行重新计算

            bool addDataToPos = true;
            if (alreadyAdd.find(ck) != alreadyAdd.end()) {
                // 之前已经添加过了，后面就不用再添加数据了，因为同一个位置可能找到两个及以上的unpair数据，处理之前的数据时候可能已经添加了这些数据
                addDataToPos = false;
            } else
                alreadyAdd.insert(ck);

            if (prevPosKey < nextFirstPos) {                   // prevpos在交叉部分之前
                auto &prevPairArr = prevUnpairInfoP->pairArr;  // prevUnpairInfoP肯定不是nullptr
                prevPairArr.push_back(prevFragEnd);
                if (nextPosKey <= prevLastPos && addDataToPos) {  // 第二种情况
                    getEqualRE(prevFragEnd, lpSM.pairs, &prevPairArr);
                }
                // 第一种情况，第二种情况下都会出现，复杂情况一
                auto gPosInfo = g.unpairedPosArr.find(prevPosKey);
                if (gPosInfo != g.unpairedPosArr.end()) {  // 可能g和p有匹配的，刚好和该位点一致
                    auto &gUnpairInfo = gPosInfo->second;
                    auto pPosInfo = pSM.unpairedPosArr.find(nextPosKey);
                    if (pPosInfo != pSM.unpairedPosArr.end()) {
                        auto &pUnpairInfo = pPosInfo->second;
                        for (auto &rn : gUnpairInfo.readNameSet) {  // 遍历每一个readname，看是否有匹配的
                            if (pUnpairInfo.readNameSet.find(rn) != pUnpairInfo.readNameSet.end()) {
                                auto pe = g.unpairedDic[rn].unpairedRE;
                                auto fe = pSM.unpairedDic[rn].unpairedRE;
                                modifyPairedEnds(fe, &pe);
                                prevPairArr.push_back(pe);
                                g.unpairedDic.erase(rn);
                                pSM.unpairedDic.erase(rn);
                            }
                        }
                    }
                }
                recalcPos[ck] = prevPosInfo.taskSeq;
                std::sort(prevPairArr.begin(), prevPairArr.end());
            } else {                                   // prevpos在交叉部分
                if (nextPosKey > prevLastPos) {        // nextpos在交叉部分之后 第三种情况
                    if (nextUnpairInfoP != nullptr) {  // 且在pos点，next task有unpair，这样才把这些数据放到next task里
                        auto &nextPairArr = nextUnpairInfoP->pairArr;
                        nextPairArr.push_back(prevFragEnd);
                        auto &prevPairArr = prevUnpairInfoP->pairArr;
                        prevPairArr.push_back(prevFragEnd);
                        if (addDataToPos) {
                            getEqualRE(prevFragEnd, pSM.pairs, &prevPairArr);
                            getEqualRE(prevFragEnd, pSM.pairs, &nextPairArr);
                        }
                        // 将数据放到next task里,（这个位点以后会可能还会计算到，目前方案是都计算，只是把冗余剔除）
                        recalcPos[ck] = nextPosInfo.taskSeq;

                        std::sort(prevPairArr.begin(), prevPairArr.end());
                        std::sort(nextPairArr.begin(), nextPairArr.end());
                    } else {  // next task在该位点没有unpair，那就把数据放到prev task里
                        auto &prevPairArr = prevUnpairInfoP->pairArr;  // prevUnpairInfoP肯定不是nullptr
                        prevPairArr.push_back(prevFragEnd);
                        if (addDataToPos)  // 第二种情况
                            getEqualRE(prevFragEnd, pSM.pairs, &prevPairArr);
                        recalcPos[ck] = prevPosInfo.taskSeq;
                        std::sort(prevPairArr.begin(), prevPairArr.end());
                    }
                } else {  // 第四种情况
                    if (prevUnpairInfoP == nullptr) {
                        prevUnpairInfoP = &lpSM.unpairedPosArr[prevPosKey];
                        prevUnpairInfoP->taskSeq = lp.taskSeq;
                    }
                    auto &prevPairArr = prevUnpairInfoP->pairArr;
                    prevPairArr.push_back(prevFragEnd);
                    if (addDataToPos) {
                        getEqualRE(prevFragEnd, lpSM.pairs, &prevPairArr);
                        getEqualRE(prevFragEnd, pSM.pairs, &prevPairArr);
                    }
                    recalcPos[ck] = prevPosInfo.taskSeq;
                    std::sort(prevPairArr.begin(), prevPairArr.end());
                }
            }
            pSM.unpairedDic.erase(readName);                               // 在next task里删除该read
        } else if (g.unpairedDic.find(readName) != g.unpairedDic.end()) {  // 在遗留数据中找到了匹配的read
            auto &remainPosInfo = g.unpairedDic[readName];
            auto remainFragEnd = remainPosInfo.unpairedRE;
            int64_t remainPosKey = remainFragEnd.posKey;
            modifyPairedEnds(prevFragEnd, &remainFragEnd);  // 在某些clip情况下，poskey可能是后面的read
            auto &remainUnpairInfo = g.unpairedPosArr[remainPosKey];
            auto &remainPairArr = remainUnpairInfo.pairArr;
            remainPairArr.push_back(remainFragEnd);
            CalcKey ck = {remainPosKey, prevFragEnd.posKey};
            recalcPos[ck] = remainPosInfo.taskSeq;
            std::sort(remainPairArr.begin(), remainPairArr.end());

            g.unpairedDic.erase(readName);
        } else {  // 都没找到，那就保存到遗留数据里
            int64_t prevPosKey = prevFragEnd.posKey;
            g.unpairedDic.insert(prevUnpair);
            addToGlobal.insert(prevPosKey);
        }
    }
    map<int64_t, TaskSeqDupInfo> taskChanged;
    MDSet<int64_t> posProcessed;
    for (auto &e : recalcPos) {
        auto posKey = e.first.read1Pos;
        if (posProcessed.find(posKey) != posProcessed.end())
            continue;
        posProcessed.insert(posKey);
        auto taskSeq = e.second;
        auto &t = taskChanged[taskSeq];
        // 在对应的任务包含的dup idx里修改结果数据
        vector<ReadEnds> *pairArrP = nullptr;
        if (taskSeq < lp.taskSeq)
            pairArrP = &g.unpairedPosArr[posKey].pairArr;
        else
            pairArrP = &lpSM.unpairedPosArr[posKey].pairArr;
        processPairs(*pairArrP, &t.dupIdx, &t.opticalDupIdx, &t.repIdx, &t.notDupIdx, &t.notOpticalDupIdx,
                     &t.notRepIdx);
        if (taskSeq < lp.taskSeq)
            g.unpairedPosArr.erase(posKey);
    }
    // 最后再添加，以防开始赋值，后来这个位置要是又添加了新的数据
    // 放在这里，因为lp中的unpairedPosArr中的readends可能会被修改（比如optical duplicate）
#if 0
    for (auto posKey : addToGlobal) {
        g.unpairedPosArr[posKey] = lpSM.unpairedPosArr[posKey];
    }
#endif
    /* 不需要把p中找不到lp的unpair，放到global中，否则最后找到pair后，还要再执行一次冗余检测，造成重复的冗余索引 */

    // 更新结果
    for (auto &e : taskChanged) {
        auto taskSeq = e.first;
        auto &t = e.second;
        if (taskSeq < lp.taskSeq) {
            refeshTaskDupInfo(t.dupIdx, t.opticalDupIdx, t.repIdx, t.notDupIdx, t.notOpticalDupIdx, t.notRepIdx,
                              g.latterDupIdxArr[taskSeq], g.latterOpticalDupIdxArr[taskSeq], g.latterRepIdxArr[taskSeq],
                              g.latterNotDupIdxArr[taskSeq], g.latterNotOpticalDupIdxArr[taskSeq],
                              g.latterNotRepIdxArr[taskSeq]);
        } else if (taskSeq == lp.taskSeq) {
            refreshMarkDupData(t.dupIdx, t.opticalDupIdx, t.repIdx, t.notDupIdx, t.notOpticalDupIdx, t.notRepIdx, lp,
                               p);
        } else {
            refreshMarkDupData(t.dupIdx, t.opticalDupIdx, t.repIdx, t.notDupIdx, t.notOpticalDupIdx, t.notRepIdx, p,
                               lp);  // 把结果放到p中
        }
    }

    // 将dupidx放进全局数据
    g.latterDupIdxArr.push_back(DPSet<DupInfo>());
    g.latterOpticalDupIdxArr.push_back(MDSet<int64_t>());
    g.latterRepIdxArr.push_back(DPSet<DupInfo>());
    g.latterNotDupIdxArr.push_back(MDSet<int64_t>());
    g.latterNotOpticalDupIdxArr.push_back(MDSet<int64_t>());
    g.latterNotRepIdxArr.push_back(MDSet<int64_t>());

    g.dupIdxArr.push_back(vector<DupInfo>());
    auto &vIdx = g.dupIdxArr.back();
    lp.pairDupIdx.insert(lp.fragDupIdx.begin(), lp.fragDupIdx.end());
    vIdx.insert(vIdx.end(), lp.pairDupIdx.begin(), lp.pairDupIdx.end());
    std::sort(vIdx.begin(), vIdx.end());

    g.opticalDupIdxArr.push_back(vector<int64_t>());
    auto &vOpticalIdx = g.opticalDupIdxArr.back();
    vOpticalIdx.insert(vOpticalIdx.end(), lp.pairOpticalDupIdx.begin(), lp.pairOpticalDupIdx.end());
    std::sort(vOpticalIdx.begin(), vOpticalIdx.end());

    g.repIdxArr.push_back(vector<DupInfo>());
    auto &vRepIdx = g.repIdxArr.back();
    vRepIdx.insert(vRepIdx.end(), lp.pairRepIdx.begin(), lp.pairRepIdx.end());
    std::sort(vRepIdx.begin(), vRepIdx.end());
}

static void *pipeRead(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;
    BamBufType inBamBuf(nsgv::gMdArg.DUPLEX_IO);
    inBamBuf.Init(nsgv::gInBamFp, nsgv::gInBamHeader, nsgv::gMdArg.MAX_MEM);
    int64_t readNumSum = 0;
    while (1) {
        PROF_START(read_wait);
        yarn::POSSESS(pipeArg.readSig);
        yarn::WAIT_FOR(pipeArg.readSig, yarn::NOT_TO_BE, 1);  // 只有一个坑，因为在bambuf内部支持异步读取
        PROF_END(gprof[GP_read_wait], read_wait);
        size_t readNum = 0;
        PROF_START(read);
        if (inBamBuf.ReadStat() >= 0)
            readNum = inBamBuf.ReadBam();  // 读入新一轮的数据
        PROF_END(gprof[GP_read], read);
        if (readNum < 1) {
            pipeArg.readFinish = 1;
            yarn::TWIST(pipeArg.readSig, yarn::BY, 1);  // 读入了一轮数据
            break;
        }
        cout << "read num: " << readNum << "\t" << readNumSum << '\t' << pipeArg.readOrder << endl;

        pipeArg.readData.bams = inBamBuf.GetBamArr();

        pipeArg.readData.bamStartIdx = readNumSum;
        pipeArg.readData.taskSeq = pipeArg.readOrder;

        readNumSum += readNum;
        inBamBuf.ClearAll();            // 清理上一轮读入的数据
        pipeArg.readOrder += 1;         // for next
        yarn::TWIST(pipeArg.readSig, yarn::BY, 1);  // 读入了一轮数据
    }
    cout << "total reads: " << readNumSum << endl;
    return 0;
}
static void *pipeGenRE(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;
    auto &genREData = pipeArg.genREData;
    // init generate read ends data by num_thread
    int genREThread = pipeArg.numThread;
    // / 4;
    for (int i = 0; i < pipeArg.GENBUFNUM; ++i) {
        genREData[i].Init(genREThread);
    }

    while (1) {
        PROF_START(gen_wait);
        yarn::POSSESS(pipeArg.readSig);
        yarn::WAIT_FOR(pipeArg.readSig, yarn::NOT_TO_BE, 0);  // 等待有数据
        yarn::POSSESS(pipeArg.genRESig);
        PROF_END(gprof[GP_gen_wait], gen_wait);

        yarn::WAIT_FOR(pipeArg.genRESig, yarn::NOT_TO_BE, pipeArg.GENBUFNUM);  // 有BUFNUM个坑
        yarn::RELEASE(pipeArg.genRESig);                                       // 因为有不止一个位置，所以要释放

        if (pipeArg.readFinish) {
            yarn::POSSESS(pipeArg.genRESig);
            pipeArg.genREFinish = 1;
            yarn::TWIST(pipeArg.genRESig, yarn::BY, 1);
            yarn::TWIST(pipeArg.readSig, yarn::BY, -1);
            break;
        }
        /* 处理bam，生成readends */
        //        cout << "genRE order: " << pipeArg.genREOrder << "\t" << pipeArg.readData.bamStartIdx << endl;
        PROF_START(gen);
        doGenRE(pipeArg);
        // usleep(200000);
        PROF_END(gprof[GP_gen], gen);

        yarn::POSSESS(pipeArg.genRESig);
        pipeArg.genREOrder += 1;
        yarn::TWIST(pipeArg.genRESig, yarn::BY, 1);
        yarn::TWIST(pipeArg.readSig, yarn::BY, -1);  // 使用了一次读入的数据
    }
    cout << "pipe gen reads" << endl;
    return 0;
}
static void *pipeSort(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;

    while (1) {
        PROF_START(sort_wait);
        yarn::POSSESS(pipeArg.genRESig);
        yarn::WAIT_FOR(pipeArg.genRESig, yarn::NOT_TO_BE, 0);  // 等待有数据
        yarn::RELEASE(pipeArg.genRESig);
        PROF_END(gprof[GP_sort_wait], sort_wait);

        yarn::POSSESS(pipeArg.sortSig);
        yarn::WAIT_FOR(pipeArg.sortSig, yarn::NOT_TO_BE, pipeArg.SORTBUFNUM);  // 有BUFNUM个位置
        yarn::RELEASE(pipeArg.sortSig);

        if (pipeArg.genREFinish) {
            // 处理完剩余数据
            while (pipeArg.sortOrder < pipeArg.genREOrder) {
                PROF_START(sort);
                doSort(pipeArg);
                PROF_END(gprof[GP_sort], sort);
                pipeArg.sortOrder += 1;
            }
            yarn::POSSESS(pipeArg.sortSig);
            pipeArg.sortFinish = 1;
            yarn::TWIST(pipeArg.sortSig, yarn::BY, 1);
            break;
        }
        /* 排序 readends */
        //        cout << "sort order: " << pipeArg.sortOrder << endl;
        PROF_START(sort);
        doSort(pipeArg);
        PROF_END(gprof[GP_sort], sort);

        yarn::POSSESS(pipeArg.genRESig);
        yarn::TWIST(pipeArg.genRESig, yarn::BY, -1);

        yarn::POSSESS(pipeArg.sortSig);
        pipeArg.sortOrder += 1;
        yarn::TWIST(pipeArg.sortSig, yarn::BY, 1);
    }
    cout << "end pipe sort" << endl;
    return 0;
}
static void *pipeMarkDup(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;

    while (1) {
        PROF_START(markdup_wait);
        yarn::POSSESS(pipeArg.sortSig);
        yarn::WAIT_FOR(pipeArg.sortSig, yarn::NOT_TO_BE, 0);  // 等待有数据
        yarn::RELEASE(pipeArg.sortSig);
        PROF_END(gprof[GP_markdup_wait], markdup_wait);

        yarn::POSSESS(pipeArg.markDupSig);
        yarn::WAIT_FOR(pipeArg.markDupSig, yarn::NOT_TO_BE, pipeArg.MARKBUFNUM);
        yarn::RELEASE(pipeArg.markDupSig);

        if (pipeArg.sortFinish) {
            // 应该还得处理剩余的数据
            while (pipeArg.markDupOrder < pipeArg.sortOrder) {
                PROF_START(markdup);
                doMarkDup(pipeArg);
                PROF_END(gprof[GP_markdup], markdup);
                pipeArg.markDupOrder += 1;
            }
            yarn::POSSESS(pipeArg.markDupSig);
            pipeArg.markDupFinish = 1;
            yarn::TWIST(pipeArg.markDupSig, yarn::BY, 2);
            break;
        }
        /* 冗余检测 readends */
        // cout << "markdup order: " << pipeArg.markDupOrder << endl;
        PROF_START(markdup);
        doMarkDup(pipeArg);
        PROF_END(gprof[GP_markdup], markdup);
        yarn::POSSESS(pipeArg.sortSig);
        yarn::TWIST(pipeArg.sortSig, yarn::BY, -1);

        yarn::POSSESS(pipeArg.markDupSig);
        pipeArg.markDupOrder += 1;
        yarn::TWIST(pipeArg.markDupSig, yarn::BY, 1);
    }
    cout << "end pipe markdup" << endl;
    return 0;
}
static void *pipeIntersect(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;
    pipeArg.intersectOrder = 1;
    while (1) {
        PROF_START(intersect_wait);
        yarn::POSSESS(pipeArg.markDupSig);
        yarn::WAIT_FOR(pipeArg.markDupSig, yarn::TO_BE_MORE_THAN, 1);  // 等待有数据
        yarn::RELEASE(pipeArg.markDupSig);
        PROF_END(gprof[GP_intersect_wait], intersect_wait);

        if (pipeArg.markDupFinish) {
            while (pipeArg.intersectOrder < pipeArg.markDupOrder) {
                PROF_START(intersect);
                doIntersect(pipeArg);
                PROF_END(gprof[GP_intersect], intersect);
                pipeArg.intersectOrder += 1;
            }
            break;
        }
        /* 交叉数据处理 readends */
        PROF_START(intersect);
        doIntersect(pipeArg);
        PROF_END(gprof[GP_intersect], intersect);

        yarn::POSSESS(pipeArg.markDupSig);
        yarn::TWIST(pipeArg.markDupSig, yarn::BY, -1);

        pipeArg.intersectOrder += 1;
    }
    cout << "end pipe intersect" << endl;
    return 0;
}

/* 当所有任务结束后，global data里还有未处理的数据 */
static void mergeAllTask(PipelineArg &pipeArg) {
    MarkDupData &lp = pipeArg.markDupData[(pipeArg.markDupOrder - 1) % pipeArg.MARKBUFNUM];
    IntersectData &g = pipeArg.intersectData;
    SortMarkData &lpSM = *(SortMarkData *)lp.dataPtr;
    // 遗留的未匹配的pair
    for (auto &prevUnpair : lpSM.unpairedDic) {  // 遍历上一个任务中的每个未匹配的read
        auto &readName = prevUnpair.first;
        auto &prevPosInfo = prevUnpair.second;
        auto prevFragEnd = prevPosInfo.unpairedRE;  // 未匹配的read end

        if (g.unpairedDic.find(readName) != g.unpairedDic.end()) {  // 在遗留数据中找到了匹配的read
            auto &remainPosInfo = g.unpairedDic[readName];
            auto remainFragEnd = remainPosInfo.unpairedRE;
            int64_t remainPosKey = remainFragEnd.posKey;
            modifyPairedEnds(prevFragEnd, &remainFragEnd);  // 在某些clip情况下，poskey可能是后面的read
            auto &remainUnpairInfo = g.unpairedPosArr[remainPosKey];

            remainUnpairInfo.pairArr.push_back(remainFragEnd);
            g.unpairedDic.erase(readName);
        } else {
            g.unpairedDic.insert(prevUnpair);  // 用来记录没有匹配的read个数
        }
    }

    map<int64_t, TaskSeqDupInfo> taskChanged;
    for (auto &e : g.unpairedPosArr) {
        auto posKey = e.first;
        auto taskSeq = e.second.taskSeq;
        auto &t = taskChanged[taskSeq];
        auto &arr = g.unpairedPosArr[posKey].pairArr;

        if (arr.size() > 1) {
            std::sort(arr.begin(), arr.end());
            processPairs(arr, &t.dupIdx, &t.opticalDupIdx, &t.repIdx, &t.notDupIdx, &t.notOpticalDupIdx, &t.notRepIdx);
        }
    }
    // 更新结果
    vector<int64_t> addDup;
    map<int64_t, int64_t> ndPosVal;
    for (auto &e : taskChanged) {
        auto taskSeq = e.first;
        auto &t = e.second;
        refeshTaskDupInfo(t.dupIdx, t.opticalDupIdx, t.repIdx, t.notDupIdx, t.notOpticalDupIdx, t.notRepIdx,
                          g.latterDupIdxArr[taskSeq], g.latterOpticalDupIdxArr[taskSeq], g.latterRepIdxArr[taskSeq],
                          g.latterNotDupIdxArr[taskSeq], g.latterNotOpticalDupIdxArr[taskSeq],
                          g.latterNotRepIdxArr[taskSeq]);
    }
    g.unpairedPosArr.clear();

    // 将dupidx放进全局数据
    vector<DupInfo> cacheDupIdx;
    vector<DupInfo> midArr;
    vector<int64_t> intCacheDupIdx;
    vector<int64_t> intMidArr;
    for (int i = 0; i < (int)g.dupIdxArr.size() - 1; ++i) {
        refeshFinalTaskDupInfo(g.latterDupIdxArr[i], g.latterNotDupIdxArr[i], g.dupIdxArr[i], cacheDupIdx, midArr);
    }
    for (int i = 0; i < (int)g.opticalDupIdxArr.size() - 1; ++i)
        refeshFinalTaskDupInfo(g.latterOpticalDupIdxArr[i], g.latterNotOpticalDupIdxArr[i], g.opticalDupIdxArr[i],
                               intCacheDupIdx, intMidArr);
    for (int i = 0; i < (int)g.repIdxArr.size() - 1; ++i)
        refeshFinalTaskDupInfo(g.latterRepIdxArr[i], g.latterNotRepIdxArr[i], g.repIdxArr[i], cacheDupIdx, midArr);

    g.dupIdxArr.push_back(vector<DupInfo>());
    auto &vIdx = g.dupIdxArr.back();
    lp.pairDupIdx.insert(lp.fragDupIdx.begin(), lp.fragDupIdx.end());
    vIdx.insert(vIdx.end(), lp.pairDupIdx.begin(), lp.pairDupIdx.end());
    std::sort(vIdx.begin(), vIdx.end());

    g.opticalDupIdxArr.push_back(vector<int64_t>());
    auto &vOpticalIdx = g.opticalDupIdxArr.back();
    vOpticalIdx.insert(vOpticalIdx.end(), lp.pairOpticalDupIdx.begin(), lp.pairOpticalDupIdx.end());
    std::sort(vOpticalIdx.begin(), vOpticalIdx.end());

    g.repIdxArr.push_back(vector<DupInfo>());
    auto &vRepIdx = g.repIdxArr.back();
    vRepIdx.insert(vRepIdx.end(), lp.pairRepIdx.begin(), lp.pairRepIdx.end());
    std::sort(vRepIdx.begin(), vRepIdx.end());
}

static void parallelPipeline() {
    PipelineArg &pipeArg = nsgv::gPipe;
    pipeArg.numThread = nsgv::gMdArg.NUM_THREADS;

    pthread_t tidArr[5];
    pthread_create(&tidArr[0], 0, pipeRead, &pipeArg);
    pthread_create(&tidArr[1], 0, pipeGenRE, &pipeArg);
    pthread_create(&tidArr[2], 0, pipeSort, &pipeArg);
    pthread_create(&tidArr[3], 0, pipeMarkDup, &pipeArg);
    pthread_create(&tidArr[4], 0, pipeIntersect, &pipeArg);
    for (int i = 0; i < 5; ++i) pthread_join(tidArr[i], 0);
    PROF_START(merge_result);
    mergeAllTask(pipeArg);
    PROF_END(gprof[GP_merge_result], merge_result);

    int64_t dupNum = 0;
    int64_t opticalDupNum = 0;
    for (auto &arr : pipeArg.intersectData.dupIdxArr) dupNum += arr.size();
    for (auto &arr : pipeArg.intersectData.opticalDupIdxArr) opticalDupNum += arr.size();

    map<int64_t, int> dup;
#if 0
    int taskSeq = 0;
    for (auto &arr : pipeArg.intersectData.dupIdxArr) {
        for (auto idx : arr) {
            if (dup.find(idx.idx) != dup.end()) {
                //if (taskSeq - 1 > dup[idx])
                cout << "dup index: " << dup[idx] << '\t' << taskSeq << '\t' << idx << endl;
            }
            dup[idx.idx] = taskSeq;
        }
        // cout << taskSeq << "\t" << arr.size() << endl;
        taskSeq++;
    }
#endif

    cout << "Final read  order: " << pipeArg.readOrder << endl;
    cout << "Final gen   order: " << pipeArg.genREOrder << endl;
    cout << "Final sort  order: " << pipeArg.sortOrder << endl;
    cout << "Final mark  order: " << pipeArg.markDupOrder << endl;
    cout << "Final inter order: " << pipeArg.intersectOrder << endl;

//    cout << "w read time     : " << tm_arr[10].acc_seconds_elapsed() << endl;
//    cout << "w gen time      : " << tm_arr[11].acc_seconds_elapsed() << endl;
//    cout << "w sort time     : " << tm_arr[12].acc_seconds_elapsed() << endl;
//    cout << "w markdup time  : " << tm_arr[13].acc_seconds_elapsed() << endl;
//    cout << "w intersect time: " << tm_arr[14].acc_seconds_elapsed() << endl;
//
//    cout << "w1 gen time      : " << tm_arr[21].acc_seconds_elapsed() << endl;
//    cout << "w1 sort time     : " << tm_arr[22].acc_seconds_elapsed() << endl;
//    cout << "w1 markdup time  : " << tm_arr[23].acc_seconds_elapsed() << endl;
//
//    cout << "read time     : " << tm_arr[0].acc_seconds_elapsed() << endl;
//    cout << "gen time      : " << tm_arr[1].acc_seconds_elapsed() << endl;
//    cout << "sort time     : " << tm_arr[2].acc_seconds_elapsed() << endl;
//    cout << "markdup time  : " << tm_arr[3].acc_seconds_elapsed() << endl;
//    cout << "intersect time: " << tm_arr[4].acc_seconds_elapsed() << endl;
//
//    cout << "copy time: " << tm_arr[5].acc_seconds_elapsed() << endl;
//    cout << "merge al6 time: " << tm_arr[6].acc_seconds_elapsed() << endl;
//
//    cout << "dup num         : " << dupNum << "\t" << dup.size() << endl;
//    cout << "optical dup num : " << opticalDupNum / 2 << "\t" << opticalDupNum << endl;

}

/* 并行流水线方式处理数据，标记冗余 */
void pipelineMarkDups() {
    if (nsgv::gMdArg.NUM_THREADS > 1)
        return parallelPipeline();

    PipelineArg &pipeArg = nsgv::gPipe;
    pipeArg.numThread = nsgv::gMdArg.NUM_THREADS;
    BamBufType inBamBuf(nsgv::gMdArg.DUPLEX_IO);
    inBamBuf.Init(nsgv::gInBamFp, nsgv::gInBamHeader, nsgv::gMdArg.MAX_MEM);
    int64_t readNumSum = 0;
    for (int i = 0; i < pipeArg.GENBUFNUM; ++i) {
        pipeArg.genREData[i].Init(pipeArg.numThread);
    }
    pipeArg.intersectOrder = 1;  // do intersect 从1开始
    while (1) {
        size_t readNum = 0;
        if (inBamBuf.ReadStat() >= 0)
            readNum = inBamBuf.ReadBam();  // 读入新一轮的数据
        if (readNum < 1) {
            break;
        }
        cout << "read num: " << readNum << "\t" << readNumSum << '\t' << pipeArg.readOrder << endl;

        pipeArg.readData.bams = inBamBuf.GetBamArr();
        pipeArg.readData.bamStartIdx = readNumSum;
        pipeArg.readData.taskSeq = pipeArg.readOrder;
        // 1. do generate read ends
        doGenRE(pipeArg);
        pipeArg.genREOrder += 1;
        // 2. do sort
        doSort(pipeArg);
        pipeArg.sortOrder += 1;
        // 3. do markduplicate
        doMarkDup(pipeArg);
        pipeArg.markDupOrder += 1;
        // 4. do intersect data
        if (pipeArg.markDupOrder > 1) {
            doIntersect(pipeArg);
            pipeArg.intersectOrder += 1;
        }

        readNumSum += readNum;
        inBamBuf.ClearAll();     // 清理上一轮读入的数据
        pipeArg.readOrder += 1;  // for next
    }
    mergeAllTask(pipeArg);

    int64_t dupNum = 0;
    int64_t opticalDupNum = 0;
    for (auto &arr : pipeArg.intersectData.dupIdxArr) dupNum += arr.size();
    for (auto &arr : pipeArg.intersectData.opticalDupIdxArr) opticalDupNum += arr.size();

    map<int64_t, int> dup;
#if 0
    int taskSeq = 0;
    for (auto &arr : pipeArg.intersectData.dupIdxArr) {
        for (auto idx : arr) {
            if (dup.find(idx.idx) != dup.end()) {
                cout << "dup index: " << dup[idx] << '\t' << taskSeq << '\t' << idx << endl;
            }
            dup[idx.idx] = taskSeq;
        }
        taskSeq++;
    }
#endif

    cout << "total reads: " << readNumSum << endl;
    cout << "Final read  order: " << pipeArg.readOrder << endl;
    cout << "Final gen   order: " << pipeArg.genREOrder << endl;
    cout << "Final sort  order: " << pipeArg.sortOrder << endl;
    cout << "Final mark  order: " << pipeArg.markDupOrder << endl;
    cout << "Final inter order: " << pipeArg.intersectOrder << endl;
    cout << "dup num         : " << dupNum << "\t" << dup.size() << endl;
    cout << "optical dup num : " << opticalDupNum / 2 << "\t" << opticalDupNum << endl;
}