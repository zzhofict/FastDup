#include "md_pipeline.h"

#include <klib/kthread.h>
#include <pthread.h>
#include <spdlog/spdlog.h>

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

using std::vector;
using namespace std;

namespace nsgv {

extern MarkDupsArg gMdArg;           // 用来测试性能
extern samFile *gInBamFp;            // 输入的bam文件
extern sam_hdr_t *gInBamHeader;      // 输入的bam文件头信息
extern DuplicationMetrics gMetrics;  // 统计信息
extern vector<ReadNameParser> gNameParsers;
extern DupResult gDupRes;
extern PipelineArg gPipe;
};  // namespace nsgv

/* 处理一组pairend的readends，标记冗余, 这个函数需要串行运行，因为需要做一些统计*/
static void markDupsForPairs(vector<const ReadEnds *> &vpRe, DPSet<DupInfo> *dupIdx, MDSet<int64_t> *opticalDupIdx,
                             DPSet<DupInfo> *repIdx, MDSet<int64_t> *notDupIdx = nullptr,
                             MDSet<int64_t> *notOpticalDupIdx = nullptr, MDSet<int64_t> *notRepIdx = nullptr) {
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
    auto range = std::equal_range(src.begin(), src.end(), re, ReadEnds::CorLittleThan);  // 只比对位点
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


/* 获取交叉部分的数据 */
static inline void getIntersectData(vector<ReadEnds> &leftArr, vector<ReadEnds> &rightArr, vector<ReadEnds> *dst,
                                    bool isPairCmp = false) {
    if (leftArr.empty() || rightArr.empty()) {
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
    if (isPairCmp)
        std::sort(dst->begin(), dst->end(), ReadEnds::PairLittleThan);
    else
        std::sort(dst->begin(), dst->end(), ReadEnds::FragLittleThan);
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

    PROF_START(gen);
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
    PROF_END(tprof[TP_gen][tid], gen);

    PROF_START(sort_frag);
    sort(frags.begin(), frags.end(), ReadEnds::FragLittleThan);
    PROF_END(tprof[TP_sort_frag][tid], sort_frag);

    PROF_START(sort_pair);
    sort(pairs.begin(), pairs.end(), ReadEnds::PairLittleThan);
    PROF_END(tprof[TP_sort_pair][tid], sort_pair);
}

static void doGenRE(PipelineArg &pipeArg) {
    GenREData &genREData = pipeArg.genREData[pipeArg.genREOrder % pipeArg.GENBUFNUM];
    ReadData &readData = pipeArg.readData;

    // generate read ends
    const int numThread = pipeArg.numThread;

    kt_for(numThread, mtGenerateReadEnds, &pipeArg, numThread);
    // 用来在genRE阶段计算一些Sort阶段应该计算的数据，保持每个step用时更平衡
    // 轮询每个线程中未找到匹配的read，找到匹配的
    genREData.unpairedDic.clear();
    vector<ReadEnds> &pairs = genREData.pairsArr[numThread];
    pairs.clear();

    int testNum = 0;
    for (int i = 0; i < numThread; ++i) {
        testNum += genREData.unpairedDicArr[i].size();
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
    sort(pairs.begin(), pairs.end(), ReadEnds::PairLittleThan);
}
// end for step 2 generate read ends

// for step-3 sort
static void doSort(PipelineArg &pipeArg) {
    // return;
    GenREData &genREData = pipeArg.genREData[pipeArg.sortOrder % pipeArg.GENBUFNUM];
    SortData &sortData = pipeArg.sortData[pipeArg.sortOrder % pipeArg.SORTBUFNUM];
    SortMarkData &smd = *(SortMarkData *)sortData.dataPtr;

    smd.unpairedDic = genREData.unpairedDic;

    smd.pairs.clear();
    smd.frags.clear();
    const ReadEnds *pRE;
    ReadEndsHeap<REFragGreaterThan> fragsHeap;
    ReadEndsHeap<REPairGreaterThan> pairsHeap;
    PROF_START(sort_pair);
    pairsHeap.Init(&genREData.pairsArr);
    while ((pRE = pairsHeap.Pop()) != nullptr) {
        smd.pairs.push_back(*pRE);
    }
    PROF_END(gprof[GP_sort_pair], sort_pair);
    PROF_START(sort_frag);
    fragsHeap.Init(&genREData.fragsArr);
    while ((pRE = fragsHeap.Pop()) != nullptr) {
        smd.frags.push_back(*pRE);
    }
    PROF_END(gprof[GP_sort_frag], sort_frag);
}
// for step-4 sort
static void doMarkDup(PipelineArg &pipeArg) {
    MarkDupData &mdData = pipeArg.markDupData[pipeArg.markDupOrder % pipeArg.MARKBUFNUM];
    SortData &sortData = pipeArg.sortData[pipeArg.markDupOrder % pipeArg.SORTBUFNUM];
    mdData.taskSeq = pipeArg.markDupOrder;
    mdData.clear();

    auto tmpPtr = mdData.dataPtr;
    mdData.dataPtr = sortData.dataPtr;
    sortData.dataPtr = tmpPtr;

    SortMarkData &smd = *(SortMarkData *)mdData.dataPtr;
    //  先处理 pair
    PROF_START(markdup_pair);
    processPairs(smd.pairs, &mdData.pairDupIdx, &mdData.pairOpticalDupIdx, &mdData.pairRepIdx);
    PROF_END(gprof[GP_markdup_pair], markdup_pair);
    // 再处理frag
    PROF_START(markdup_frag);
    processFrags(smd.frags, &mdData.fragDupIdx);
    PROF_END(gprof[GP_markdup_frag], markdup_frag);
}

template <typename T>
static void refreshDupIdx(T &srcArr, T &insArr) {
    for (auto dup : srcArr) {
        insArr.insert(dup);
    }
}
template <typename T1, typename T2>
static void refreshNotDupIdx(T1 &srcArr, T2 &delArr) {
    for (auto dup : srcArr) {
        delArr.erase(dup);
    }
}

static void refreshMarkDupData(DPSet<DupInfo> &dupIdx, MDSet<int64_t> &opticalDupIdx, DPSet<DupInfo> &repIdx,
                               MDSet<int64_t> &notDupIdx, MDSet<int64_t> &notOpticalDupIdx, MDSet<int64_t> &notRepIdx,
                               MarkDupData &lp) {
    refreshDupIdx(dupIdx, lp.pairDupIdx);
    refreshDupIdx(opticalDupIdx, lp.pairOpticalDupIdx);
    refreshDupIdx(repIdx, lp.pairRepIdx);
    refreshNotDupIdx(notDupIdx, lp.pairDupIdx);
    refreshNotDupIdx(notOpticalDupIdx, lp.pairOpticalDupIdx);
    refreshNotDupIdx(notRepIdx, lp.pairRepIdx);
}

static void refreshMarkDupData(DPSet<DupInfo> &dupIdx, MDSet<int64_t> &opticalDupIdx, DPSet<DupInfo> &repIdx,
                               MDSet<int64_t> &notDupIdx, MDSet<int64_t> &notOpticalDupIdx, MDSet<int64_t> &notRepIdx,
                               MarkDupData &lp, MarkDupData &p) {
    refreshDupIdx(dupIdx, lp.pairDupIdx);
    refreshDupIdx(opticalDupIdx, lp.pairOpticalDupIdx);
    refreshDupIdx(repIdx, lp.pairRepIdx);
    refreshNotDupIdx(notDupIdx, lp.pairDupIdx);
    refreshNotDupIdx(notOpticalDupIdx, lp.pairOpticalDupIdx);
    refreshNotDupIdx(notRepIdx, lp.pairRepIdx);
    refreshNotDupIdx(notDupIdx, p.pairDupIdx);
    refreshNotDupIdx(notOpticalDupIdx, p.pairOpticalDupIdx);
    refreshNotDupIdx(notRepIdx, p.pairRepIdx);
    refreshNotDupIdx(dupIdx, p.pairDupIdx);
    refreshNotDupIdx(opticalDupIdx, p.pairOpticalDupIdx);
    refreshNotDupIdx(repIdx, p.pairRepIdx);
}

// 处理相邻数据块之间重叠的部分
static void processIntersectFragPairs(MarkDupData &lp, MarkDupData &cp) {
    SortMarkData &lsm = *(SortMarkData *)lp.dataPtr;
    SortMarkData &csm = *(SortMarkData *)cp.dataPtr;

    vector<ReadEnds> reArr;
    DPSet<DupInfo> dupIdx;
    MDSet<int64_t> opticalDupIdx;
    DPSet<DupInfo> repIdx;
    MDSet<int64_t> notOpticalDupIdx;
    MDSet<int64_t> notDupIdx;
    MDSet<int64_t> notRepIdx;
    // 先处理重叠的frags
    getIntersectData(lsm.frags, csm.frags, &reArr);
    processFrags(reArr, &dupIdx, &notDupIdx);
    refreshDupIdx(dupIdx, lp.fragDupIdx);
    refreshNotDupIdx(dupIdx, cp.fragDupIdx);
    refreshNotDupIdx(notDupIdx, lp.fragDupIdx);
    refreshNotDupIdx(notDupIdx, cp.fragDupIdx);

    // 再处理重叠的pairs
    reArr.clear();
    dupIdx.clear();
    notDupIdx.clear();

    getIntersectData(lsm.pairs, csm.pairs, &reArr, true);
    processPairs(reArr, &dupIdx, &opticalDupIdx, &repIdx, &notDupIdx, &notOpticalDupIdx, &notRepIdx);
    refreshMarkDupData(dupIdx, opticalDupIdx, repIdx, notDupIdx, notOpticalDupIdx, notRepIdx, cp, lp); // 放在cp里，因为后面global里可能有相同的dup，防止多次出现
}

// 在相邻的数据块之间寻找未匹配的readends, 找到匹配的放到lp里
static void findUnpairedInDatas(MarkDupData &lp, MarkDupData &cp) {
    auto &interPairedData = lp.ckeyReadEndsMap;
    SortMarkData &lsm = *(SortMarkData *)lp.dataPtr;
    SortMarkData &csm = *(SortMarkData *)cp.dataPtr;
    for (auto itr = lsm.unpairedDic.begin(); itr != lsm.unpairedDic.end(); ) {  // 遍历上一个任务中的每个未匹配的read
        auto &lastUnpair = *itr;
        auto &readName = lastUnpair.first;
        auto &lastUnpairInfo = lastUnpair.second;
        auto lastRE = lastUnpairInfo.unpairedRE;                        // 未匹配的read end
        if (csm.unpairedDic.find(readName) != csm.unpairedDic.end()) {  // 在当前这个任务里找到了这个未匹配的read
            auto &curUnpairInfo = csm.unpairedDic[readName];
            auto &curRE = curUnpairInfo.unpairedRE;
            modifyPairedEnds(curRE, &lastRE); // lastRE当做找到匹配后，完整的ReadEnds
            CalcKey ck(lastRE);
            auto &pairArr = interPairedData[ck];
            pairArr.push_back(lastRE);
            // 从dict中清除配对后的readends
            csm.unpairedDic.erase(readName);
            itr = lsm.unpairedDic.erase(itr);
        } else {
            ++itr;
        }
    }
}

// 在global和lp中寻找未匹配的readends, 找到匹配的放到global里
static void findUnpairedInGlobal(IntersectData &g, MarkDupData &lp) {
    auto &interPairedData = g.ckeyReadEndsMap;
    SortMarkData &lsm = *(SortMarkData *)lp.dataPtr;
    for (auto itr = lsm.unpairedDic.begin(); itr != lsm.unpairedDic.end();) {  // 遍历上一个任务中的每个未匹配的read
        auto &lastUnpair = *itr;
        auto &readName = lastUnpair.first;
        auto &lastUnpairInfo = lastUnpair.second;
        auto lastRE = lastUnpairInfo.unpairedRE;                    // 未匹配的read end
        if (g.unpairedDic.find(readName) != g.unpairedDic.end()) {  // 在global里找到了这个未匹配的read
            auto &gUnpairInfo = g.unpairedDic[readName];
            auto &gRE = gUnpairInfo.unpairedRE;
            modifyPairedEnds(lastRE, &gRE);  // gRE当做找到匹配后，完整的ReadEnds
            CalcKey ck(gRE);
            auto &pairArr = interPairedData[ck];
            pairArr.push_back(gRE);
            // 从dict中清除配对后的readends
            g.unpairedDic.erase(readName);
            itr = lsm.unpairedDic.erase(itr);
        } else {
            ++itr;
        }
    }
}

static void putDupinfoToGlobal(IntersectData &g, MarkDupData &lp) {
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

// for step-5 handle intersect data
static void doIntersect(PipelineArg &pipeArg) {
    // spdlog::info("intersect order: {}", pipeArg.intersectOrder);
    const int kInitIntersectOrder = 1;
    IntersectData &g = pipeArg.intersectData;
    MarkDupData &lp = pipeArg.markDupData[(pipeArg.intersectOrder - 1) % pipeArg.MARKBUFNUM];
    MarkDupData &cp = pipeArg.markDupData[(pipeArg.intersectOrder) % pipeArg.MARKBUFNUM];

    SortMarkData &lsm = *(SortMarkData *)lp.dataPtr;
    SortMarkData &csm = *(SortMarkData *)cp.dataPtr;

    // 处理相邻数据块之间重叠的部分
    processIntersectFragPairs(lp, cp);

    // 检查确保lp和np之间没有数据交叉
    int64_t lastLeft = INT64_MIN, lastRight = INT64_MAX, curLeft = INT64_MAX, curRight = INT64_MAX;
    if (lsm.pairs.size() > 0) {
        lastLeft = lsm.frags[0].Left();
        lastRight = lsm.frags.back().Left();
    }
    if (csm.pairs.size() > 0) {
        curLeft = csm.frags[0].Left();
        curRight = csm.frags.back().Left();
    }


    if (g.rightPos >= curLeft) {
        spdlog::error("previous data can not contain readends included by next data block! {} {} {} {} {} {}",
                      lp.taskSeq, cp.taskSeq, g.rightPos, curLeft, lsm.pairs.size(), csm.pairs.size());
    }
    g.rightPos = lastRight;

    findUnpairedInDatas(lp, cp);  // 找到的匹配放到lp里
    findUnpairedInGlobal(g, cp);  // 把cp中未匹配的放到global里保存

    // 处理lp中的新找到的匹配
    TaskSeqDupInfo t;
    for (auto itr = lp.ckeyReadEndsMap.begin(); itr != lp.ckeyReadEndsMap.end();) {
        auto &ckVal = *itr;
        auto &ck = ckVal.first;
        auto &pairArr = ckVal.second;
        getEqualRE(pairArr[0], lsm.pairs, &pairArr); // 如果不在计算范围内，会放在global里
        if (ck.Right() <= lastRight) { // 必须在当前数据块的范围内, 才进行处理
            if (ck.Left() >= curLeft) {  // 在交叉的范围内才去加上这些在cp中的数据
                getEqualRE(pairArr[0], csm.pairs, &pairArr);
            }
            // 在global里找一找ck
            auto gitr = g.ckeyReadEndsMap.find(ck);
            if (gitr != g.ckeyReadEndsMap.end()) {
                auto &gPairArr = gitr->second;
                pairArr.insert(pairArr.end(), gPairArr.begin(), gPairArr.end());
                g.ckeyReadEndsMap.erase(gitr);
            }
            sort(pairArr.begin(), pairArr.end(), ReadEnds::PairLittleThan);
            processPairs(pairArr, &t.dupIdx, &t.opticalDupIdx, &t.repIdx, &t.notDupIdx, &t.notOpticalDupIdx, &t.notRepIdx);
            itr = lp.ckeyReadEndsMap.erase(itr);
        } else {
            ++itr;
        }
    }

    // 处理找到匹配的global数据
    for (auto itr = g.ckeyReadEndsMap.begin(); itr != g.ckeyReadEndsMap.end();) {
        auto &ckVal = *itr;
        auto &ck = ckVal.first;
        auto &pairArr = ckVal.second;
        if (ck.Left() >= lastLeft) {
            getEqualRE(pairArr[0], lsm.pairs, &pairArr);
        }
        if (ck.Right() <= lastRight) {  // 只有在这个范围内，对应位点的所有reads才完全都包含了
            sort(pairArr.begin(), pairArr.end(), ReadEnds::PairLittleThan);
            processPairs(pairArr, &t.dupIdx, &t.opticalDupIdx, &t.repIdx, &t.notDupIdx, &t.notOpticalDupIdx, &t.notRepIdx);
            itr = g.ckeyReadEndsMap.erase(itr);
        } else {
            ++itr;
        }
    }

    // 剩余的在lp中没处理的放到global里
    for (auto &ckVal : lp.ckeyReadEndsMap) {
        auto &pairArr = g.ckeyReadEndsMap[ckVal.first];
        pairArr.insert(pairArr.end(), ckVal.second.begin(), ckVal.second.end());

    }
    lp.ckeyReadEndsMap.clear();
    // 更新一下冗余结果
    refreshMarkDupData(t.dupIdx, t.opticalDupIdx, t.repIdx, t.notDupIdx, t.notOpticalDupIdx, t.notRepIdx, lp, cp);
    // 处理g中新找到的匹配
    putDupinfoToGlobal(g, lp);

    for (auto &unPair : lsm.unpairedDic) {
        g.unpairedDic.insert(unPair);
    }
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
        spdlog::info("{} reads processed in {} round", readNum, pipeArg.readOrder);

        pipeArg.readData.bams = inBamBuf.GetBamArr();

        pipeArg.readData.bamStartIdx = readNumSum;
        pipeArg.readData.taskSeq = pipeArg.readOrder;

        readNumSum += readNum;
        inBamBuf.ClearAll();                        // 清理上一轮读入的数据
        pipeArg.readOrder += 1;                     // for next
        yarn::TWIST(pipeArg.readSig, yarn::BY, 1);  // 读入了一轮数据
    }
    spdlog::info("total reads processed {}, last order {}", readNumSum, pipeArg.readOrder);
    return 0;
}
static void *pipeGenRE(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;
    auto &genREData = pipeArg.genREData;
    // init generate read ends data by num_thread
    int genREThread = pipeArg.numThread;
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

        if (pipeArg.readFinish) { // 读取结束的时候，没有新的数据需要处理了
            yarn::POSSESS(pipeArg.genRESig);
            pipeArg.genREFinish = 1;
            yarn::TWIST(pipeArg.genRESig, yarn::BY, 1);
            yarn::TWIST(pipeArg.readSig, yarn::BY, -1);
            break;
        }
        /* 处理bam，生成readends */
        PROF_START(gen);
        doGenRE(pipeArg);
        PROF_END(gprof[GP_gen], gen);

        yarn::POSSESS(pipeArg.genRESig);
        pipeArg.genREOrder += 1;
        yarn::TWIST(pipeArg.genRESig, yarn::BY, 1);
        yarn::TWIST(pipeArg.readSig, yarn::BY, -1);  // 使用了一次读入的数据
    }
    spdlog::info("genRE last order {}", pipeArg.genREOrder);
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
            cout << "zzh pipeSort: " << pipeArg.genREOrder << '\t' << pipeArg.sortOrder << endl;
            while (pipeArg.sortOrder < pipeArg.genREOrder) {
                yarn::POSSESS(pipeArg.sortSig);
                yarn::WAIT_FOR(pipeArg.sortSig, yarn::NOT_TO_BE, pipeArg.SORTBUFNUM);  // 有BUFNUM个位置
                yarn::RELEASE(pipeArg.sortSig);

                PROF_START(sort);
                doSort(pipeArg);
                PROF_END(gprof[GP_sort], sort);

                yarn::POSSESS(pipeArg.sortSig);
                pipeArg.sortOrder += 1;
                yarn::TWIST(pipeArg.sortSig, yarn::BY, 1);
            }
            yarn::POSSESS(pipeArg.sortSig);
            pipeArg.sortFinish = 1;
            yarn::TWIST(pipeArg.sortSig, yarn::BY, 1);
            break;
        }
        /* 排序 readends */
        PROF_START(sort);
        doSort(pipeArg);
        PROF_END(gprof[GP_sort], sort);

        yarn::POSSESS(pipeArg.genRESig);
        yarn::TWIST(pipeArg.genRESig, yarn::BY, -1);

        yarn::POSSESS(pipeArg.sortSig);
        pipeArg.sortOrder += 1;
        yarn::TWIST(pipeArg.sortSig, yarn::BY, 1);
    }
    spdlog::info("sort last order {}", pipeArg.sortOrder);
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
            cout << "zzh pipeMarkDup: " << pipeArg.sortOrder << '\t' << pipeArg.markDupOrder << endl;
            // 应该还得处理剩余的数据
            while (pipeArg.markDupOrder < pipeArg.sortOrder) {
                cout << "zzh pipeMarkDup: " << pipeArg.sortOrder << '\t' << pipeArg.markDupOrder << endl;
                yarn::POSSESS(pipeArg.markDupSig);
                yarn::WAIT_FOR(pipeArg.markDupSig, yarn::NOT_TO_BE, pipeArg.MARKBUFNUM);
                yarn::RELEASE(pipeArg.markDupSig);

                PROF_START(markdup);
                doMarkDup(pipeArg);
                PROF_END(gprof[GP_markdup], markdup);

                yarn::POSSESS(pipeArg.markDupSig);
                pipeArg.markDupOrder += 1;
                yarn::TWIST(pipeArg.markDupSig, yarn::BY, 1);
            }
            cout << "zzh pipeMarkDup: " << pipeArg.sortOrder << '\t' << pipeArg.markDupOrder << endl;
            yarn::POSSESS(pipeArg.markDupSig);
            pipeArg.markDupFinish = 1;
            yarn::TWIST(pipeArg.markDupSig, yarn::TO, 2 + pipeArg.MARKBUFNUM);
            break;
        }
        /* 冗余检测 readends */
        PROF_START(markdup);
        doMarkDup(pipeArg);
        PROF_END(gprof[GP_markdup], markdup);
        yarn::POSSESS(pipeArg.sortSig);
        yarn::TWIST(pipeArg.sortSig, yarn::BY, -1);

        yarn::POSSESS(pipeArg.markDupSig);
        pipeArg.markDupOrder += 1;
        yarn::TWIST(pipeArg.markDupSig, yarn::BY, 1);
    }
    spdlog::info("markdup last order {}", pipeArg.markDupOrder);
    return 0;
}
static void *pipeIntersect(void *data) {
    PipelineArg &pipeArg = *(PipelineArg *)data;
    const int kInitIntersectOrder = 1;
    pipeArg.intersectOrder = kInitIntersectOrder;
    while (1) {
        PROF_START(intersect_wait);
        yarn::POSSESS(pipeArg.markDupSig);
        yarn::WAIT_FOR(pipeArg.markDupSig, yarn::TO_BE_MORE_THAN, kInitIntersectOrder);  // 等待有数据
        yarn::RELEASE(pipeArg.markDupSig);
        PROF_END(gprof[GP_intersect_wait], intersect_wait);

        if (pipeArg.markDupFinish) {
            cout << "zzh pipeIntersect: " << pipeArg.markDupOrder << '\t' << pipeArg.intersectOrder << endl;
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
    spdlog::info("intersect last order {}", pipeArg.intersectOrder);
    return 0;
}

/* 当所有任务结束后，global data里还有未处理的数据 */
static void processLastData(PipelineArg &pipeArg) {
    IntersectData &g = pipeArg.intersectData;
    MarkDupData &lp = pipeArg.markDupData[(pipeArg.intersectOrder - 1) % pipeArg.MARKBUFNUM];
    SortMarkData &lsm = *(SortMarkData *)lp.dataPtr;
    int64_t lastLeft = INT64_MIN;
    if (lsm.pairs.size() > 0) {
        lastLeft = lsm.frags[0].Left();
    }
    // 处理找到匹配的global数据
    TaskSeqDupInfo t;
    for (auto itr = g.ckeyReadEndsMap.begin(); itr != g.ckeyReadEndsMap.end();) {
        auto &ckVal = *itr;
        auto &ck = ckVal.first;
        auto &pairArr = ckVal.second;
        if (ck.Left() >= lastLeft) {
            getEqualRE(pairArr[0], lsm.pairs, &pairArr);
        }
        sort(pairArr.begin(), pairArr.end(), ReadEnds::PairLittleThan );
        processPairs(pairArr, &t.dupIdx, &t.opticalDupIdx, &t.repIdx);
        itr = g.ckeyReadEndsMap.erase(itr);
    }
    // 更新一下冗余结果
    refreshMarkDupData(t.dupIdx, t.opticalDupIdx, t.repIdx, t.notDupIdx, t.notOpticalDupIdx, t.notRepIdx, lp);
    // 处理g中新找到的匹配
    putDupinfoToGlobal(g, lp);
}

static void parallelPipeline() {
    PipelineArg pipeArg(&nsgv::gDupRes);
    pipeArg.numThread = nsgv::gMdArg.NUM_THREADS;

    pthread_t tidArr[5];
    pthread_create(&tidArr[0], 0, pipeRead, &pipeArg);
    pthread_create(&tidArr[1], 0, pipeGenRE, &pipeArg);
    pthread_create(&tidArr[2], 0, pipeSort, &pipeArg);
    pthread_create(&tidArr[3], 0, pipeMarkDup, &pipeArg);
    pthread_create(&tidArr[4], 0, pipeIntersect, &pipeArg);
    for (int i = 0; i < 5; ++i) pthread_join(tidArr[i], 0);
    PROF_START(merge_result);
    processLastData(pipeArg);
    PROF_END(gprof[GP_merge_result], merge_result);

    spdlog::info("pipeArg size : {} GB", pipeArg.byteSize() / 1024.0 / 1024 / 1024);

    size_t repNum = 0;
    for (auto &v : pipeArg.intersectData.repIdxArr) repNum += v.size();
    spdlog::info("rep num : {}", repNum);

    spdlog::info("result size : {} GB", nsgv::gDupRes.byteSize() / 1024.0 / 1024 / 1024);
}

/* 并行流水线方式处理数据，标记冗余 */
void PipelineMarkDups() {
    if (nsgv::gMdArg.NUM_THREADS > 1)
        return parallelPipeline();

    PipelineArg pipeArg(&nsgv::gDupRes);
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
        spdlog::info("{} reads processed in {} round", readNum, pipeArg.readOrder);

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
    processLastData(pipeArg);
}