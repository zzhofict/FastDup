/*
Description:
bam，bam，bam

Copyright : All right reserved by ICT

Author : Zhang Zhonghai
Date : 2023/10/23
*/
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <spdlog/spdlog.h>

#include <iomanip>
#include <vector>

#include "dup_metrics.h"
#include "fastdup_version.h"
#include "md_args.h"
#include "md_funcs.h"
#include "md_pipeline.h"
#include "read_name_parser.h"
#include "util/profiling.h"

#define BAM_BLOCK_SIZE 16L * 1024 * 1024

namespace nsgv {

MarkDupsArg gMdArg;                   // 
std::vector<ReadNameParser> gNameParsers;  // read name parser
samFile *gInBamFp;                    // bam
sam_hdr_t *gInBamHeader;              // bam
samFile *gOutBamFp;                   // , sambam
sam_hdr_t *gOutBamHeader;             // header
DuplicationMetrics gMetrics;          // 
DupResult gDupRes;
PipelineArg gPipe(&gDupRes);
};  // namespace nsgv

// 
struct ByteBuf {
    uint8_t *buf = nullptr;
    int size = 0;      // 
    int capacity = 0;  // 
};

/*
 * 
 */
static string getFileExtension(const string &filename) {
    auto last_dot = filename.find_last_of('.');
    if (last_dot == string::npos) {
        return "";
    }
    return filename.substr(last_dot + 1);
}

// entrance of mark duplicates
int MarkDuplicates() {
    PROF_START(whole_process);
    /* */
    nsgv::gNameParsers.resize(nsgv::gMdArg.NUM_THREADS);
    for (auto &parser : nsgv::gNameParsers)
        parser.SetReadNameRegex(nsgv::gMdArg.READ_NAME_REGEX);  // read nametile，x，y

    /* bam */
    nsgv::gInBamFp = sam_open_format(nsgv::gMdArg.INPUT_FILE.c_str(), "r", nullptr);
    if (!nsgv::gInBamFp) {
        spdlog::error("[{}] load sam/bam file failed.\n", __func__);
        return -1;
    }
    hts_set_opt(nsgv::gInBamFp, HTS_OPT_BLOCK_SIZE, BAM_BLOCK_SIZE);
    nsgv::gInBamHeader = sam_hdr_read(nsgv::gInBamFp);  // header
    // (libraryId)
    nsgv::gMetrics.LIBRARY = sam_hdr_line_name(nsgv::gInBamHeader, "RG", 0);

    /*  */
    htsThreadPool htsPoolRead = {NULL, 0};   // ，
    htsThreadPool htsPoolWrite = {NULL, 0};  // 
    htsPoolRead.pool = hts_tpool_init(nsgv::gMdArg.NUM_THREADS);
    htsPoolWrite.pool = hts_tpool_init(nsgv::gMdArg.NUM_THREADS);
    if (!htsPoolRead.pool || !htsPoolWrite.pool) {
        spdlog::error("[{}] failed to set up thread pool", __LINE__);
        sam_close(nsgv::gInBamFp);
        return -1;
    }
    hts_set_opt(nsgv::gInBamFp, HTS_OPT_THREAD_POOL, &htsPoolRead);

    /*  */
    PROF_START(markdup_all);
    PipelineMarkDups();
    PROF_END(gprof[GP_markdup_all], markdup_all);

    /* ,  */
    char modeout[12] = "wb";
    sam_open_mode(modeout + 1, nsgv::gMdArg.OUTPUT_FILE.c_str(), NULL);
    nsgv::gOutBamFp = sam_open(nsgv::gMdArg.OUTPUT_FILE.c_str(), modeout);
    if (!nsgv::gOutBamFp) {
        spdlog::error("[{}] create output sam/bam file failed.\n", __func__);
        return -1;
    }
    nsgv::gOutBamHeader = sam_hdr_dup(nsgv::gInBamHeader);
    // header
    sam_hdr_add_line(nsgv::gOutBamHeader, "PG", "ID", nsgv::gMdArg.PROGRAM_RECORD_ID.c_str(), "VN", FASTDUP_VERSION,
                     "CL", nsgv::gMdArg.CLI_STR.c_str(), NULL);

    // 
    hts_set_opt(nsgv::gOutBamFp, HTS_OPT_BLOCK_SIZE, BAM_BLOCK_SIZE);
    hts_set_opt(nsgv::gOutBamFp, HTS_OPT_THREAD_POOL, &htsPoolWrite);
    sam_close(nsgv::gInBamFp);  // bam
    nsgv::gInBamFp = sam_open_format(nsgv::gMdArg.INPUT_FILE.c_str(), "r", nullptr);
    if (!nsgv::gInBamFp) {
        spdlog::error("[{}] load sam/bam file failed.\n", __func__);
        return -1;
    }
    hts_set_opt(nsgv::gInBamFp, HTS_OPT_BLOCK_SIZE, BAM_BLOCK_SIZE);
    hts_set_opt(nsgv::gInBamFp, HTS_OPT_THREAD_POOL, &htsPoolRead);
    nsgv::gInBamHeader = sam_hdr_read(nsgv::gInBamFp);
    if (sam_hdr_write(nsgv::gOutBamFp, nsgv::gOutBamHeader) != 0) {
        spdlog::error("failed writing header to \"{}\"", nsgv::gMdArg.OUTPUT_FILE);
        sam_close(nsgv::gOutBamFp);
        sam_close(nsgv::gInBamFp);
        return -1;
    }
    // index
    string indexFn = nsgv::gMdArg.OUTPUT_FILE + ".bai";  // min_shift = 0 bai
    if ("sam" == getFileExtension(nsgv::gMdArg.OUTPUT_FILE) || !nsgv::gMdArg.CREATE_INDEX) {
        indexFn = "";
    }
    if (!indexFn.empty()) {
        int index_min_shift = 0;
        if (nsgv::gMdArg.INDEX_FORMAT == nsmd::IndexFormat::CSI) {
            indexFn = nsgv::gMdArg.OUTPUT_FILE + ".csi";
            index_min_shift = 14;
        }
        if (sam_idx_init(nsgv::gOutBamFp, nsgv::gOutBamHeader, 0 /*csi 14*/, indexFn.c_str()) < 0) {
            spdlog::error("failed to open index \"{}\" for writing", indexFn);
            sam_close(nsgv::gOutBamFp);
            sam_close(nsgv::gInBamFp);
            return -1;
        }
    }

    // ,
    BamBufType inBuf(nsgv::gMdArg.DUPLEX_IO);
    inBuf.Init(nsgv::gInBamFp, nsgv::gInBamHeader, nsgv::gMdArg.MAX_MEM);

    DupIdxQueue<DupInfo> dupIdxQue, repIdxQue;
    DupIdxQueue<int64_t> opticalIdxQue;
    dupIdxQue.Init(&nsgv::gDupRes.dupIdxArr);
    repIdxQue.Init(&nsgv::gDupRes.repIdxArr);
    opticalIdxQue.Init(&nsgv::gDupRes.opticalDupIdxArr);
    spdlog::info("{} duplicate reads found", dupIdxQue.Size());
    spdlog::info("{} optical reads found", opticalIdxQue.Size());
    // spdlog::info("{} represent reads found", repIdxQue.Size());
    // dupIdxQue.RealSize("na12878.dup");
    // opticalIdxQue.RealSize("normal.odup");

    // return 0;

    uint64_t bamIdx = 0;
    DupInfo dupIdx = dupIdxQue.Pop();
    DupInfo repIdx = repIdxQue.Pop();
    uint64_t opticalIdx = opticalIdxQue.Pop();

    ByteBuf bytes;
    bam1_t *bp = bam_init1();
    bool isDup = false;
    bool isOpticalDup = false;
    bool isInDuplicateSet = false;
    uint32_t representativeReadIndexInFile = 0;
    uint32_t duplicateSetSize = 0;

    int64_t realDupSize = 0;

    PROF_START(write);
    while (inBuf.ReadStat() >= 0) {
        PROF_START(final_read);
        size_t readNum = inBuf.ReadBam();
        PROF_END(gprof[GP_final_read], final_read);
        // PROF_PRINT_START(read_write);
        for (size_t i = 0; i < readNum; ++i) {
            BamWrap *bw = inBuf[i];
            if (bam_copy1(bp, bw->b) == nullptr) {
                spdlog::error("Can not copy sam record!");
                return -1;
            }
            bw->b = bp;
            isDup = false;
            isOpticalDup = false;
            isInDuplicateSet = false;

            // duplicate tag
            if (nsgv::gMdArg.CLEAR_DT) {
                uint8_t *oldTagVal = bam_aux_get(bw->b, nsgv::gMdArg.DUPLICATE_TYPE_TAG.c_str());
                if (oldTagVal != NULL) bam_aux_del(bw->b, oldTagVal);
            }

            // 
            if (bw->GetReadUnmappedFlag()) {
                ++nsgv::gMetrics.UNMAPPED_READS;
            } else if (bw->IsSecondaryOrSupplementary()) {
                ++nsgv::gMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS;
            } else if (!bw->GetReadPairedFlag() || bw->GetMateUnmappedFlag()) {
                ++nsgv::gMetrics.UNPAIRED_READS_EXAMINED;
            } else {
                ++nsgv::gMetrics.READ_PAIRS_EXAMINED;  // will need to be divided by 2 at the end
            }

            /*  */
            if (bamIdx == dupIdx) {
                ++realDupSize;  // for test
                isDup = true;
                if (nsgv::gMdArg.TAG_DUPLICATE_SET_MEMBERS && dupIdx.dupSet != 0) {
                    isInDuplicateSet = true;
                    representativeReadIndexInFile = dupIdx.GetRepIdx();
                    duplicateSetSize = dupIdx.dupSet;
                }

                // ，dupidx，duprepIdxdupSetSize
                while ((dupIdx = dupIdxQue.Pop()) == bamIdx);
                if (opticalIdx == bamIdx)
                    isOpticalDup = true;
                else if (opticalIdx < bamIdx) {
                    while ((opticalIdx = opticalIdxQue.Pop()) < bamIdx);
                    if (opticalIdx == bamIdx)
                        isOpticalDup = true;
                }

                // 
                bw->SetDuplicateReadFlag(true);

                uint8_t *oldTagVal = bam_aux_get(bw->b, nsgv::gMdArg.DUPLICATE_TYPE_TAG.c_str());
                if (oldTagVal != NULL) bam_aux_del(bw->b, oldTagVal);

                if (isOpticalDup)
                    bam_aux_append(bw->b, nsgv::gMdArg.DUPLICATE_TYPE_TAG.c_str(), 'Z',
                                   nsgv::gMdArg.DUPLICATE_TYPE_SEQUENCING.size() + 1,
                                   (const uint8_t *)nsgv::gMdArg.DUPLICATE_TYPE_SEQUENCING.c_str());
                else
                    bam_aux_append(bw->b, nsgv::gMdArg.DUPLICATE_TYPE_TAG.c_str(), 'Z',
                                   nsgv::gMdArg.DUPLICATE_TYPE_LIBRARY.size() + 1,
                                   (const uint8_t *)nsgv::gMdArg.DUPLICATE_TYPE_LIBRARY.c_str());

                // 
                if (!bw->IsSecondaryOrSupplementary() && !bw->GetReadUnmappedFlag()) {
                    // Update the duplication metrics
                    if (!bw->GetReadPairedFlag() || bw->GetMateUnmappedFlag()) {
                        ++nsgv::gMetrics.UNPAIRED_READ_DUPLICATES;
                    } else {
                        ++nsgv::gMetrics.READ_PAIR_DUPLICATES;  // will need to be divided by 2 at the end
                    }
                }
            } else {
                bw->SetDuplicateReadFlag(false);
            }
            if (nsgv::gMdArg.TAG_DUPLICATE_SET_MEMBERS && bamIdx == repIdx) {  // repressent
                isInDuplicateSet = true;
                representativeReadIndexInFile = repIdx.GetRepIdx();
                duplicateSetSize = repIdx.dupSet;
                while (repIdx == bamIdx || repIdx.dupSet == 0) {
                    if (repIdxQue.Size() > 0)
                        repIdx = repIdxQue.Pop();
                    else {
                        repIdx = -1;
                        break;
                    }
                }
            }

            if (nsgv::gMdArg.TAG_DUPLICATE_SET_MEMBERS && isInDuplicateSet) {
                if (!bw->IsSecondaryOrSupplementary() && !bw->GetReadUnmappedFlag()) {
                    uint8_t *oldTagVal = bam_aux_get(bw->b, nsgv::gMdArg.DUPLICATE_SET_INDEX_TAG.c_str());
                    if (oldTagVal != NULL)
                        bam_aux_del(bw->b, oldTagVal);
                    bam_aux_append(bw->b, nsgv::gMdArg.DUPLICATE_SET_INDEX_TAG.c_str(), 'i',
                                   sizeof(representativeReadIndexInFile), (uint8_t *)&representativeReadIndexInFile);
                    oldTagVal = bam_aux_get(bw->b, nsgv::gMdArg.DUPLICATE_SET_SIZE_TAG.c_str());
                    if (oldTagVal != NULL)
                        bam_aux_del(bw->b, oldTagVal);
                    bam_aux_append(bw->b, nsgv::gMdArg.DUPLICATE_SET_SIZE_TAG.c_str(), 'i', sizeof(duplicateSetSize),
                                   (const uint8_t *)&duplicateSetSize);
                }
            }
            // readoutput，，record
            ++bamIdx;
            if (isDup && nsgv::gMdArg.REMOVE_DUPLICATES) {
                continue;
            }
            if (isOpticalDup && nsgv::gMdArg.REMOVE_SEQUENCING_DUPLICATES) {
                continue;
            }
            if (!nsgv::gMdArg.PROGRAM_RECORD_ID.empty() && nsgv::gMdArg.ADD_PG_TAG_TO_READS) {
                uint8_t *oldTagVal = bam_aux_get(bw->b, "PG");
                if (oldTagVal != NULL)
                    bam_aux_del(bw->b, oldTagVal);
                bam_aux_append(bw->b, "PG", 'Z', nsgv::gMdArg.PROGRAM_RECORD_ID.size() + 1,
                               (const uint8_t *)nsgv::gMdArg.PROGRAM_RECORD_ID.c_str());
            }
#if 1
            if (sam_write1(nsgv::gOutBamFp, nsgv::gOutBamHeader, bw->b) < 0) {
                spdlog::error("failed writing sam record to \"{}\"", nsgv::gMdArg.OUTPUT_FILE.c_str());
                sam_close(nsgv::gOutBamFp);
                sam_close(nsgv::gInBamFp);
                return -1;
            }
#endif
        }
        inBuf.ClearAll();
        // PROF_PRINT_END(read_write);
        spdlog::info("write {} reads to output", readNum);
    }
    bam_destroy1(bp);

    // 
    nsgv::gMetrics.READ_PAIRS_EXAMINED /= 2;
    nsgv::gMetrics.READ_PAIR_DUPLICATES /= 2;
    for (auto &arr : nsgv::gDupRes.opticalDupIdxArr) nsgv::gMetrics.READ_PAIR_OPTICAL_DUPLICATES += arr.size();
    nsgv::gMetrics.READ_PAIR_OPTICAL_DUPLICATES = nsgv::gMetrics.READ_PAIR_OPTICAL_DUPLICATES / 2;
    nsgv::gMetrics.ESTIMATED_LIBRARY_SIZE =
        estimateLibrarySize(nsgv::gMetrics.READ_PAIRS_EXAMINED - nsgv::gMetrics.READ_PAIR_OPTICAL_DUPLICATES,
                            nsgv::gMetrics.READ_PAIRS_EXAMINED - nsgv::gMetrics.READ_PAIR_DUPLICATES);
    if (nsgv::gMetrics.UNPAIRED_READS_EXAMINED + nsgv::gMetrics.READ_PAIRS_EXAMINED != 0) {
        nsgv::gMetrics.PERCENT_DUPLICATION =
            (nsgv::gMetrics.UNPAIRED_READ_DUPLICATES + nsgv::gMetrics.READ_PAIR_DUPLICATES * 2) /
            (double)(nsgv::gMetrics.UNPAIRED_READS_EXAMINED + nsgv::gMetrics.READ_PAIRS_EXAMINED * 2);
    } else {
        nsgv::gMetrics.PERCENT_DUPLICATION = 0;
    }
    calculateRoiHistogram(nsgv::gMetrics);

    // 
    if (!nsgv::gMdArg.METRICS_FILE.empty()) {
        ofstream ofsM(nsgv::gMdArg.METRICS_FILE);
        ofsM << "## StringHeader" << endl;
        ofsM << "# " << nsgv::gMdArg.CLI_STR << endl;
        ofsM << "## StringHeader" << endl;
        ofsM << "# Started on: " << nsgv::gMdArg.START_TIME << endl << endl;
        ofsM << "## METRICS"
             << endl;
        ofsM << "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_"
                "READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_"
                "DUPLICATION\tESTIMATED_LIBRARY_SIZE"
             << endl;
        ofsM << nsgv::gMetrics.LIBRARY << "\t" << nsgv::gMetrics.UNPAIRED_READS_EXAMINED << "\t" << nsgv::gMetrics.READ_PAIRS_EXAMINED
             << "\t" << nsgv::gMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS << "\t" << nsgv::gMetrics.UNMAPPED_READS << "\t"
             << nsgv::gMetrics.UNPAIRED_READ_DUPLICATES << "\t" << nsgv::gMetrics.READ_PAIR_DUPLICATES << "\t"
             << nsgv::gMetrics.READ_PAIR_OPTICAL_DUPLICATES << "\t" << nsgv::gMetrics.PERCENT_DUPLICATION << "\t"
             << nsgv::gMetrics.ESTIMATED_LIBRARY_SIZE << endl
             << endl;
        ofsM << "## HISTOGRAM\tDouble" << endl;
        ofsM << "BIN CoverageMult" << endl;
        for (int i = 1; i <= 100; ++i) {
            ofsM << i << "\t" << std::fixed << std::setprecision(6) << nsgv::gMetrics.CoverageMult[i] << endl;
        }
        ofsM.close();
    }
    PROF_END(gprof[GP_write], write);

    if (!indexFn.empty() && sam_idx_save(nsgv::gOutBamFp) < 0) {
        spdlog::error("writing index failed");
        sam_close(nsgv::gOutBamFp);
        sam_close(nsgv::gInBamFp);
        return -1;
    }

    /* , */
    sam_close(nsgv::gOutBamFp);
    sam_close(nsgv::gInBamFp);

    PROF_END(gprof[GP_whole_process], whole_process);

    return 0;
}
