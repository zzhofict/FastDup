/*
Description:
标记bam文件中的冗余信息，只处理按照坐标排序后的bam，且bam为单一样本数据

Copyright : All right reserved by ICT

Author : Zhang Zhonghai
Date : 2023/10/23
*/
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <spdlog/spdlog.h>

#include <vector>

#include "dup_metrics.h"
#include "fastdup_version.h"
#include "md_args.h"
#include "md_funcs.h"
#include "pipeline_md.h"
#include "read_name_parser.h"
#include "util/profiling.h"

#define BAM_BLOCK_SIZE 16L * 1024 * 1024

namespace nsgv {

MarkDupsArg gMdArg;                   // 用来测试性能
std::vector<ReadNameParser> gNameParsers;  // 每个线程一个read name parser
samFile *gInBamFp;                    // 输入的bam文件
sam_hdr_t *gInBamHeader;              // 输入的bam文件头信息
samFile *gOutBamFp;                   // 输出文件, sam或者bam格式
sam_hdr_t *gOutBamHeader;             // 输出文件的header
DuplicationMetrics gMetrics;          // 统计信息
PipelineArg gPipe;

};

/*
 * 获取文件名后缀
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

    /* 初始化一些参数和变量*/
    nsgv::gNameParsers.resize(nsgv::gMdArg.NUM_THREADS);
    for (auto &parser : nsgv::gNameParsers)
        parser.SetReadNameRegex(nsgv::gMdArg.READ_NAME_REGEX);  // 用来解析read name中的tile，x，y信息

    /* 打开输入bam文件 */
    nsgv::gInBamFp = sam_open_format(nsgv::gMdArg.INPUT_FILE.c_str(), "r", nullptr);
    if (!nsgv::gInBamFp) {
        spdlog::error("[{}] load sam/bam file failed.\n", __func__);
        return -1;
    }
    hts_set_opt(nsgv::gInBamFp, HTS_OPT_BLOCK_SIZE, BAM_BLOCK_SIZE);
    nsgv::gInBamHeader = sam_hdr_read(nsgv::gInBamFp);  // 读取header
    // 获取样本名称(libraryId)
    nsgv::gMetrics.LIBRARY = sam_hdr_line_name(nsgv::gInBamHeader, "RG", 0);

    /* 利用线程池对输入输出文件进行读写 */
    htsThreadPool htsPoolRead = {NULL, 0};   // 多线程读取，创建线程池
    htsThreadPool htsPoolWrite = {NULL, 0};  // 读写用不同的线程池
    htsPoolRead.pool = hts_tpool_init(nsgv::gMdArg.NUM_THREADS);
    htsPoolWrite.pool = hts_tpool_init(nsgv::gMdArg.NUM_THREADS);
    // htsPoolRead.pool = hts_tpool_init(8);
    // htsPoolWrite.pool = hts_tpool_init(32);
    if (!htsPoolRead.pool || !htsPoolWrite.pool) {
        spdlog::error("[{}] failed to set up thread pool", __LINE__);
        sam_close(nsgv::gInBamFp);
        return -1;
    }
    hts_set_opt(nsgv::gInBamFp, HTS_OPT_THREAD_POOL, &htsPoolRead);

    // 测试读写速度
#if 0
    bam1_t *bamp = bam_init1();
    while (sam_read1(nsgv::gInBamFp, nsgv::gInBamHeader, bamp) >= 0);
    DisplayProfiling(nsgv::gMdArg.NUM_THREADS);
    exit(0);
#endif

    /* 冗余检查和标记 */
    pipelineMarkDups();

    /* 初始化输出文件 */
    char modeout[12] = "wb";
    sam_open_mode(modeout + 1, nsgv::gMdArg.OUTPUT_FILE.c_str(), NULL);
    nsgv::gOutBamFp = sam_open(nsgv::gMdArg.OUTPUT_FILE.c_str(), modeout);
    nsgv::gOutBamHeader = sam_hdr_dup(nsgv::gInBamHeader);
    // 修改输出文件的header
    // sam_hdr_add_line(nsgv::gOutBamHeader, "PG", "ID", nsgv::gMdArg.PROGRAM_RECORD_ID.c_str(), "VN", FASTDUP_VERSION,
    //                  "CL",
    //                  (nsgv::gMdArg.PROGRAM_RECORD_ID + " " + nsgv::gMdArg.GetArgValueString() + " " +
    //                   nsgv::gMdArg.GetArgValueString())
    //                      .c_str(),
    //                  NULL);

    // 用同样的线程池处理输出文件
    hts_set_opt(nsgv::gOutBamFp, HTS_OPT_BLOCK_SIZE, BAM_BLOCK_SIZE);
    hts_set_opt(nsgv::gOutBamFp, HTS_OPT_THREAD_POOL, &htsPoolWrite);


    return 0;
}