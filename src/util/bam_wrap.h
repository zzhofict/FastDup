/*
     Description: 读入sam/bam时，开辟一个大的buf，存放这些数据

     Copyright : All right reserved by ICT

     Author : Zhang Zhonghai
     Date : 2019/11/27
*/
#pragma once

#include <htslib/sam.h>
#include <limits.h>
#include <math.h>

#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

/*
    这里的成员函数命名有点混乱，特此说明，小写加下划线的函数命名，无论是静态函数，还是普通成员函数，更侧重说明
    这是类似bam的一个属性，而大写加驼峰命名的函数，更侧重说明这是通过计算得出的。
*/
/*
 * sam read的封装
 */
struct BamWrap {
    // 将contig左移后加上pos作为全局位置
    const static int MAX_CONTIG_LEN_SHIFT = 40;  // 将染色体id左移多少位，和位点拼合在一起
    const static int READ_MAX_LENGTH = 200;
    const static int READ_MAX_DEPTH = 1000;  // 这只是用来初始化空间用的，深度大于这个值也没关系

    // 成员变量尽量少，减少占用内存空间
    bam1_t *b;
    int64_t end_pos_;  // bam的全局结束位置, 闭区间

    // 全局开始位置
    inline int64_t start_pos() { return bam_global_pos(b); }
    // 全局结束位置
    inline int64_t end_pos() { return end_pos_; }
    // 和reference对应的序列长度
    inline int16_t read_len() { return (end_pos_ - start_pos() + 1); }

    // 在contig内的开始位置
    inline int32_t contig_pos() { return b->core.pos; }
    // 在contig内部的结束位置
    inline int32_t contig_end_pos() { return bam_pos(end_pos_); }
    // 序列的长度（AGTC字母个数）
    inline int16_t seq_len() { return b->core.l_qseq; }

    // 算上开头的softclip
    inline int32_t softclip_start() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        const char c = bam_cigar_opchr(cigar[0]);
        const int len = bam_cigar_oplen(cigar[0]);
        if (c == 'S')
            return bc.pos - len;
        return bc.pos;
    }

    // 算上结尾的softclip
    inline int32_t softclip_end() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        const char c = bam_cigar_opchr(cigar[bc.n_cigar - 1]);
        const int len = bam_cigar_oplen(cigar[bc.n_cigar - 1]);
        if (c == 'S')
            return bam_pos(end_pos_) + len;
        return bam_pos(end_pos_);
    }

    // 算上结尾的softclip
    inline int32_t right_softclip_len() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        const char c = bam_cigar_opchr(cigar[bc.n_cigar - 1]);
        const int len = bam_cigar_oplen(cigar[bc.n_cigar - 1]);
        if (c == 'S')
            return len;
        return 0;
    }

    // 获取序列
    inline std::string sequence() {
        ostringstream oss;
        char *seq = (char *)bam_get_seq(b);
        const bam1_core_t &bc = b->core;
        const char base_to_char[16] = {'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};
        for (int i = 0; i < bc.l_qseq; ++i) {
            char base = base_to_char[bam_seqi(seq, i)];
            oss << base;
        }
        return std::move(oss.str());
    }

    // 获取名字
    inline const char *query_name() { return bam_get_qname(b); }
    // 获取cigar 字符串
    inline string cigar_str() {
        ostringstream oss;
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            oss << len << c;
        }
        return std::move(oss.str());
    }

    // 占用的内存大小
    inline int16_t length() { return sizeof(*this) + sizeof(bam1_t) + b->l_data; }

    // 获取cigar中insert的总长度
    inline int32_t insert_cigar_len() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        int ret = 0;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'I')
                ret += len;
        }
        return ret;
    }

    // 获取cigar中delete的总长度
    inline int32_t del_cigar_len() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        int ret = 0;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'D')
                ret += len;
        }
        return ret;
    }

    // 计算sam read的终点位置
    static inline int64_t BamEndPos(const bam1_t *b) {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        int start_offset = -1;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'D' || c == 'N' || c == 'M' || c == '=' || c == 'X')
                start_offset += len;
        }
        return (((int64_t)b->core.tid << MAX_CONTIG_LEN_SHIFT) | (int64_t)(b->core.pos + start_offset));
    };

    bool HasWellDefinedFragmentSize() {
        const bam1_core_t &bc = b->core;
        bool hasWellDefinedFragmentSize = true;
        if (bc.isize == 0 || !(bc.flag & BAM_FPAIRED) || ((bc.flag & BAM_FUNMAP) || (bc.flag & BAM_FMUNMAP)) ||
            ((bool)(bc.flag & BAM_FREVERSE) == (bool)(bc.flag & BAM_FMREVERSE))) {
            hasWellDefinedFragmentSize = false;
        } else if (bc.flag & BAM_FREVERSE) {
            hasWellDefinedFragmentSize = contig_end_pos() > bc.mpos ? true : false;
        } else {
            hasWellDefinedFragmentSize = bc.pos <= bc.mpos + bc.isize ? true : false;
        }
        return hasWellDefinedFragmentSize;
    }

    // 计算bam的adapterBoundary
    int GetAdapterBoundary() {
        const bam1_core_t &bc = b->core;
        int adapterBoundary;
        if (!HasWellDefinedFragmentSize())
            adapterBoundary = INT_MIN;
        else if (bc.flag & BAM_FREVERSE)
            adapterBoundary = bc.mpos - 1;
        else
            adapterBoundary = bc.pos + abs(bc.isize);  // GATK4.0 和 GATK3.5不一样，3.5的这里+1
        return adapterBoundary;
    }

    // 获取开头的I的长度
    inline int GetHeadInsertLen() {
        int insLen = 0;
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'I') {
                insLen = len;
                break;
            } else if (c != 'H' && c != 'S')
                break;
        }
        return insLen;
    }

    // 获取soft clip开始位置(能处理H和S相连的情况，有这种情况么？,
    // 注意开头的I要当做S？)
    inline int64_t GetSoftStart() {
        int64_t softStart = b->core.pos;
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'S' || c == 'I')
                softStart -= len;
            else if (c != 'H')
                break;
        }
        return softStart;
    }

    // 获取unclipped开始位置(包括hardclip)
    inline int64_t GetUnclippedStart() {
        int64_t start = b->core.pos;
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'S' || c == 'H')
                start -= len;
            else
                break;
        }
        return start;
    }

    // 获取unclipped结束位置(包括hardclip)
    inline int64_t GetUnclippedEnd() {
        int64_t end_pos = bam_endpos(b);
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        for (int i = bc.n_cigar - 1; i >= 0; --i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            if (c == 'S' || c == 'H')
                end_pos += len;
            else
                break;
        }
        return end_pos - 1;
    }

    /* 获取碱基质量分数的加和 */
    /** Calculates a score for the read which is the sum of scores over Q15. */
    inline int GetSumOfBaseQualities() {
        int score = 0;
        uint8_t *qual = bam_get_qual(b);
        for (int i = 0; i < b->core.l_qseq; ++i) {
            if (qual[i] >= 15)
                score += qual[i];
        }

        return score;
    }

    /* 与flag相关的检测 */

    /* 没有比对上 unmapped */
    inline bool GetReadUnmappedFlag() { return b->core.flag & BAM_FUNMAP; }

    /* Template having multiple segments in sequencing  */
    inline bool GetReadPairedFlag() { return b->core.flag & BAM_FPAIRED; }

    /**
     * the read fails platform/vendor quality checks.
     */
    inline bool GetReadFailsVendorQualityCheckFlag() { return b->core.flag & BAM_FQCFAIL; }

    /**
     * the mate is unmapped.
     */
    bool GetMateUnmappedFlag() { return b->core.flag & BAM_FMUNMAP; }

    /**
     * @return whether the alignment is secondary (an alternative alignment of
     * the read).
     */
    bool IsSecondaryAlignment() { return b->core.flag & BAM_FSECONDARY; }

    /**
     * @return whether the alignment is supplementary (a split alignment such as
     * a chimeric alignment).
     */
    bool GetSupplementaryAlignmentFlag() { return b->core.flag & BAM_FSUPPLEMENTARY; }

    /*
     * Tests if this record is either a secondary and/or supplementary
     * alignment;
     */
    bool IsSecondaryOrSupplementary() { return IsSecondaryAlignment() || GetSupplementaryAlignmentFlag(); }

    /**
     * the read is the first read in a pair.
     */
    bool GetFirstOfPairFlag() { return b->core.flag & BAM_FREAD1; }

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    bool GetReadNegativeStrandFlag() { return b->core.flag & BAM_FREVERSE; }

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    bool GetMateNegativeStrandFlag() { return b->core.flag & BAM_FMREVERSE; }

    /* 其他的一些信息 */
    inline int GetReferenceLength() {
        int length = 0;
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        for (int i = 0; i < bc.n_cigar; ++i) {
            const char c = bam_cigar_opchr(cigar[i]);
            const int len = bam_cigar_oplen(cigar[i]);
            switch (c) {
            case 'M':
            case 'D':
            case 'N':
            case '=':
            case 'X':
                length += len;
                break;
            default:
                break;
            }
        }
        return length;
    }

    // 计算bam的全局位置，算上染色体序号和比对位置
    static inline int64_t bam_global_pos(bam1_t *b) {
        return (((int64_t)b->core.tid << MAX_CONTIG_LEN_SHIFT) | (int64_t)b->core.pos);
    }
    static inline int64_t bam_global_pos(int tid, int pos) {
        return (((int64_t)tid << MAX_CONTIG_LEN_SHIFT) | (int64_t)pos);
    }
    // 根据全局位置获取bam的染色体序号
    static inline int32_t bam_tid(int64_t global_pos) {
        const int64_t mask = ~(((int64_t)1 << MAX_CONTIG_LEN_SHIFT) - 1);
        const int64_t high_tid = global_pos & mask;
        return (int32_t)(high_tid >> MAX_CONTIG_LEN_SHIFT);
    }
    // 根据全局位置获取bam的比对位置(染色体内)
    static inline int32_t bam_pos(int64_t global_pos) {
        const int64_t mask = ((int64_t)1 << MAX_CONTIG_LEN_SHIFT) - 1;
        return (int32_t)(global_pos & mask);
    }

    // 设置是否冗余的标记
    void SetDuplicateReadFlag(bool flag) { setFlag(flag, BAM_FDUP); }

    void setFlag(bool flag, int bit) {
        if (flag)
            this->b->core.flag |= bit;
        else
            this->b->core.flag &= ~bit;
    }
};

typedef std::map<const std::string, std::vector<BamWrap *>> SampleBamMap;