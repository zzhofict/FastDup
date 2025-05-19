/*
     Description: sam/bam，buf，

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
    bam
*/
/*
 * sam read
 */
struct BamWrap {
    // contigpos
    const static int MAX_CONTIG_LEN_SHIFT = 40;  // id，
    const static int READ_MAX_LENGTH = 200;
    const static int READ_MAX_DEPTH = 1000;  // ，

    // 
    bam1_t *b;
    int64_t end_pos_;  // bam, 

    // 
    inline int64_t start_pos() { return bam_global_pos(b); }
    // 
    inline int64_t end_pos() { return end_pos_; }
    // reference
    inline int16_t read_len() { return (end_pos_ - start_pos() + 1); }

    // contig
    inline int32_t contig_pos() { return b->core.pos; }
    // contig
    inline int32_t contig_end_pos() { return bam_pos(end_pos_); }
    // （AGTC）
    inline int16_t seq_len() { return b->core.l_qseq; }

    // softclip
    inline int32_t softclip_start() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        const char c = bam_cigar_opchr(cigar[0]);
        const int len = bam_cigar_oplen(cigar[0]);
        if (c == 'S')
            return bc.pos - len;
        return bc.pos;
    }

    // softclip
    inline int32_t softclip_end() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        const char c = bam_cigar_opchr(cigar[bc.n_cigar - 1]);
        const int len = bam_cigar_oplen(cigar[bc.n_cigar - 1]);
        if (c == 'S')
            return bam_pos(end_pos_) + len;
        return bam_pos(end_pos_);
    }

    // softclip
    inline int32_t right_softclip_len() {
        const uint32_t *cigar = bam_get_cigar(b);
        const bam1_core_t &bc = b->core;
        const char c = bam_cigar_opchr(cigar[bc.n_cigar - 1]);
        const int len = bam_cigar_oplen(cigar[bc.n_cigar - 1]);
        if (c == 'S')
            return len;
        return 0;
    }

    // 
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

    // 
    inline const char *query_name() { return bam_get_qname(b); }
    // cigar 
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

    // 
    inline int16_t length() { return sizeof(*this) + sizeof(bam1_t) + b->l_data; }

    // cigarinsert
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

    // cigardelete
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

    // sam read
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

    // bamadapterBoundary
    int GetAdapterBoundary() {
        const bam1_core_t &bc = b->core;
        int adapterBoundary;
        if (!HasWellDefinedFragmentSize())
            adapterBoundary = INT_MIN;
        else if (bc.flag & BAM_FREVERSE)
            adapterBoundary = bc.mpos - 1;
        else
            adapterBoundary = bc.pos + abs(bc.isize);  // GATK4.0  GATK3.5，3.5+1
        return adapterBoundary;
    }

    // I
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

    // soft clip(
    // IS?)
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

    // unclipped(hardclip)
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

    // unclipped(hardclip)
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

    /*  */
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

    /* flag */

    /*  unmapped */
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

    /*  */
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

    // bam，
    static inline int64_t bam_global_pos(bam1_t *b) {
        return (((int64_t)b->core.tid << MAX_CONTIG_LEN_SHIFT) | (int64_t)b->core.pos);
    }
    static inline int64_t bam_global_pos(int tid, int pos) {
        return (((int64_t)tid << MAX_CONTIG_LEN_SHIFT) | (int64_t)pos);
    }
    // bam
    static inline int32_t bam_tid(int64_t global_pos) {
        const int64_t mask = ~(((int64_t)1 << MAX_CONTIG_LEN_SHIFT) - 1);
        const int64_t high_tid = global_pos & mask;
        return (int32_t)(high_tid >> MAX_CONTIG_LEN_SHIFT);
    }
    // bam()
    static inline int32_t bam_pos(int64_t global_pos) {
        const int64_t mask = ((int64_t)1 << MAX_CONTIG_LEN_SHIFT) - 1;
        return (int32_t)(global_pos & mask);
    }

    // 
    void SetDuplicateReadFlag(bool flag) { setFlag(flag, BAM_FDUP); }

    void setFlag(bool flag, int bit) {
        if (flag)
            this->b->core.flag |= bit;
        else
            this->b->core.flag &= ~bit;
    }
};

typedef std::map<const std::string, std::vector<BamWrap *>> SampleBamMap;