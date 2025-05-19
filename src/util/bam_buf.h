/*
     Description: sam/bam buf

     Copyright : All right reserved by ICT

     Author : Zhang Zhonghai
     Date : 2019/11/27
*/
#pragma once

#include <htslib/sam.h>
#include <pthread.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <vector>

#include "bam_wrap.h"

using std::vector;
using namespace std;

/*
 * bam
 */
struct BamBuf {
    sam_hdr_t *hdr;            // samheader
    samFile *fp;               // sam
    BamWrap *bw = nullptr;     // bam
    uint8_t *mem = nullptr;    // bam,
                               // ，
    int64_t mem_offset = 0;    // 
    int64_t mem_size;          // 
    int read_stat_ = 0;        // ，
    vector<BamWrap *> bv;      // bam
    int64_t legacy_start = 0;  // bammem, 
    int64_t legacy_end = 0;    // bammem, 
    bool handle_last = false;  // bam

    // 
    void Init(samFile *fp, sam_hdr_t *hdr, int64_t mem_size);
    // ，
    int ReadBam();
    // , 
    void ClearBeforeIdx(size_t idxInBv);
    // 
    void ClearAll();
    inline int64_t Size() { return bv.size(); }  // read
    inline int ReadStat() { return read_stat_; }  // ，（）
    ~BamBuf() {
        if (this->mem != nullptr) {
            free(this->mem);
        }
        if (this->bw != nullptr) {
            bam_destroy1(bw->b);
            free(bw);
        }
    }
    void FreeMemory()  // 
    {
        if (this->mem != nullptr) {
            free(this->mem);
        }
        if (this->bw != nullptr) {
            bam_destroy1(bw->b);
            free(bw);
        }
        this->mem = nullptr;
        this->bw = nullptr;
    }
    void prepare_read();
    // 
    bool has_enough_space();
    // bam
    void append_one_bam();
    // read
    bool handle_last_read();

    // bv
    inline BamWrap *operator[](int64_t pos) { return bv[pos]; }
    inline void push_back(BamWrap *val) { bv.push_back(val); }
    inline void clear() { bv.clear(); }
    inline void resize(int64_t s) { bv.resize(s); }
};

/*
 * io
 */
struct AsyncIoBamBuf {
    BamBuf buf1_;
    BamBuf buf2_;
    BamBuf *pi_;  // buf
    BamBuf *po_;  // buf
    pthread_t *tid_ = NULL;
    bool hasThread = false;
    int64_t legacy_size_ = 0;  // ，read
    bool first_read_ = true;
    int last_read_num_ = 0;  // reads
    bool need_read_ = true;
    bool use_async_io_ = true;
    int64_t clear_before_idx_ = 0;  // ，clear_before_idx_reads
    bool clear_all_ = false;        // ，reads

    vector<BamWrap *> bam_arr_;  // bufbam

    AsyncIoBamBuf() {}
    AsyncIoBamBuf(bool use_async) : use_async_io_(use_async) {}
    // 
    ~AsyncIoBamBuf() {
        if (tid_ != NULL) {
            if (hasThread)
                pthread_join(*tid_, 0);
            free(tid_);
        }
        // 
        // buf
    }

    // 
    void Init(samFile *fp, sam_hdr_t *hdr, int64_t mem_size);

    // 
    int ReadBam();
    // , 
    void ClearBeforeIdx(size_t idxInBv);
    vector<BamWrap *> &GetBamArr() { return bam_arr_; }  // bam array
    // 
    void ClearAll();
    // read
    inline int64_t Size() { return legacy_size_ + pi_->Size(); }
    inline int ReadStat() { return pi_->read_stat_; }
    inline BamWrap *operator[](int64_t pos) { return bam_arr_[pos]; }
    // reads
    inline vector<BamWrap *> Slice(size_t startIdx, size_t endIdx) {
        if (endIdx > startIdx) {
            auto begItr = bam_arr_.begin();
            return std::move(vector<BamWrap *>(begItr + startIdx, begItr + endIdx));
        }
        return std::move(vector<BamWrap *>());
    }
    void FreeMemory() {
        buf1_.FreeMemory();
        buf2_.FreeMemory();
    }

    // 
    int sync_read_bam();
    // 
    int async_read_bam();
    // 
    static void *async_read(void *data);
    void resize_buf();
    inline void refresh_bam_arr() {
        bam_arr_.resize(this->Size());
        for (int i = 0; i < bam_arr_.size(); ++i) {
            if (i < legacy_size_)
                bam_arr_[i] = (*po_)[i];
            else
                bam_arr_[i] = (*pi_)[i - legacy_size_];
        }
    }
};

typedef AsyncIoBamBuf BamBufType;

typedef vector<BamWrap *> BamArray;