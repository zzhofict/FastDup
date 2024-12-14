/*
     Description: 读入sam/bam时，开辟一个大的buf，存放这些数据

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
 * 存放读入的bam数据
 */
struct BamBuf {
    sam_hdr_t *hdr;            // sam文件的header信息
    samFile *fp;               // sam文件指针
    BamWrap *bw = nullptr;     // 用来循环读入bam
    uint8_t *mem = nullptr;    // 用来存放bam的数据,
                               // 程序结束后自动释放，所以没在析构函数里释放
    int64_t mem_offset = 0;    // 下一次要存放的位置
    int64_t mem_size;          // 缓存大小
    int read_stat_ = 0;        // 读取状态，是否读完
    vector<BamWrap *> bv;      // 方便对bam数据的访问
    int64_t legacy_start = 0;  // 没处理完的bam在mem中的起始位置, 闭区间
    int64_t legacy_end = 0;    // 没处理完的bam在mem中的结束位置, 开区间
    bool handle_last = false;  // 上次最后读入的bam是否需要处理

    // 初始化缓存
    void Init(samFile *fp, sam_hdr_t *hdr, int64_t mem_size);
    // 读取数据直到读完，或者缓冲区满
    int ReadBam();
    // 为下一次读取做准备, 计算一些边界条件
    void ClearBeforeIdx(size_t idxInBv);
    // 清空上一次所有读入的数据
    void ClearAll();
    inline int64_t Size() { return bv.size(); }  // 包含多少个read
    inline int ReadStat() { return read_stat_; }  // 文件的读取状态，是否可读（读取完全）
    ~BamBuf() {
        if (this->mem != nullptr) {
            free(this->mem);
        }
        if (this->bw != nullptr) {
            bam_destroy1(bw->b);
            free(bw);
        }
    }
    void FreeMemory()  // 释放开辟的内存
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
    // 检查缓存是否还有空间
    bool has_enough_space();
    // 处理一个读取后的bam
    void append_one_bam();
    // 处理上次读入的最后一个read
    bool handle_last_read();

    // 针对bv的操作
    inline BamWrap *operator[](int64_t pos) { return bv[pos]; }
    inline void push_back(BamWrap *val) { bv.push_back(val); }
    inline void clear() { bv.clear(); }
    inline void resize(int64_t s) { bv.resize(s); }
};

/*
 * io异步缓冲区
 */
struct AsyncIoBamBuf {
    BamBuf buf1_;
    BamBuf buf2_;
    BamBuf *pi_;  // 当前用的buf
    BamBuf *po_;  // 后台在读取的buf
    pthread_t *tid_ = NULL;
    bool hasThread = false;
    int64_t legacy_size_ = 0;  // 上一轮运算之后，缓存中还剩余的上次读取的read数量
    bool first_read_ = true;
    int last_read_num_ = 0;  // 上一次读取了多少reads
    bool need_read_ = true;
    bool use_async_io_ = true;
    int64_t clear_before_idx_ = 0;  // 用户异步读取，下一轮读取之前清理掉clear_before_idx_之前的所有reads
    bool clear_all_ = false;        // 用于异步读取，下一轮读取之前清理掉之前的所有reads

    vector<BamWrap *> bam_arr_;  // 用来访问buf中的bam

    AsyncIoBamBuf() {}
    AsyncIoBamBuf(bool use_async) : use_async_io_(use_async) {}
    // 析构
    ~AsyncIoBamBuf() {
        if (tid_ != NULL) {
            if (hasThread)
                pthread_join(*tid_, 0);
            free(tid_);
        }
        // 其他的内存就等程序结束自动释放
        // buf的析构函数会自动调用
    }

    // 初始化缓存
    void Init(samFile *fp, sam_hdr_t *hdr, int64_t mem_size);

    // 读取数据
    int ReadBam();
    // 为下一次读取做准备, 计算一些边界条件
    void ClearBeforeIdx(size_t idxInBv);
    vector<BamWrap *> &GetBamArr() { return bam_arr_; }  // 获取bam array
    // 清空上一次所有读入的数据
    void ClearAll();
    // 包含的read数量
    inline int64_t Size() { return legacy_size_ + pi_->Size(); }
    inline int ReadStat() { return pi_->read_stat_; }
    inline BamWrap *operator[](int64_t pos) { return bam_arr_[pos]; }
    // 获取某一段reads
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

    // 同步读取
    int sync_read_bam();
    // 异步读取
    int async_read_bam();
    // 异步读取线程函数
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