/*
     Description: 读入sam/bam时，开辟一个大的buf，存放这些数据

     Copyright : All right reserved by ICT

     Author : Zhang Zhonghai
     Date : 2019/11/27
*/

#include "bam_buf.h"

/*
 * BamBuf类
 */
// 读取数据直到读完，或者缓冲区满
int BamBuf::ReadBam() {
    int read_num = 0;
    if (handle_last) {             // 处理上次读入的最后一个bam
        if (has_enough_space()) {  // 必须调用，在边界处调整memffset
            ++read_num;
            append_one_bam();
        } else {
            return read_num;  // 还是没空间
        }
    }
    while (read_stat_ >= 0 && (read_stat_ = sam_read1(fp, hdr, bw->b)) >= 0) {
        bw->end_pos_ = BamWrap::BamEndPos(bw->b);
        if (has_enough_space()) {  // 还有空间
                                   // if (true) {  // 还有空间
            append_one_bam();
            ++read_num;  // 放进缓存才算读取到
        } else {
            break;
        }
    }
    if (read_stat_ >= 0) {
        handle_last = true;
    } else {
        handle_last = false;
    }
    return read_num;
}

// 初始化缓存
void BamBuf::Init(samFile *fp, sam_hdr_t *hdr, int64_t mem_size) {
    this->fp = fp;
    this->hdr = hdr;
    this->mem_size = mem_size;
    this->mem = (uint8_t *)malloc(mem_size);
    this->bw = (BamWrap *)malloc(sizeof(BamWrap));
    this->bw->b = bam_init1();
    if (bw == NULL || this->mem == NULL || this->bw->b == NULL) {
        fprintf(stderr, "allocate memory failed! Abort\n");
        exit(-1);
    }
}

void BamBuf::ClearBeforeIdx(size_t idxInBv) {
    if (idxInBv < 1)
        return;
    int i = 0, j = idxInBv;
    for (; j < bv.size(); ++i, ++j) {
        bv[i] = bv[j];
    }
    bv.resize(i);
    prepare_read();
}

void BamBuf::ClearAll() {
    bv.clear();
    prepare_read();
}

// 为下一次读取做准备, 计算一些边界条件
inline void BamBuf::prepare_read() {
    // 计算余留的下次计算可能用到的bam所占的位置
    if (bv.size() > 0) {
        BamWrap *bw = bv[0];
        legacy_start = (int64_t)bw - (int64_t)mem;
        bw = bv.back();
        legacy_end = (int64_t)bw + bw->length() - (int64_t)mem;
    } else {
        legacy_start = legacy_end = 0;
        mem_offset = 0;  // 上次没剩下，那就从头存储
    }
}

// 检查缓存是否还有空间
inline bool BamBuf::has_enough_space() {
    const uint32_t bam_len = bw->length();
    int64_t potential_end = mem_offset + bam_len;
    if (legacy_start <= legacy_end)
        legacy_start += mem_size;
    if (potential_end >= legacy_start) {
        return false;
    }
    if (potential_end >= mem_size) {
        mem_offset = 0;
    }
    int64_t virtual_offset = mem_offset;
    if (virtual_offset < legacy_end)
        virtual_offset += mem_size;
    potential_end = virtual_offset + bam_len;
    return potential_end < legacy_start;
}

// 处理一个读取后的bam
inline void BamBuf::append_one_bam() {
    BamWrap *bwp = (BamWrap *)(mem + mem_offset);
    *bwp = *bw;
    bwp->b = (bam1_t *)((char *)bwp + sizeof(*bwp));
    bam1_t *bp = bwp->b;
    *bp = *bw->b;
    bp->data = (uint8_t *)((char *)bwp->b + sizeof(bam1_t));
    memcpy(bp->data, bw->b->data, bw->b->l_data);
    // 更新下次存储的位置
    mem_offset = (mem_offset + bw->length() + 8 - 1) & ~((size_t)(8 - 1));
    bv.push_back(bwp);
}

// 处理上次读入的最后一个read
inline bool BamBuf::handle_last_read() {
    if (handle_last) {             // 处理上次读入的最后一个bam
        if (has_enough_space()) {  // 必须调用，在边界处调整memffset
            append_one_bam();
            handle_last = false;
            return true;
        }
    }
    return false;
}

/*
 * AsyncIoBamBuf 类
 */
// 初始化缓存
void AsyncIoBamBuf::Init(samFile *fp, sam_hdr_t *hdr, int64_t mem_size) {
    if (use_async_io_) {
        buf1_.Init(fp, hdr, mem_size >> 1);
        buf2_.Init(fp, hdr, mem_size >> 1);
        pi_ = &buf1_;
        po_ = &buf2_;
        tid_ = (pthread_t *)malloc(sizeof(pthread_t));
    } else {
        buf1_.Init(fp, hdr, mem_size);
        pi_ = &buf1_;
    }
}

// 读取数据
int AsyncIoBamBuf::ReadBam() {
    if (use_async_io_) {
        hasThread = true;
        return async_read_bam();
    } else {
        return sync_read_bam();
    }
}

int AsyncIoBamBuf::sync_read_bam() {
    int read_num = 0;
    if (clear_all_) {
        clear_all_ = false;
        pi_->ClearAll();
    } else if (clear_before_idx_ > 0) {
        pi_->ClearBeforeIdx(clear_before_idx_);
        clear_before_idx_ = 0;
    }
    read_num = pi_->ReadBam();
    refresh_bam_arr();
    return read_num;
}

int AsyncIoBamBuf::async_read_bam() {
    int read_num = 0;
    if (first_read_) {
        read_num = pi_->ReadBam();
        first_read_ = false;
        refresh_bam_arr();
    } else {
        // join, 交换缓冲区指针
        pthread_join(*tid_, 0);
        resize_buf();

        if (need_read_) {  // 需要交换指针
            BamBuf *tmp = pi_;
            pi_ = po_;
            po_ = tmp;
        }
        read_num = last_read_num_;
        refresh_bam_arr();
    }
    // 异步读
    pthread_create(tid_, 0, async_read, this);
    return read_num;
}

void *AsyncIoBamBuf::async_read(void *data) {
    AsyncIoBamBuf *ab = (AsyncIoBamBuf *)data;
    if (ab->need_read_ && ab->ReadStat() >= 0) {  // 需要读取
        ab->last_read_num_ = ab->po_->ReadBam();
    } else {
        ab->last_read_num_ = 0;
    }
    pthread_exit(0);
}

// 为下一次读取做准备,
// 计算一些边界条件，延迟操作，因为此时可能po_对应的buf正在读取
void AsyncIoBamBuf::ClearBeforeIdx(size_t idxInBv) { clear_before_idx_ = idxInBv; }

// 清空上一次所有读入的数据，延迟操作，因为此时可能po_对应的buf正在读取
void AsyncIoBamBuf::ClearAll() { clear_all_ = true; }

inline void AsyncIoBamBuf::resize_buf() {
    if (clear_all_) {  // 清理上一轮的数据
        clear_all_ = false;
        po_->ClearBeforeIdx(legacy_size_);
        pi_->ClearAll();
        if (pi_->handle_last_read()) {  // 上次读取有一个read没放入缓存
            last_read_num_ += 1;
            legacy_size_ = pi_->Size();  // 应该只有一个read
            need_read_ = true;
        } else {  // 没空间存放，则不交换指针，或者文件已经读取完毕
            legacy_size_ = 0;
            need_read_ = false;
        }
    } else if (clear_before_idx_ > 0) {
        if (clear_before_idx_ < legacy_size_) {
            po_->ClearBeforeIdx(clear_before_idx_);
            legacy_size_ -= clear_before_idx_;
            // 不需要交换指针，不需要读取
            need_read_ = false;
        } else {
            po_->ClearBeforeIdx(legacy_size_);
            pi_->ClearBeforeIdx(clear_before_idx_ - legacy_size_);
            if (pi_->handle_last_read()) {  // 上次读取有一个read没放入缓存
                last_read_num_ += 1;
                legacy_size_ = pi_->Size();  // 应该只有一个read
                need_read_ = true;
            } else {  // 没空间存放，则不交换指针，或者文件已经读取完毕
                legacy_size_ = 0;
                need_read_ = false;
            }
        }
        clear_before_idx_ = 0;
    }
}