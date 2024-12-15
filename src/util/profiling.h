#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>

// #define SHOW_PERF

#ifdef __cplusplus
extern "C" {
#endif

#define LIM_THREAD 128
#define LIM_THREAD_PROF_TYPE 128
#define LIM_GLOBAL_PROF_TYPE 128
#define LIM_THREAD_DATA_TYPE 128
#define LIM_GLOBAL_DATA_TYPE 128

#ifdef SHOW_PERF
extern uint64_t proc_freq;
extern uint64_t tprof[LIM_THREAD_PROF_TYPE][LIM_THREAD];
extern uint64_t gprof[LIM_GLOBAL_PROF_TYPE];
#endif

#ifdef SHOW_PERF
#define PROF_START(tmp_time) uint64_t prof_tmp_##tmp_time = RealtimeMsec()
#define PROF_START_AGAIN(tmp_time) prof_tmp_##tmp_time = RealtimeMsec()
#define PROF_END(result, tmp_time) result += RealtimeMsec() - prof_tmp_##tmp_time
#define PROF_PRINT_START(tmp_time) uint64_t tmp_time = RealtimeMsec()
#define PROF_PRINT_END(tmp_time)          \
    tmp_time = RealtimeMsec() - tmp_time; \
    fprintf(stderr, "time %-15s:     %0.2lfs\n", #tmp_time, tmp_time * 1.0 / proc_freq)
#else
#define PROF_START(tmp_time)
#define PROF_END(result, tmp_time)
#define PROF_PRINT_START(tmp_time)
#define PROF_PRINT_END(tmp_time)
#endif

// GLOBAL
enum { GP_0 = 0, GP_1, GP_2, GP_3, GP_4, GP_5, GP_6, GP_7, GP_8, GP_9, GP_10 };
enum {
    GP_read_wait = 11,
    GP_gen_wait,
    GP_sort_wait,
    GP_markdup_wait,
    GP_intersect_wait,
    GP_read,
    GP_gen,
    GP_sort,
    GP_markdup,
    GP_intersect,
    GP_merge_result,
    GP_markdup_pair,
    GP_markdup_frag,
    GP_sort_pair,
    GP_sort_frag,
    GP_merge_match,
    GP_merge_markdup,
    GP_merge_update,
    GP_merge_add,
    GP_markdup_all,
    GP_final_read,
    GP_write
};
// THREAD
enum { TP_0 = 0, TP_1, TP_2, TP_3, TP_4, TP_5, TP_6, TP_7, TP_8, TP_9, TP_10 };
enum { TP_gen = 11, TP_sort, TP_sort_frag, TP_sort_pair};

uint64_t RealtimeMsec(void);

int DisplayProfiling(int);

#ifdef __cplusplus
}
#endif