#include "profiling.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#ifdef SHOW_PERF
uint64_t tprof[LIM_THREAD_PROF_TYPE][LIM_THREAD] = {0};
uint64_t proc_freq = 1000;
uint64_t gprof[LIM_GLOBAL_PROF_TYPE] = {0};
#endif

uint64_t RealtimeMsec(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (uint64_t)1000 * (tv.tv_sec + ((1e-6) * tv.tv_usec));
}

static int CalcThreadTime(uint64_t *a, int len, double *max, double *min, double *avg) {
#ifdef SHOW_PERF
    int i = 0;
    uint64_t umax = 0, umin = UINT64_MAX, uavg = 0;
    for (i = 0; i < len; i++) {
        if (a[i] > umax)
            umax = a[i];
        if (a[i] < umin)
            umin = a[i];
        uavg += a[i];
    }
    *avg = uavg * 1.0 / len / proc_freq;
    *max = umax * 1.0 / proc_freq;
    *min = umin * 1.0 / proc_freq;
#endif
    return 1;
}

#define PRINT_GP(gpname) \
    fprintf(stderr, "time G %-15s:     %0.2lfs\n", #gpname, gprof[GP_##gpname] * 1.0 / proc_freq);

#define PRINT_TP(tpname, nthread)                                                                                 \
    {                                                                                                             \
        double maxTime, minTime, avgTime;                                                                         \
        CalcThreadTime(tprof[TP_##tpname], nthread, &maxTime, &minTime, &avgTime);                                \
        fprintf(stderr, "time T %-15s:  avg %0.2lfs min %0.2lfs max %0.2lfs\n", #tpname, avgTime, minTime, maxTime); \
    }

int DisplayProfiling(int nthread) {

#ifdef SHOW_PERF
    fprintf(stderr, "\n");
    PRINT_GP(read_wait);
    PRINT_GP(gen_wait);
    PRINT_GP(sort_wait);
    PRINT_GP(markdup_wait);
    PRINT_GP(intersect_wait);
    PRINT_GP(read);
    PRINT_GP(gen);
    PRINT_GP(sort);
    PRINT_GP(markdup);
    PRINT_GP(intersect);
    PRINT_GP(merge_result);

    PRINT_TP(gen, nthread);
    PRINT_TP(sort, nthread);
    fprintf(stderr, "\n");
#endif

    return 0;
}