#ifndef VRNA_PARALLELIZATION_HELPERS
#define VRNA_PARALLELIZATION_HELPERS

#if VRNA_WITH_PTHREADS

#include <pthread.h>
#include "thpool.h"

pthread_mutex_t output_mutex;
pthread_mutex_t output_file_mutex;
unsigned int    max_threads;
threadpool      worker_pool;

#define ATOMIC_BLOCK(a) { \
    if (max_threads > 1) { \
      pthread_mutex_lock(&output_file_mutex); \
      (a); \
      pthread_mutex_unlock(&output_file_mutex); \
    } else { \
      (a); \
    } \
}

#define THREADSAFE_FILE_OUTPUT(a) { \
    if (max_threads > 1) { \
      pthread_mutex_lock(&output_file_mutex); \
      (a); \
      pthread_mutex_unlock(&output_file_mutex); \
    } else { \
      (a); \
    } \
}

#define THREADSAFE_STREAM_OUTPUT(a) { \
    if (max_threads > 1) { \
      pthread_mutex_lock(&output_mutex); \
      (a); \
      pthread_mutex_unlock(&output_mutex); \
    } else { \
      (a); \
    } \
}

#define INIT_PARALLELIZATION(a) { \
    max_threads = ((a) > 1) ? (unsigned int)(a) : 1; \
    /* initialize semaphores and thread pool */ \
    if (max_threads > 1) { \
      pthread_mutex_init(&output_mutex, NULL); \
      pthread_mutex_init(&output_file_mutex, NULL); \
      worker_pool = thpool_init(max_threads); \
    } \
}

#define UNINIT_PARALLELIZATION  { \
    if (max_threads > 1) \
      thpool_wait(worker_pool); \
    pthread_mutex_destroy(&output_mutex); \
    pthread_mutex_destroy(&output_file_mutex); \
    if (max_threads > 1) \
      thpool_destroy(worker_pool); \
}

#define RUN_IN_PARALLEL(fun, data)  { \
    if (max_threads > 1) { thpool_add_work(worker_pool, (void *)&fun, (void *)data); } \
    else { fun(data); } \
}

#define WAIT_FOR_FREE_SLOT(a) { \
    if (max_threads > 1) { \
      while (thpool_num_threads_working(worker_pool) == (a)) \
        sleep(1); \
    } \
}

#else

#define ATOMIC_BLOCK(a)             { (a); }
#define THREADSAFE_FILE_OUTPUT(a)   { (a); }
#define THREADSAFE_STREAM_OUTPUT(a)   { (a); }
#define INIT_PARALLELIZATION(a)
#define UNINIT_PARALLELIZATION
#define RUN_IN_PARALLEL(fun, data)  { fun(data); }
#define WAIT_FOR_FREE_SLOT(a)

#endif

int
num_proc_cores(int  *num_cores,
               int  *num_cores_conf);


int
max_user_threads(void);


#endif
