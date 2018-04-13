#ifndef VRNA_PARALLELIZATION_HELPERS
#define VRNA_PARALLELIZATION_HELPERS

#if VRNA_WITH_PTHREADS

#include <pthread.h>
#include "thpool.h"

pthread_mutex_t output_mutex;
pthread_mutex_t output_file_mutex;

#define THREADSAFE_FILE_OUTPUT(a) { \
    pthread_mutex_lock(&output_file_mutex); \
    (a); \
    pthread_mutex_unlock(&output_file_mutex); \
}

#define THREADSAFE_STREAM_OUTPUT(a) { \
    pthread_mutex_lock(&output_mutex); \
    (a); \
    pthread_mutex_unlock(&output_mutex); \
}

#define INIT_PARALLELIZATION(a) \
  pthread_mutex_init(&output_mutex, NULL); \
  pthread_mutex_init(&output_file_mutex, NULL); \
  /* initialize thread pool */ \
  threadpool worker_pool = thpool_init(a);

#define UNINIT_PARALLELIZATION  \
  thpool_wait(worker_pool); \
  pthread_mutex_lock(&output_mutex); \
  pthread_mutex_unlock(&output_mutex); \
  pthread_mutex_destroy(&output_mutex); \
  pthread_mutex_lock(&output_file_mutex); \
  pthread_mutex_unlock(&output_file_mutex); \
  pthread_mutex_destroy(&output_file_mutex); \
  thpool_destroy(worker_pool);

#define RUN_IN_PARALLEL(fun, data)  { \
    thpool_add_work(worker_pool, (void *)&fun, (void *)data); \
}

#define WAIT_FOR_FREE_SLOT(a) { \
    while (thpool_num_threads_working(worker_pool) == (a)) \
      sleep(1); \
}

#else

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


#endif
