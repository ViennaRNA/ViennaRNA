#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if VRNA_WITH_PTHREADS
# include <pthread.h>
#endif

#include "ViennaRNA/utils.h"
#include "ViennaRNA/stream_output.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define QUEUE_OVERHEAD  32


struct vrna_ordered_stream_s {
  unsigned int                start;
  unsigned int                end;
  unsigned int                size;
  unsigned int                shift;
  unsigned int                cancel;

  vrna_callback_stream_output *output;
  void                        **data;
  void                        *auxdata;
#if VRNA_WITH_PTHREADS
  pthread_mutex_t             mtx;
  pthread_cond_t              flush_cond;
  pthread_t                   flush_thread;
#endif
};

PRIVATE INLINE void
flush_output(struct vrna_ordered_stream_s *queue)
{
  unsigned int j;

  /* process all consecutive blocks available from the start */
  for (j = queue->start; (j <= queue->end) && (queue->data[j]); j++) {
    /* process output callback */
    if (queue->output)
      queue->output(queue->auxdata, j, queue->data[j]);

    queue->data[j] = NULL;
    queue->start++;
  }

  if (queue->end < queue->start) {
    queue->data[j]  = NULL;
    queue->end      = queue->start;
  }
}


#if VRNA_WITH_PTHREADS

PRIVATE void *
flush_thread(void *data)
{
  struct vrna_ordered_stream_s *queue;

  queue = (struct vrna_ordered_stream_s *)data;

  pthread_mutex_lock(&queue->mtx);

  do {
    pthread_cond_wait(&queue->flush_cond, &queue->mtx);
    flush_output(queue);
  } while (!queue->cancel);

  pthread_mutex_unlock(&queue->mtx);
  pthread_exit(NULL);
}


#endif


struct vrna_ordered_stream_s *
vrna_ostream_init(vrna_callback_stream_output *output,
                  void                        *auxdata)
{
  struct vrna_ordered_stream_s *queue;

  queue = (struct vrna_ordered_stream_s *)vrna_alloc(sizeof(struct vrna_ordered_stream_s));

  queue->start    = 0;
  queue->end      = 0;
  queue->size     = QUEUE_OVERHEAD;
  queue->shift    = 0;
  queue->output   = output;
  queue->auxdata  = auxdata;
  queue->cancel   = 0;
  queue->data     = (void **)vrna_alloc(sizeof(void *) * QUEUE_OVERHEAD);

#if VRNA_WITH_PTHREADS
  pthread_attr_t attr;
  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* init semaphores */
  pthread_mutex_init(&queue->mtx, NULL);
  pthread_cond_init(&queue->flush_cond, NULL);
  /* create flush thread */
  pthread_create(&queue->flush_thread, &attr, flush_thread, (void *)queue);
  pthread_attr_destroy(&attr);
#endif
  return queue;
}


void
vrna_ostream_free(struct vrna_ordered_stream_s *queue)
{
  if (queue) {
#if VRNA_WITH_PTHREADS
    void *status;

    /* flush buffer and signal flush thread to cancel */
    pthread_mutex_lock(&queue->mtx);
    queue->cancel = 1;
    pthread_cond_signal(&queue->flush_cond);
    pthread_mutex_unlock(&queue->mtx);

    /* wait for the thread to actually being cancelled */
    int rc = pthread_join(queue->flush_thread, &status);

    /* destroy semaphores */
    pthread_cond_destroy(&queue->flush_cond);
    pthread_mutex_destroy(&queue->mtx);
#endif
    /* free remaining memory */
    queue->data += queue->shift;
    free(queue->data);
    /* free ostream itself */
    free(queue);
  }
}


void
vrna_ostream_request(struct vrna_ordered_stream_s *queue,
                     unsigned int                 num)
{
  unsigned int i;

  if (queue) {
    if (num >= queue->end) {
#if VRNA_WITH_PTHREADS
      pthread_mutex_lock(&queue->mtx);
#elif defined (_OPENMP)
#     pragma omp critical (vrna_ostream_access)
      {
#endif
      /* check whether we have to increase memory */
      unsigned int new_size = num - queue->shift + 1;

      if (queue->size < new_size) {
        /*
         *  Check whether we can simply move data around to get more space.
         *  We do this only, if more than half of the first buffer is empty.
         *  Otherwise, we simply increase the buffer to fit incoming data.
         */
        if ((queue->start - queue->shift > (queue->size / 2)) &&
            ((new_size - queue->start) < queue->size)) {
          /* reset pointer shift */
          queue->data += queue->shift;

          /* move remaining data to the front */
          queue->data = memmove(queue->data,
                                queue->data + queue->start - queue->shift,
                                sizeof(void *) * (queue->end - queue->start + 1));

          /* reset start and pointer shifts */
          queue->shift  = queue->start;
          queue->data   -= queue->start;
        } else {
          /* increase stream buffer size */
          new_size += QUEUE_OVERHEAD;
          /* reset pointer shift */
          queue->data += queue->shift;

          /* reallocate memory blocks */
          queue->data = (void **)vrna_realloc(queue->data, sizeof(void *) * new_size);
          queue->size = new_size;

          /* restore pointer shift */
          queue->data -= queue->shift;
        }
      }

      for (i = queue->end + 1; i <= num; i++)
        queue->data[i] = NULL;

      queue->end = num;
#if VRNA_WITH_PTHREADS
      pthread_mutex_unlock(&queue->mtx);
#elif defined (_OPENMP)
    }

#endif
    }
  }
}


void
vrna_ostream_provide(struct vrna_ordered_stream_s *queue,
                     unsigned int                 i,
                     void                         *data)
{
  unsigned int j;

  if (queue) {
#if VRNA_WITH_PTHREADS
    pthread_mutex_lock(&queue->mtx);
#elif defined (_OPENMP)
#   pragma omp critical (vrna_ostream_access)
    {
#endif
    if ((queue->end < i) || (i < queue->start)) {
      vrna_message_warning(
        "vrna_ostream_provide(): data position (%d) out of range [%d:%d]!",
        i,
        queue->start,
        queue->end);
      return;
    }

    if (data) {
      /* store data */
      queue->data[i] = data;

      if (i == queue->start) {
#if VRNA_WITH_PTHREADS
        pthread_cond_signal(&queue->flush_cond);
#else
        /* process all consecutive blocks available from the start */
        flush_output(queue);
#endif
      }
    }

#if VRNA_WITH_PTHREADS
    pthread_mutex_unlock(&queue->mtx);
#elif defined (_OPENMP)
  }

#endif
  }
}
