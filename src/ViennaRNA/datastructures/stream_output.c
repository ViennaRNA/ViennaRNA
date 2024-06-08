#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if VRNA_WITH_PTHREADS
# include <pthread.h>
#endif

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/datastructures/stream_output.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define QUEUE_OVERHEAD  32


struct vrna_ordered_stream_s {
  unsigned int                start;      /* first element index in queue, i.e. start of queue */
  unsigned int                end;        /* last element index in queue */
  unsigned int                size;       /* available memory size for 'data' and 'provided' */
  unsigned int                shift;      /* pointer offset for 'data' and 'provided' */

  vrna_stream_output_f output;    /* callback to execute if consecutive elements from head are available */
  void                        **data;     /* actual data passed to the callback */
  unsigned char               *provided;  /* for simplicity we use unsigned char instead of single bits per element */
  void                        *auxdata;   /* auxiliary data passed to the callback */
#if VRNA_WITH_PTHREADS
  pthread_mutex_t             mtx;        /* semaphore to provide concurrent access */
#endif
};


PRIVATE INLINE void
flush_output(struct vrna_ordered_stream_s *queue)
{
  unsigned int j;

  /* flush all consecutive blocks available from the start of queue */

  /* 1. process output callback */
  if (queue->output)
    for (j = queue->start; (j <= queue->end) && (queue->provided[j]); j++)
      queue->output(queue->auxdata, j, queue->data[j]);

  /* 2. move start of queue */
  while ((queue->start <= queue->end) && (queue->provided[queue->start]))
    queue->start++;

  /* 3. set empty queue condition if necessary */
  if (queue->end < queue->start) {
    queue->provided[queue->start] = 0;
    queue->end                    = queue->start;
  }
}


PUBLIC struct vrna_ordered_stream_s *
vrna_ostream_init(vrna_stream_output_f output,
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
  queue->data     = (void **)vrna_alloc(sizeof(void *) * QUEUE_OVERHEAD);
  queue->provided = (unsigned char *)vrna_alloc(sizeof(unsigned char) * QUEUE_OVERHEAD);

#if VRNA_WITH_PTHREADS
  pthread_mutex_init(&queue->mtx, NULL);
#endif

  return queue;
}


PUBLIC void
vrna_ostream_free(struct vrna_ordered_stream_s *queue)
{
  if (queue) {
#if VRNA_WITH_PTHREADS
    pthread_mutex_lock(&queue->mtx);
#endif

    flush_output(queue);

#if VRNA_WITH_PTHREADS
    pthread_mutex_unlock(&queue->mtx);
#endif

    /* free remaining memory */
    queue->data     += queue->shift;
    queue->provided += queue->shift;
    free(queue->data);
    free(queue->provided);

    /* free ostream itself */
    free(queue);
  }
}


PUBLIC int
vrna_ostream_threadsafe(void)
{
#if VRNA_WITH_PTHREADS
  return 1;
#else
  return 0;
#endif
}


PUBLIC void
vrna_ostream_request(struct vrna_ordered_stream_s *queue,
                     unsigned int                 num)
{
  unsigned int i;

  if (queue) {
#if VRNA_WITH_PTHREADS
    pthread_mutex_lock(&queue->mtx);
#endif
    if (num >= queue->end) {
      /* check whether we have to increase memory */
      unsigned int new_size = num - queue->shift + 1;

      if (queue->size < new_size + 1) {
        /*
         *  Check whether we can simply move data around to get more space.
         *  We do this only, if more than half of the first buffer is empty.
         *  Otherwise, we simply increase the buffer to fit incoming data.
         */
        unsigned int mem_unavail = queue->start - queue->shift;
        if ((mem_unavail > (queue->size / 2)) &&
            (new_size - mem_unavail < queue->size + 1)) {
          /* reset pointer shift */
          queue->data     += queue->shift;
          queue->provided += queue->shift;

          /* move remaining data to the front */
          queue->data = memmove(queue->data,
                                queue->data + mem_unavail,
                                sizeof(void *) * (queue->end - queue->start + 1));

          /* move provider bitmask to the front */
          queue->provided = memmove(queue->provided,
                                    queue->provided + mem_unavail,
                                    sizeof(unsigned char) * (queue->end - queue->start + 1));

          /* reset start and pointer shifts */
          queue->shift    = queue->start;
          queue->data     -= queue->start;
          queue->provided -= queue->start;
        } else {
          /* increase stream buffer size */
          new_size += QUEUE_OVERHEAD;
          /* reset pointer shift */
          queue->data     += queue->shift;
          queue->provided += queue->shift;

          /* reallocate memory blocks */
          queue->data     = (void **)vrna_realloc(queue->data, sizeof(void *) * new_size);
          queue->provided =
            (unsigned char *)vrna_realloc(queue->provided, sizeof(void *) * new_size);
          queue->size = new_size;

          /* restore pointer shift */
          queue->data     -= queue->shift;
          queue->provided -= queue->shift;
        }
      }

      for (i = queue->end + 1; i <= num; i++)
        queue->provided[i] = 0;

      queue->end = num;
    }

#if VRNA_WITH_PTHREADS
    pthread_mutex_unlock(&queue->mtx);
#endif
  }
}


PUBLIC void
vrna_ostream_provide(struct vrna_ordered_stream_s *queue,
                     unsigned int                 i,
                     void                         *data)
{
  if (queue) {
#if VRNA_WITH_PTHREADS
    pthread_mutex_lock(&queue->mtx);
#endif
    if ((queue->end < i) || (i < queue->start)) {
      vrna_log_warning(
        "vrna_ostream_provide(): data position (%d) out of range [%d:%d]!",
        i,
        queue->start,
        queue->end);
      return;
    }

    /* store data */
    queue->data[i]      = data;
    queue->provided[i]  = 1;

    /* process all consecutive blocks available from the start */
    if (i == queue->start)
      flush_output(queue);

#if VRNA_WITH_PTHREADS
    pthread_mutex_unlock(&queue->mtx);
#endif
  }
}
