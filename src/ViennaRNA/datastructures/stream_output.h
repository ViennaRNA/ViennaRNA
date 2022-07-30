#ifndef VIENNA_RNA_PACKAGE_STREAM_OUTPUT_H
#define VIENNA_RNA_PACKAGE_STREAM_OUTPUT_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file     ViennaRNA/datastructures/stream_output.h
 *  @ingroup  utils, buffer_utils
 *  @brief    An implementation of a buffered, ordered stream output data structure
 */

/**
 *  @addtogroup   buffer_utils
 *  @{
 */

/**
 *  @brief  An ordered output stream structure with unordered insert capabilities
 */
typedef struct vrna_ordered_stream_s *vrna_ostream_t;

/**
 *  @brief  Ordered stream processing callback
 *
 *  This callback will be processed in sequential order as soon as sequential
 *  data in the output stream becomes available.
 *
 *  @note The callback must also release the memory occupied by the
 *        data passed since the stream will lose any reference to it
 *        after the callback has been executed.
 *
 *  @param  auxdata   A shared pointer for all calls, as provided by the second argument to vrna_ostream_init()
 *  @param  i         The index number of the data passed to @p data
 *  @param  data      A block of data ready for processing
 */
typedef void (*vrna_stream_output_f)(void        *auxdata,
                                            unsigned int i,
                                            void         *data);
DEPRECATED(typedef void (vrna_callback_stream_output)(void         *auxdata,
                                           unsigned int i,
                                           void         *data),
           "Use vrna_stream_output_f instead!");


/**
 *  @brief  Get an initialized ordered output stream
 *
 *  @see  vrna_ostream_free(), vrna_ostream_request(), vrna_ostream_provide()
 *
 *  @param  output    A callback function that processes and releases data in the stream
 *  @param  auxdata   A pointer to auxiliary data passed as first argument to the @p output callback
 *  @return           An initialized ordered output stream
 */
vrna_ostream_t
vrna_ostream_init(vrna_stream_output_f output,
                  void                        *auxdata);


/**
 *  @brief  Free an initialized ordered output stream
 *
 *  @see vrna_ostream_init()
 *
 *  @param  dat   The output stream for which occupied memory should be free'd
 */
void
vrna_ostream_free(vrna_ostream_t dat);


int
vrna_ostream_threadsafe(void);


/**
 *  @brief  Request index in ordered output stream
 *
 *  This function must be called prior to vrna_ostream_provide() to
 *  indicate that data associted with a certain index number is expected
 *  to be inserted into the stream in the future.
 *
 *  @see vrna_ostream_init(), vrna_ostream_provide(), vrna_ostream_free()
 *
 *  @param  dat   The output stream for which the index is requested
 *  @param  num   The index to request data for
 */
void
vrna_ostream_request(vrna_ostream_t dat,
                     unsigned int   num);


/**
 *  @brief  Provide output stream data for a particular index
 *
 *  @pre  The index data is provided for must have been requested using
 *        vrna_ostream_request() beforehand.
 *
 *  @see  vrna_ostream_request()
 *
 *  @param  dat   The output stream for which data is provided
 *  @param  i     The index of the provided data
 *  @param  data  The data provided
 */
void
vrna_ostream_provide(vrna_ostream_t dat,
                     unsigned int   i,
                     void           *data);


/**
 *  @}
 */


#endif
