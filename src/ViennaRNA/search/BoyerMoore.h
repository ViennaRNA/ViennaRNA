#ifndef VIENNA_RNA_PACKAGE_SEARCH_BOYER_MOORE_H
#define VIENNA_RNA_PACKAGE_SEARCH_BOYER_MOORE_H

/**
 *  @file     ViennaRNA/search/BoyerMoore.h
 *  @ingroup  utils, search_utils
 *  @brief    Variants of the Boyer-Moore string search algorithm
 */

/**
 *  @addtogroup   search_utils
 *  @{
 */

/**
 *  @brief  Search for a string of elements in a larger string of elements using
 *          the Boyer-Moore-Horspool algorithm
 *
 *  To speed-up subsequent searches with this function, the Bad Character Table
 *  should be precomputed and passed as argument @p badchars.
 *
 *  @see    vrna_search_BM_BCT_num(), vrna_search_BMH()
 *
 *  @param  needle        The pattern of object representations to search for
 *  @param  needle_size   The size (length) of the pattern provided in @p needle
 *  @param  haystack      The string of objects the search will be performed on
 *  @param  haystack_size The size (length) of the @p haystack string
 *  @param  start         The position within @p haystack where to start the search
 *  @param  badchars      A pre-computed Bad Character Table obtained from vrna_search_BM_BCT_num()
 *                        (If NULL, a Bad Character Table will be generated automatically)
 *  @param  cyclic        Allow for cyclic matches if non-zero, stop search at end of haystack otherwise
 *  @return               A pointer to the first occurence of @p needle within @p haystack after position @p start
 */
const unsigned int *
vrna_search_BMH_num(const unsigned int  *needle,
                    size_t              needle_size,
                    const unsigned int  *haystack,
                    size_t              haystack_size,
                    size_t              start,
                    size_t              *badchars,
                    unsigned char       cyclic);


/**
 *  @brief  Search for an ASCII pattern within a larger ASCII string using the
 *          Boyer-Moore-Horspool algorithm
 *
 *  To speed-up subsequent searches with this function, the Bad Character Table
 *  should be precomputed and passed as argument @p badchars. Furthermore, both,
 *  the lengths of @p needle and the length of @p haystack should be pre-computed
 *  and must be passed along with each call.
 *
 *  @see    vrna_search_BM_BCT(), vrna_search_BMH_num()
 *
 *  @param  needle        The NULL-terminated ASCII pattern to search for
 *  @param  needle_size   The size (length) of the pattern provided in @p needle
 *  @param  haystack      The NULL-terminated ASCII string of the search will be performed on
 *  @param  haystack_size The size (length) of the @p haystack string
 *  @param  start         The position within @p haystack where to start the search
 *  @param  badchars      A pre-computed Bad Character Table obtained from vrna_search_BM_BCT()
 *                        (If NULL, a Bad Character Table will be generated automatically)
 *  @param  cyclic        Allow for cyclic matches if non-zero, stop search at end of haystack otherwise
 *  @return               A pointer to the first occurence of @p needle within @p haystack after position @p start
 */
const char *
vrna_search_BMH(const char    *needle,
                size_t        needle_size,
                const char    *haystack,
                size_t        haystack_size,
                size_t        start,
                size_t        *badchars,
                unsigned char cyclic);


/**
 *  @brief  Retrieve a Boyer-Moore Bad Character Table for a pattern of elements
 *          represented by natural numbers
 *
 *  @note   We store the maximum number representation of an element @p num_max
 *          at position @p 0. So the actual bad character table @p T starts at
 *          @p T[1] for an element represented by number @p 0.
 *
 *  @see  vrna_search_BMH_num(), vrna_search_BM_BCT()
 *
 *  @param  pattern       The pattern of element representations used in the subsequent search
 *  @param  pattern_size  The size (length) of the pattern provided in @p pattern
 *  @param  num_max       The maximum number representation of an element, i.e. the size of the alphabet
 *  @return               A Bad Character Table for use in our Boyer-Moore search algorithm implementation(s)
 */
size_t *
vrna_search_BM_BCT_num(const unsigned int *pattern,
                       size_t             pattern_size,
                       unsigned int       num_max);


/**
 *  @brief  Retrieve a Boyer-Moore Bad Character Table for a NULL-terminated pattern of ASCII characters
 *
 *  @note   We store the maximum number representation of an element, i.e. @p 127
 *          at position @p 0. So the actual bad character table @p T starts at
 *          @p T[1] for an element represented by ASCII code @p 0.
 *
 *  @see  vrna_search_BMH(), vrna_search_BM_BCT_num()
 *
 *  @param  pattern       The NULL-terminated pattern of ASCII characters used in the subsequent search
 *  @return               A Bad Character Table for use in our Boyer-Moore search algorithm implementation(s)
 */
size_t *
vrna_search_BM_BCT(const char *pattern);


/**
 *  @}
 */
#endif
