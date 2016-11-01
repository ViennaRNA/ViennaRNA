#ifndef _MISC_H
#define _MISC_H

#define DELETE(T)   if(T) \
                    { \
                      delete T; \
                      T=NULL; \
                    }

#define DELETE_ARRAY(T)  if(T) \
                         { \
                           delete[] T; \
                           T=NULL; \
                         }

template <class T>
void showArray(T *array, unsigned int m,unsigned int n);

#endif
