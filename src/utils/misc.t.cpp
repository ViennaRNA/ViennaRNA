#ifndef _MISC_T_CPP_
#define _MISC_T_CPP_

// #include <stdio.h>

#include "misc.h"

template <class T>
void showArray(T *array, unsigned int m,unsigned int n) {
//  char buf[20];
    unsigned int i,maxElemLen;
    T maxValue;

    // check for the largest value
    maxValue=0;
    for (i=0; i<m*n; i++)
        maxValue = max(maxValue,array[i]);

    //snprintf(buf,20,"%.02f",maxValue);
    //maxElemLen=strlen(buf)+1;
    maxElemLen=5;

    for (i=0; i<m*n; i++) {
        if (i%n==0)
            std::cout << std::endl;

        //sprintf(buf,"%*.02f ",maxElemLen,array[i]);
        std::cout << "[" << i << "]" << array[i] << " ";
//    std::cout << buf;
    }
    std::cout << std::endl;
}


#endif
