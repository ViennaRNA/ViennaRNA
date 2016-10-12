#include "rna_alphabet.h"

#include <cassert>

int alpha2RNA_Alpha(char c) {
    switch (c) {
    case ALPHA_BASE_A:
        return ALPHA_PRO_BASE_A;
    case ALPHA_BASE_C:
        return ALPHA_PRO_BASE_C;
    case ALPHA_BASE_G:
        return ALPHA_PRO_BASE_G;
    case ALPHA_BASE_U:
        return ALPHA_PRO_BASE_U;
    case ALPHA_GAP:
        return ALPHA_PRO_GAP;
    case ALPHA_BASEPAIR:
        return ALPHA_PRO_BASEPAIR;
    case ALPHA_BASE:
        return ALPHA_PRO_BASE;
    default:
        assert(true);
    }

    return 0;	// never reached !!
}





