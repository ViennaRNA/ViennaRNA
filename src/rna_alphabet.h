#ifndef _RNA_ALPHABET_
#define _RNA_ALPHABET_

#include <string>



typedef char RNA_Alphabet;

/** pair alphabet of RNA_Alphabet */
typedef struct {
    RNA_Alphabet a;
    RNA_Alphabet b;
} RNA_AlphaPair;

const int ALPHA_BASE_A='a';
const int ALPHA_BASE_C='c';
const int ALPHA_BASE_G='g';
const int ALPHA_BASE_U='u';
const int ALPHA_BASEPAIR='P';
const int ALPHA_GAP='-';
const int ALPHA_BASE='B';


const int RNA_ALPHABET_SIZE=8;

const int ALPHA_PRO_BASE_A=0;
const int ALPHA_PRO_BASE_C=1;
const int ALPHA_PRO_BASE_G=2;
const int ALPHA_PRO_BASE_U=3;
const int ALPHA_PRO_BASEPAIR=4;
const int ALPHA_PRO_GAP=5;
const int ALPHA_PRO_BASE=6;
const int ALPHA_PRO_ANCHOR=7;

struct RNA_Alphabet_Profile {
    RNA_Alphabet_Profile() {};

    double p[RNA_ALPHABET_SIZE];
    std::string columnStr;
};

int alpha2RNA_Alpha(char c);

#endif
