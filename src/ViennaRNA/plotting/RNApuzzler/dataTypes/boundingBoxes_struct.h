#ifndef BOUNDING_BOXES_STRUCT_H
#define BOUNDING_BOXES_STRUCT_H

typedef struct boundingboxLoop {
    struct configtree* parent;

    // circle definition
    double c[2]; // center
    double r;    // radius
} loopBox;

typedef struct boundingboxStem {
    struct configtree* parent;

    // rectangle definition
    double a[2]; // direction 1 (unit vector) // direction from stem center to loop center
    double b[2]; // direction 2 (unit vector) // points to the left of vector a
    double c[2]; // center
    double e[2]; // half width extension of a and b

    // additional information on RNA
    int bulgeCount;
    double  bulgeDist;
    double** bulges;
} stemBox;

#endif
