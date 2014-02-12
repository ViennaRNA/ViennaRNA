/**
 *  \file edit_cost.h
 *  \brief global variables for Edit Costs included by treedist.c and stringdist.c
 */

#define PRIVATE static

PRIVATE char   sep    = ':';
PRIVATE char  *coding = "Null:U:P:H:B:I:M:S:E:R";

#define DIST_INF 10000  /* infinity */

typedef int CostMatrix[10][10];

PRIVATE CostMatrix *EditCost;  /* will point to UsualCost or ShapiroCost */

PRIVATE CostMatrix  UsualCost =
{

/*    Null,       U,        P,        H,        B,        I,        M,        S,        E,        R     */

   {        0,        1,        2,        2,        2,        2,        2,        1,        1, DIST_INF},   /* Null replaced */
   {        1,        0,        1, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF},   /* U    replaced */
   {        2,        1,        0, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF},   /* P    replaced */
   {        2, DIST_INF, DIST_INF,        0,        2,        2,        2, DIST_INF, DIST_INF, DIST_INF},   /* H    replaced */
   {        2, DIST_INF, DIST_INF,        2,        0,        1,        2, DIST_INF, DIST_INF, DIST_INF},   /* B    replaced */
   {        2, DIST_INF, DIST_INF,        2,        1,        0,        2, DIST_INF, DIST_INF, DIST_INF},   /* I    replaced */
   {        2, DIST_INF, DIST_INF,        2,        2,        2,        0, DIST_INF, DIST_INF, DIST_INF},   /* M    replaced */
   {        1, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF,        0, DIST_INF, DIST_INF},   /* S    replaced */
   {        1, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF,        0, DIST_INF},   /* E    replaced */
   { DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF,        0},   /* R    replaced */

};


PRIVATE CostMatrix ShapiroCost =
{

/*    Null,       U,        P,        H,        B,        I,        M,        S,        E,        R     */

   {        0,        1,        2,      100,        5,        5,       75,        5,        5, DIST_INF},   /* Null replaced */
   {        1,        0,        1, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF},   /* U    replaced */
   {        2,        1,        0, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF},   /* P    replaced */
   {      100, DIST_INF, DIST_INF,        0,        8,        8,        8, DIST_INF, DIST_INF, DIST_INF},   /* H    replaced */
   {        5, DIST_INF, DIST_INF,        8,        0,        3,        8, DIST_INF, DIST_INF, DIST_INF},   /* B    replaced */
   {        5, DIST_INF, DIST_INF,        8,        3,        0,        8, DIST_INF, DIST_INF, DIST_INF},   /* I    replaced */
   {       75, DIST_INF, DIST_INF,        8,        8,        8,        0, DIST_INF, DIST_INF, DIST_INF},   /* M    replaced */
   {        5, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF,        0, DIST_INF, DIST_INF},   /* S    replaced */
   {        5, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF,        0, DIST_INF},   /* E    replaced */
   { DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF, DIST_INF,        0},   /* R    replaced */

};

