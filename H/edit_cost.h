/*   

cost.h   ::  global variables for Edit Costs
             included by treedist.c and stringdist.c

*/
#define PRIVATE static

PRIVATE char   sep    = ':';
PRIVATE char  *coding = "Null:U:P:H:B:I:M:S:E:R";

#define  INF 10000  /* infinity */  

typedef int CostMatrix[10][10];

PRIVATE CostMatrix *EditCost;  /* will point to UsualCost or ShapiroCost */

PRIVATE CostMatrix  UsualCost = 
{

/*  Null,   U,   P,   H,   B,   I,   M,   S,   E,   R     */

   {   0,   1,   2,   2,   2,   2,   2,   1,  1,  INF},   /* Null replaced */
   {   1,   0,   1, INF, INF, INF, INF, INF, INF, INF},   /* U    replaced */
   {   2,   1,   0, INF, INF, INF, INF, INF, INF, INF},   /* P    replaced */
   {   2, INF, INF,   0,   2,   2,   2, INF, INF, INF},   /* H    replaced */
   {   2, INF, INF,   2,   0,   1,   2, INF, INF, INF},   /* B    replaced */
   {   2, INF, INF,   2,   1,   0,   2, INF, INF, INF},   /* I    replaced */
   {   2, INF, INF,   2,   2,   2,   0, INF, INF, INF},   /* M    replaced */
   {   1, INF, INF, INF, INF, INF, INF,   0, INF, INF},   /* S    replaced */
   {   1, INF, INF, INF, INF, INF, INF, INF,   0, INF},   /* E    replaced */
   { INF, INF, INF, INF, INF, INF, INF, INF, INF,   0},   /* R    replaced */

};


PRIVATE CostMatrix ShapiroCost = 
{

/*  Null,   U,   P,   H,   B,   I,   M,   S,   E,  R     */

   {   0,   1,   2, 100,   5,   5,  75,   5,   5, INF},   /* Null replaced */
   {   1,   0,   1, INF, INF, INF, INF, INF, INF, INF},   /* U    replaced */
   {   2,   1,   0, INF, INF, INF, INF, INF, INF, INF},   /* P    replaced */
   { 100, INF, INF,   0,   8,   8,   8, INF, INF, INF},   /* H    replaced */
   {   5, INF, INF,   8,   0,   3,   8, INF, INF, INF},   /* B    replaced */
   {   5, INF, INF,   8,   3,   0,   8, INF, INF, INF},   /* I    replaced */
   {  75, INF, INF,   8,   8,   8,   0, INF, INF, INF},   /* M    replaced */
   {   5, INF, INF, INF, INF, INF, INF,   0, INF, INF},   /* S    replaced */
   {   5, INF, INF, INF, INF, INF, INF, INF,   0, INF},   /* E    replaced */
   { INF, INF, INF, INF, INF, INF, INF, INF, INF,   0},   /* R    replaced */

};

