#ifndef DIST_VARS_H
#define DIST_VARS_H
/*  Global variables for Distance-Package */

extern int   edit_backtrack;  /* set to 1 if you want backtracking */ 
   
extern char *aligned_line[4]; /* containes alignment after backtracking */

extern int  cost_matrix;     /* 0 usual costs (default), 1 Shapiro's costs */

/*  Global type defs for Distance-Package */

typedef struct {
                 int  type; 
                 int  weight;
                 int  father;
                 int  sons;
                 int  leftmostleaf;                 
               } Postorder_list;

typedef struct {
                 Postorder_list *postorder_list; 
                 int            *keyroots; 
               } Tree;

typedef struct {
                 int    type;
                 int    sign;
                 float  weight; 
               } swString;
#endif
