typedef struct {
        int   *split_list[2];
        int    split_size;
        double   isolation_index; } Split;

extern Split   *split_decomposition(float **dist);
extern void     free_Split(Split *x);
extern void     sort_Split(Split *x);
extern void     print_Split(Split *x);


/*     Auxiliary information contained in datatype 'Split': 
       Split[0].splitsize   .... contains the total number of splits.
*/
