typedef struct{
        int   set1;
        int   set2;
        float distance;
        float distance2;
        } Union;

extern Union *wards_cluster(float **clmat);
extern Union *neighbour_joining(float **clmat);
extern void   printf_phylogeny(Union *tree, char *type);


/* Auxiliary information in a Union tree: 
      tree[0].set1 contains the number of elements of tree, i.e.,
      tree[tree[0].set1-1]] is the last one !!
*/
