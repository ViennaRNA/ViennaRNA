/*
		Tree edit distances for RNA secondary structures
		Walter Fontana, Ivo L Hofacker, Peter F Stadler
			     Vienna RNA Package
*/
/* Last changed Time-stamp: <97/10/27 15:23:48 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include  "edit_cost.h"
#include  "dist_vars.h"
#include  "utils.h"
static char rcsid[] = "$Id: treedist.c,v 1.3 1997/11/03 10:39:43 ivo Rel $";

#define PRIVATE  static
#define PUBLIC

#define MNODES	  4000	  /* Maximal number of nodes for alignment    */

PUBLIC  Tree     *make_tree(char *struc);
PUBLIC  float     tree_edit_distance(Tree *T1, Tree *T2);
PUBLIC  void      print_tree(Tree *t);
PUBLIC  void      free_tree(Tree *t);


PRIVATE void            tree_dist(int i, int j);
PRIVATE int             edit_cost(int i, int j);
PRIVATE int            *make_keyroots(Postorder_list *pl);
PRIVATE void            sort(int n, int *ra);
PRIVATE Postorder_list *make_postorder_list(char *struc);
PRIVATE int             decode(char *id);
PRIVATE int             number_of_nodes(char *struc);
PRIVATE void            encode(int type, char *label);
PRIVATE void            print_keyroots(int *keyroots);
PRIVATE void            print_postorder_list(Postorder_list *pl);
PRIVATE void            backtracking(void);
PRIVATE void            sprint_aligned_trees(void);


PRIVATE  Tree  *tree1, *tree2;
PRIVATE  int  **tdist;       /* contains distances between subtrees */
PRIVATE  int  **fdist;       /* contains distances between forests */
PRIVATE  int  *alignment[2]; /* contains numeric information on the alignment:
				alignment[0][p], aligment[1][p] are aligned postions.
				INDELs have one 0.
				alignment[0][0] contains the length of the alignment. */

/*---------------------------------------------------------------------------*/

PUBLIC float tree_edit_distance(Tree *T1, Tree *T2)
{
   int i1,j1,i,j, dist;
   int n1, n2;

   if (cost_matrix==0) EditCost = &UsualCost;
   else EditCost = &ShapiroCost;
   
   n1 = T1->postorder_list[0].sons;
   n2 = T2->postorder_list[0].sons;
   
   tdist = (int **) space(sizeof(int *) * (n1+1));
   fdist = (int **) space(sizeof(int *) * (n1+1));
   for (i=0; i<=n1; i++) { 
      tdist[i] = (int *) space(sizeof(int) * (n2+1));
      fdist[i] = (int *) space(sizeof(int) * (n2+1));
   }
   
   tree1 = T1;  tree2 = T2;
   
   for (i1 = 1; i1 <= T1->keyroots[0]; i1++) {
      i = T1->keyroots[i1];
      for (j1 = 1; j1 <= T2->keyroots[0]; j1++) {
	 j = T2->keyroots[j1];

	 tree_dist(i,j);
      }
   }

   if (edit_backtrack) {

      if ((n1>MNODES)||(n2>MNODES)) nrerror("tree too large for alignment");
      
      alignment[0] = (int *) space((n1+1)*sizeof(int));
      alignment[1] = (int *) space((n2+1)*sizeof(int));
      
      backtracking();
      sprint_aligned_trees();
      free(alignment[0]);
      free(alignment[1]);
   }
   dist = tdist[n1][n2];
   for (i=0; i<=n1; i++) { 
      free(tdist[i]);
      free(fdist[i]);
   }
   free(tdist);
   free(fdist);
   
   return (float) dist;
}

/*---------------------------------------------------------------------------*/

PRIVATE void tree_dist(int i, int j)
{
   int li,lj,i1,j1,i1_1,j1_1,li1_1,lj1_1,f1,f2,f3,f;
   int cost, lleaf_i1, lleaf_j1;
   
   fdist[0][0] = 0;
   
   li = tree1->postorder_list[i].leftmostleaf;
   lj = tree2->postorder_list[j].leftmostleaf;
   
   for (i1 = li; i1 <= i; i1++) {
      i1_1 = (li == i1 ? 0 : i1-1);
      fdist[i1][0] = fdist[i1_1][0] + edit_cost(i1, 0);
   }
   
   for (j1 = lj; j1 <= j; j1++) {
      j1_1 = (lj == j1 ? 0 : j1-1);
      fdist[0][j1] = fdist[0][j1_1] + edit_cost(0, j1);
   }
   
   for (i1 = li; i1 <= i; i1++) {
      
      lleaf_i1 = tree1->postorder_list[i1].leftmostleaf;
      li1_1 = (li > lleaf_i1-1 ? 0 : lleaf_i1-1);
      i1_1 = (i1 == li ? 0: i1-1);
      cost = edit_cost(i1, 0);
      
      for (j1 = lj; j1 <= j; j1++) {
	 
	 lleaf_j1 = tree2->postorder_list[j1].leftmostleaf;
	 j1_1 = (j1 == lj ? 0: j1-1);
	 
	 f1 = fdist[i1_1][j1] + cost;
	 f2 = fdist[i1][j1_1] + edit_cost(0, j1);
	 
	 f = f1 < f2 ? f1 : f2;
	 
	 if (lleaf_i1 == li && lleaf_j1 == lj)   {
	    
	    f3 = fdist[i1_1][j1_1] + edit_cost(i1, j1);
	    
	    fdist[i1][j1] = f3 < f ? f3 : f;
	    
	    tdist[i1][j1] = fdist[i1][j1];   /* store in array permanently */
	 }
	 else {
	    lj1_1 = (lj > lleaf_j1-1 ? 0 : lleaf_j1-1);
	    
	    f3 = fdist[li1_1][lj1_1] + tdist[i1][j1];
	    
	    fdist[i1][j1] = f3 < f ? f3 : f;
	 }
      }
   }
}

/*---------------------------------------------------------------------------*/

PRIVATE int edit_cost(int i, int j)
{
   int  c, diff, cd, min, a, b;
   
   c = (*EditCost)[tree1->postorder_list[i].type][tree2->postorder_list[j].type];
   
   diff = abs((a=tree1->postorder_list[i].weight) - (b=tree2->postorder_list[j].weight));
   
   min = (a < b ? a: b);
   if (min == a) cd = (*EditCost)[0][tree2->postorder_list[j].type];
   else          cd = (*EditCost)[0][tree1->postorder_list[i].type];
   
   return (c * min + cd * diff);
   
}

/*---------------------------------------------------------------------------*/

PUBLIC Tree *make_tree(char *struc)
{
   Tree *tree;
   
   tree = (Tree *) space(sizeof(Tree));
   
   tree->postorder_list = make_postorder_list(struc);
   tree->keyroots       = make_keyroots(tree->postorder_list);
   
   return (tree);
}


/*---------------------------------------------------------------------------*/

PRIVATE int  *make_keyroots(Postorder_list *pl)
{
   int   i, k, keys;
   int  *keyroots;
   
   keyroots = (int *) space(sizeof(int)*(pl[0].sons+1));
   keys      = 0;
   
   for (i = 1; i <= pl[0].sons; i++) {
      if (!pl[i].sons) {
	 
	 /* leaf */
	 
	 k = pl[0].sons;
	 while (pl[k].leftmostleaf != i) k--;
	 keyroots[++keys] = k;
      }
   }
   
   sort(keys, keyroots);
   keyroots[0] = keys;
   
   return (keyroots);
}

/*---------------------------------------------------------------------------*/

PRIVATE void sort(int n, int *ra)       /* heap sort,  indices are 1..n !!! */
{
   int l,j,ir,i;
   int rra;
   
   if (n == 1) return;
   
   l = (n >> 1)+1;
   ir = n;
   for (;;) {
      if (l > 1)
	 rra = ra[--l];
      else {
	 rra = ra[ir];
	 ra[ir] = ra[1];
	 if (--ir == 1) {
	    ra[1] = rra;
	    return;
	 }
      }
      i = l;
      j = l << 1;
      while (j <= ir) {
	 if (j < ir && ra[j] < ra[j+1]) ++j;
	 if (rra < ra[j]) {
	    ra[i] = ra[j];
	    j += (i = j);
	 }
	 else j = ir+1;
      }
      ra[i] = rra;
   }
}

/*---------------------------------------------------------------------------*/

PRIVATE Postorder_list *make_postorder_list(char *struc)

     /*
       Convention for structure representation "struc":
       Nodes are one pair of matching parentheses, with the type and possibly
       a weight of the node immediately preceding the closing parentheses.
       
       Types: 
       
       U....unpaired
       P....paired
       H....hairpin loop
       B....bulge loop
       I....internal loop
       M....multiloop
       S....stack
       R....virtual root
       
       Example:
       
       .((..(((...)))..((..)))). in usual notation becomes:
       
       full tree:
       ((U)(((U)(U)((((U)(U)(U)P)P)P)(U)(U)(((U)(U)P)P)P)P)(U)R)
       HIT:
       ((U1)((U2)((U3)P3)(U2)((U2)P2)P2)(U1)R)
       Shapiro:
       (((((H)S)((H)S)M)S)R)
       
       */
{
   int  paren, i, l, order, local_order, w, sons, count;
   int  n_nodes, p;
   char id[100];
   Postorder_list *pl;
   int  match_pos[MNODES], match_order[MNODES];
   
   
   n_nodes = number_of_nodes(struc);
   if (n_nodes>MNODES) nrerror("structure too long in make_postorder_list");
   pl = (Postorder_list *) space(sizeof(Postorder_list)*(n_nodes+1));
   pl[0].sons = n_nodes;
   
   paren = 1;
   match_pos[paren]   = 0;
   match_order[paren] = 0;
   
   i     = 1;
   l     = 0;
   order = 0;
   
   while (paren) {
      switch (struc[i]) {
       case '(': 
	 match_pos[++paren] = i;
	 match_order[paren] = order;
	 break;
       case ')':
	 order++;
	 id[l] = '\0';
	 l = 0;
	 while (isalpha((int) id[l])) l++;
	 if (id[l]) sscanf(id+l, "%d", &w);
	 else  w = 1;
	 id[l] = '\0';
	 pl[order].type         = decode(id);
	 pl[order].weight       = w;
	 pl[order].leftmostleaf = match_order[paren]+1;
	 
	 sons = count = 0;
	 local_order  = match_order[paren];
	 for (p = match_pos[paren]+1; p < i; p++) {
	    if (struc[p] == '(') count++;
	    else if (struc[p] == ')') {
	       local_order++;
	       if (count == 1) {
		  sons++;
		  pl[local_order].father = order;
	       }
	       count--;
	    }
	 }
	 
	 pl[order].sons         = sons;
	 paren--;
	 l = 0;
	 break;
       default:
	 id[l++] = struc[i];
	 break;
      }
      i++;
   }
   
   return (pl);
}

/*---------------------------------------------------------------------------*/

PRIVATE int decode(char *id)
{
   int   n, quit, i;
   char  label[100], *code;
   
   n = 0;
   
   quit = 0;
   code = coding;
   
   while (!quit) {
      for (i = 0; code[i] != sep; i++) {
	 if (code[i] == '\0') {
	    quit = 1;
	    break;
	 }
	 label[i] = code[i];
      }
      label[i] = '\0';
      if (strcmp(id, label) == 0) return (n);
      code += (i+1);
      n++;
   }
   
   fprintf(stderr,"Syntax error: node identifier \"%s\" not found "
		  "in coding string \"%s\"\n", id, coding);
   fprintf(stderr, "Exiting...");
   exit(0);
}

/*---------------------------------------------------------------------------*/

PRIVATE void encode(int type, char label[])
{
   int   i, l;
   
   l = 0;
   for (i = 0; i < type; i++) {
      while (coding[l] != sep && coding[l]) l++;
      l++;
   }
   
   for (i = 0; coding[l+i] != sep; i++) {
      if (coding[l+i] == '\0') break;
      label[i] = coding[l+i];
   }
   label[i] = '\0';
}

/*---------------------------------------------------------------------------*/

PRIVATE int number_of_nodes(char *struc)
{
   int  l, c, i;
   
   l = strlen(struc);
   for (c = 0, i = 0; i < l; i++) if (struc[i] == ')') c++;
   return (c);
}

/*---------------------------------------------------------------------------*/

PRIVATE void print_keyroots(int *keyroots)
{
   int i;
   
   printf("--->  key roots  <---\n\n");
   
   printf("entries: %d\n", keyroots[0]);
   printf("{");
   for (i = 1; i <= keyroots[0]; i++) printf(" %d", keyroots[i]);
   printf(" }\n\n");
}

/*---------------------------------------------------------------------------*/

PRIVATE void print_postorder_list(Postorder_list *pl)
{
   register i;
   char     label[100];
   
   printf("--->  postorder list  <---\n\n");
   
   for (i = 1; i <= pl[0].sons; i++) {
      printf("    postorder: %3d\n", i);
      *label = '\0';
      encode(pl[i].type, label);
      printf("         type: %3d (%s)\n", pl[i].type, label);
      printf("       weight: %3d\n", pl[i].weight);
      printf("       father: %3d\n", pl[i].father);
      printf("         sons: %3d\n", pl[i].sons);
      printf("leftmost leaf: %3d\n", pl[i].leftmostleaf);
      printf("\n");
   }
}

/*---------------------------------------------------------------------------*/

PUBLIC void print_tree(Tree *t)
{
   print_postorder_list(t->postorder_list);
   print_keyroots(t->keyroots);
   fflush(stdout);
}

/*---------------------------------------------------------------------------*/

PUBLIC void free_tree(Tree *t)
{
   free(t->postorder_list);
   free(t->keyroots);
   free(t);
}

/*---------------------------------------------------------------------------*/


PRIVATE void backtracking(void)
{
   int li,lj,i1,j1,i1_1,j1_1,li1_1,lj1_1,f;
   int cost, lleaf_i1, lleaf_j1, ss, i,j,k;
   struct {int i,j;} sector[MNODES/2];
   
   ss=0;
   
   i=i1=tree1->postorder_list[0].sons;
   j=j1=tree2->postorder_list[0].sons;
   
 start:
   li = tree1->postorder_list[i].leftmostleaf;
   lj = tree2->postorder_list[j].leftmostleaf;
   
   
   while ((i1>=li)&&(j1>=lj)) {
      
      lleaf_i1 = tree1->postorder_list[i1].leftmostleaf;
      li1_1 = (li > lleaf_i1-1 ? 0 : lleaf_i1-1);
      i1_1 = (i1 == li ? 0: i1-1);
      lleaf_j1 = tree2->postorder_list[j1].leftmostleaf;
      lj1_1 = (lj > lleaf_j1-1 ? 0 : lleaf_j1-1);       
      j1_1 = (j1 == lj ? 0: j1-1);
      
      f = fdist[i1][j1];
      
      cost = edit_cost(i1, 0);
      if (f == fdist[i1_1][j1] + cost) {
	 alignment[0][i1]=0;
	 i1=i1_1;
      }
      else {
	 if (f ==  fdist[i1][j1_1] + edit_cost(0, j1)) {
	    alignment[1][j1]=0;
	    j1=j1_1;
	 }
	 else if (lleaf_i1 == li && lleaf_j1 == lj) {
	    alignment[0][i1] = j1;
	    alignment[1][j1] = i1;
	    i1=i1_1; j1=j1_1;
	 } else  {
	    sector[ss].i=i1;
	    sector[ss++].j=j1;
	    i1=li1_1;
	    j1=lj1_1;
	 }
      }
   }  
   for (; i1>=li; ) {
      alignment[0][i1]=0;
      i1 = (i1 == li ? 0: i1-1);
   }
   for (; j1>=lj; ) {
      alignment[1][j1]=0;
      j1 = (j1 == lj ? 0: j1-1);
   }
   while (ss>0) {
      i1=sector[--ss].i;
      j1=sector[ss].j;
      for (k=1; 1; k++) {
	 i = tree1->keyroots[k];
	 if (tree1->postorder_list[i].leftmostleaf ==
	     tree1->postorder_list[i1].leftmostleaf) break;
      }
      for (k=1; 1; k++) {
	 j = tree2->keyroots[k];
	 if (tree2->postorder_list[j].leftmostleaf ==
	     tree2->postorder_list[j1].leftmostleaf) break;
      }
      tree_dist(i,j);
      goto start;
   }
}

/*---------------------------------------------------------------------------*/

PRIVATE void sprint_aligned_trees(void)
{
   int i,j,n1,n2,k,l,p, ni, nj, weights;
   char t1[2*MNODES+1], t2[2*MNODES+1], a1[8*MNODES], a2[8*MNODES], ll[20], ll1[20];

   weights=0;
   n1=tree1->postorder_list[0].sons;
   n2=tree2->postorder_list[0].sons;
   for (i=1; i<=n1; i++) weights |= (tree1->postorder_list[i].weight!=1);
   for (i=1; i<=n2; i++) weights |= (tree2->postorder_list[i].weight!=1);
   
   for (i=n1, l=2*n1-1; i>0; i--) {
      if (alignment[0][i]!=0) t1[l--]=']';
      else t1[l--]=')';
      p=i;
      while(i==tree1->postorder_list[p].leftmostleaf) {
	 if (alignment[0][p]!=0) t1[l--]='[';
	 else t1[l--]='(';
	 p=tree1->postorder_list[p].father;
      }
   }
   t1[2*n1]='\0';
   for (j=n2, l=2*n2-1; j>0; j--) {
      if (alignment[1][j]!=0) t2[l--]=']';
      else t2[l--]=')';
      p=j;
      while(j==tree2->postorder_list[p].leftmostleaf) {
	 if (alignment[1][p]!=0) t2[l--]='[';
	 else t2[l--]='(';
	 p=tree2->postorder_list[p].father;
      }
   }
   t2[2*n2]='\0';
   
   ni=nj=l=i=j=0;
   while (t1[i]||t2[j]) {
      while ((t1[i]=='(')||(t1[i]==')')) {
	 if (t1[i]==')') {
	    ni++;
	    encode(tree1->postorder_list[ni].type, ll);
	    if (weights)
	       sprintf(ll+strlen(ll), "%d", tree1->postorder_list[ni].weight);
	    for (k=0; k< strlen(ll); k++) {
	       a1[l]=ll[k];
	       a2[l++]='_';
	    }
	    a1[l]=')'; a2[l++]='_';
	 } else {
	    a1[l] = t1[i];
	    a2[l++] ='_';
	 }
	 i++;
      }
      
      while ((t2[j]=='(')||(t2[j]==')')) {
	 if (t2[j]==')') {
	    nj++;
	    encode(tree2->postorder_list[nj].type, ll);
	    if (weights)
	       sprintf(ll+strlen(ll), "%d", tree2->postorder_list[nj].weight);
	    for (k=0; k< strlen(ll); k++) {
	       a2[l]=ll[k];
	       a1[l++]='_';
	    }
	    a2[l]=')'; a1[l++]='_';
	 } else {
	    a2[l] = t2[j];
	    a1[l++] ='_';
	 }
	 j++;
      }
      
      if (t2[j]==']') {
	 ni++; nj++;
	 encode(tree2->postorder_list[nj].type, ll);
	 if (weights)
	    sprintf(ll+strlen(ll), "%d", tree2->postorder_list[nj].weight);
	 encode(tree1->postorder_list[ni].type, ll1);
	 if (weights)
	    sprintf(ll1+strlen(ll1), "%d", tree1->postorder_list[ni].weight);
	 if (strlen(ll)>strlen(ll1)) 
	    for (k=0; k<strlen(ll)-strlen(ll1); k++) strcat(ll1,"_");
	 if (strlen(ll)<strlen(ll1))
	    for (k=0; k<strlen(ll1)-strlen(ll); k++) strcat(ll,"_");
	 for (k=0; k< strlen(ll); k++) a2[l+k]=ll[k];
	 for (k=0; k< strlen(ll); k++) a1[l+k]=ll1[k];
	 l+=k;
	 a1[l]=a2[l]=')'; l++;
	 i++; j++;
      } else if (t2[j]=='[') {
	 a1[l]=a2[l]='('; l++;
	 i++; j++;
      }
   }
   a1[l]=a2[l]='\0';
   if (l>8*MNODES) nrerror("structure too long in sprint_aligned_trees");
   if (aligned_line[0]!= NULL)  free(aligned_line[0]); 
   if (aligned_line[1]!= NULL)  free(aligned_line[1]);
   aligned_line[0] = (char *) space((l+1)*sizeof(char));
   aligned_line[1] = (char *) space((l+1)*sizeof(char));
   strcpy(aligned_line[0], a1);
   strcpy(aligned_line[1], a2);
}

/*---------------------------------------------------------------------------*/
