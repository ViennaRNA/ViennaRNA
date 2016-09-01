
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "distance_matrix.h"
#include "ViennaRNA/utils.h"

#define PUBLIC
#define PRIVATE  static

#define MAXTAXA_FOR_LABELS  200

typedef struct{
        int   set1;
        int   set2;
        float distance;
        float distance2;
        } Union;

typedef struct _Node_ {
                 float  height;
                 float  brr;
                 float  brl;
                 int    whoami;
                 int    size;
                 struct _Node_  *father;
                 struct _Node_  *left;
                 struct _Node_  *right;                
               } Node;

PUBLIC void  PSplot_phylogeny(Union *cluster, char *filename, char *type);

PRIVATE Node *W2Phylo(Union *cluster);
PRIVATE Node *Nj2Phylo(Union *cluster);
PRIVATE void  free_phylo_tree(Node *root);
PRIVATE void fill_height_into_tree(Node *root);
PRIVATE void fill_br_from_height(Node *root);
PRIVATE void plot_branch(Node *root, FILE *fp);
PRIVATE void format_number(float x);

PRIVATE char  str[10];
PRIVATE float threshold;
PRIVATE int   print_labels=1;

/* --------------------------------------------------------------------------*/

PUBLIC void  PSplot_phylogeny(Union *cluster, char *filename, char *type)
{
   int   n;
   char  outfile[50];
   char  *tmp, *tfont;
   Node  *root;
   Node  *tempnode;
   FILE  *fp;
   float xsize, ysize, tfontsize, lfontsize, lwidth;
      
   n = cluster[0].set1;
   switch(type[0]) {
      case 'W' :
        root = W2Phylo(cluster);
        break;
      case 'N' : 
        root = Nj2Phylo(cluster);
        break;
      default :
        return;
   }

   outfile[0] = '\0';
   tmp = get_taxon_label(-1);   /* retrieve data set identifier */
   if(tmp) {strcat(outfile,tmp); strcat(outfile,"_"); free(tmp);}
   strcat(outfile,filename);
   
   fp = fopen(outfile,"w");
   if (fp==NULL) {
      fprintf(stderr,"couldn't open %s -- not doing treeplot\n", outfile);
      return;
   }
   xsize = root->size;
   tempnode = root;
   while(tempnode->right) 
      tempnode = tempnode->right;   
   ysize = tempnode->height;
   threshold = ysize*0.06;

   lwidth = MIN2(1.5, 8/sqrt((double)n));
   tfontsize = MIN2(15,550./(10+n));
   lfontsize = 2./3.*tfontsize;
   if (n>30) tfont = "Helvetica";
   else tfont = "Times-Roman";
   if(n > MAXTAXA_FOR_LABELS) print_labels=0;
   else print_labels=1;
   fprintf(fp,"%%!PS-Adobe-2.0 EPSF-1.2\n");
   fprintf(fp,"%%%%Title: TreePlot (%s)\n",type);
   fprintf(fp,"%%%%Creator: AnalyseDists\n");
   fprintf(fp,"%%%%CreationDate: %s", vrna_time_stamp());
   /* BoundingBox is only approximate */
   fprintf(fp,"%%%%BoundingBox: 35 45 535 640\n");
   fprintf(fp,"%%%%Pages: 1\n");
   fprintf(fp,"%%%%EndComments: No comment!\n");
   fprintf(fp,"288.5 50 translate\n");
   fprintf(fp,"%3.1f setlinewidth\n"
	      "/cmtx matrix currentmatrix def\n", lwidth);
   fprintf(fp,"500 %g div 360 %g div scale\n", xsize, ysize); 
   fprintf(fp,"/rotshow {gsave cmtx setmatrix\n"); 
   fprintf(fp,"          90 rotate 5 %4.1f rmoveto show grestore} def\n",
	   -tfontsize/4);
   fprintf(fp,"/cshow {gsave cmtx setmatrix\n"
              "          /Helvetica findfont %4.1f scalefont setfont\n"
              "          90 rotate 0 %3.1f rmoveto\n"
              "          dup stringwidth pop 2 div neg 0 rmoveto show\n"
              "        grestore} def\n", lfontsize, lfontsize/5);
   fprintf(fp,"/%s findfont %4.1f scalefont setfont\n", tfont, tfontsize);

   fprintf(fp,"0 0 moveto\n");
   plot_branch(root,fp);
   
   fprintf(fp,"cmtx setmatrix stroke\nshowpage\n");
   free_phylo_tree(root);
   fclose(fp);

}

/* --------------------------------------------------------------------------*/

PRIVATE Node *W2Phylo(Union *cluster)
{
   int i,n;
   float b;
   Node **taxa;
   Node  *father;
   
   n=cluster[0].set1;
   taxa = (Node **) vrna_alloc(sizeof(Node*)*(n+1));
   
   b = sqrt(MAX2(cluster[n-1].distance,0.));
   for(i=1;i<=n;i++){
      taxa[i] = (Node *) vrna_alloc(sizeof(Node));
      taxa[i]->whoami = i;
      taxa[i]->size   = 1;
      taxa[i]->height = b;
   }
   
   for(i=1;i<n;i++) {
      father = (Node *) vrna_alloc(sizeof(Node));
      father->whoami = 0;
      father->left  = taxa[cluster[i].set1];
      father->right = taxa[cluster[i].set2];
      father->size  = father->left->size+father->right->size;
      b = sqrt(MAX2(cluster[n-1].distance-cluster[i].distance,0.));
      father->height   = b;
      father->left->father   = father;
      father->right->father  = father;
      taxa[cluster[i].set1] = father;
   }
   free(taxa);
   father->whoami = -1;
   fill_br_from_height(father);
   return father;
}


/* --------------------------------------------------------------------------*/

PRIVATE void fill_br_from_height(Node *root)
{
   root->brr  = root->right->height - root->height;
   root->brl  = root->left->height  - root->height;
   if(!root->left->whoami)  fill_br_from_height(root->left);
   if(!root->right->whoami)  fill_br_from_height(root->right);
   return;
}

/* -------------------------------------------------------------------------- */
#define ORDER(X) \
if ((X)->brl>(X)->brr) { \
   xnode=(X)->right; \
   (X)->right=(X)->left; \
   (X)->left=xnode; \
   br = (X)->brr; \
   (X)->brr=(X)->brl; \
   (X)->brl=br; \
}



PRIVATE Node *Nj2Phylo(Union *cluster)
{
   int i,n, br;
   float h1,h2,maxdist,dist;
   Node **taxa;
   Node  *father;
   Node  *topnode;   /* of the longest path in the tree */
   Node  *tempnode, *xnode;
   Node  *root;

   maxdist = 0.;
   n=cluster[0].set1;
   taxa = (Node **) vrna_alloc(sizeof(Node*)*(n+1));
   
   for(i=1;i<=n;i++){
      taxa[i] = (Node *) vrna_alloc(sizeof(Node));
      taxa[i]->whoami = i;
      taxa[i]->size   = 1;
   }
 
   for(i=1;i<n;i++) {
      father = (Node *) vrna_alloc(sizeof(Node));
      father->whoami = 0;
      h1 = cluster[i].distance +taxa[cluster[i].set1]->height;
      h2 = cluster[i].distance2+taxa[cluster[i].set2]->height;  
      dist = h1 + h2;
      if(dist > maxdist) {
         maxdist = dist;
         topnode = father;
      }
      if(h1<=h2) {
         father->brr   = cluster[i].distance2;
         father->brl   = cluster[i].distance;
         father->left  = taxa[cluster[i].set1];
         father->right = taxa[cluster[i].set2];
      }
      else {
         father->brr   = cluster[i].distance;
         father->brl   = cluster[i].distance2;
         father->left  = taxa[cluster[i].set2];
         father->right = taxa[cluster[i].set1];
      }
      father->height        = MAX2(h1,h2);
      father->size          = father->left->size + father->right->size;
      father->left->father   = father;
      father->right->father  = father;
      taxa[cluster[i].set1] = father;
   }

   /* new root is inserted in the middle of the longest path */
   dist = maxdist/2.;
   tempnode = topnode;
   while (tempnode->height > dist)
      tempnode = tempnode->right;
   
   topnode  = tempnode->father;
   /* new root must go between topnode and tempnode */
   h1 = dist - tempnode->height;
   h2 = topnode->height - dist ;

   root = (Node *) vrna_alloc(sizeof(Node));
   root->father = NULL;
   root->whoami = -1;
   root->right  = tempnode;
   root->left   = topnode;
   root->height = 0.0;
   root->brr    = h1;
   root->brl    = h2;
   root->size   = n;

   topnode->size = n - tempnode->size;

   tempnode = root;
   while(topnode->father) {
      topnode->right  = topnode->father;
      topnode->brr    = topnode->father->brr;
      topnode->father = tempnode;
      tempnode = topnode;
      topnode  = topnode->right;
      topnode->size   = tempnode->size - tempnode->left->size;
      ORDER(tempnode);
   }
   /* topnode is now old root, remove it */
   if (tempnode->right==topnode) {
      tempnode->right = topnode->left;
      tempnode->brr   = topnode->brr + topnode->brl;
      tempnode->right->father = tempnode;
   } else {
      tempnode->left = topnode->left;
      tempnode->brl  = topnode->brr + topnode->brl;
      tempnode->left->father = tempnode;
   }
   free(topnode);
   ORDER(root);
   
   fill_height_into_tree(root);
   
   return root;
}

/* -------------------------------------------------------------------------- */

PRIVATE void fill_height_into_tree(Node *root)
{
   if(root->whoami>0) return;
   root->left->height  = root->height + root->brl;
   root->right->height = root->height + root->brr;
   fill_height_into_tree(root->left);
   fill_height_into_tree(root->right);
   return;
}
   
/* -------------------------------------------------------------------------- */

PRIVATE void free_phylo_tree(Node *root) {
   if(root->left)  free_phylo_tree(root->left);
   if(root->right) free_phylo_tree(root->right);
   free(root);
}  


/* -------------------------------------------------------------------------- */

PRIVATE void plot_branch(Node *root, FILE *fp)
{
   char *label;

   fprintf(fp,"currentpoint %g 0  rlineto 0 %g rlineto \n",
                     -(float)(root->left->size)/2., root->brr);
   if((print_labels)&&(root->brr > threshold)) {
      format_number(root->brr);
      fprintf(fp,"currentpoint 0 %g rlineto \n", -root->brr/2.);
      fprintf(fp,"(%s) cshow moveto\n", str);
   }
   if(root->right->whoami==0) plot_branch(root->right,fp);
   else if(print_labels) {
      label = get_taxon_label(root->right->whoami);
      fprintf(fp, "(%s) rotshow\n", label);
      free(label);
   }
   fprintf(fp,"moveto\n");
   fprintf(fp,"currentpoint %g 0  rlineto 0 %g rlineto \n",
                     +(float)(root->right->size)/2., root->brl);
   if((print_labels)&&(root->brl > threshold)) {
      format_number(root->brl);
      fprintf(fp,"currentpoint 0 %g rlineto \n", -root->brl/2.);
      fprintf(fp,"(%s) cshow moveto\n",str);
   }
   if(root->left->whoami==0) plot_branch(root->left,fp);
   else if(print_labels) {
      label = get_taxon_label(root->left->whoami);
      fprintf(fp, "(%s) rotshow\n", label);
      free(label);
   }
   fprintf(fp,"moveto\n");  
}

/* -------------------------------------------------------------------------- */

PRIVATE void format_number(float x)
{
   if(x >= 1.e5) { sprintf(str, "%4.1g",x ); return;}
   if(x >=  100) { sprintf(str, "%.0f" ,x ); return;}
   if(x >=   10) { sprintf(str, "%.1f", x ); return;}
   if(x >=    1) { sprintf(str, "%.2f", x ); return;}
   if(x >=  0.1) { sprintf(str, "%.2f", x ); return;}
   if(x >= 0.01) { sprintf(str, "%.3f", x ); return;}
   sprintf(str, "%4.1g",x ); 
   return;
}
   
/* -------------------------------------------------------------------------- */
