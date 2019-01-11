#include "ViennaRNA/RNApuzzler/output/configtree_debug.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/definitions.h"

#include "ViennaRNA/utils.h"

#include <stdio.h>
#include <stdlib.h>

int indexNR = -1;

void setRGB(
        double rgb[3],
        const RGB_SELECTOR color
) {
    switch (color) {
    case RGB_BLACK:
        rgb[0] =   0; rgb[1] =   0; rgb[2] =   0;
        break;
    case RGB_GREY:
        rgb[0] = 189; rgb[1] = 189; rgb[2] = 189;
        break;
    case RGB_RED:
        rgb[0] = 228; rgb[1] =  26; rgb[2] =  28;
        break;
    case RGB_BLUE:
        rgb[0] =  55; rgb[1] = 126; rgb[2] = 184;
        break;
    case RGB_GREEN:
        rgb[0] =  77; rgb[1] = 175; rgb[2] =  74;
        break;
    case RGB_ORANGE:
        rgb[0] = 255; rgb[1] = 127; rgb[2] =   0;
        break;
    case RGB_PURPLE:
        rgb[0] = 152; rgb[1] =  78; rgb[2] = 163;
        break;
    }
}

void printNode(treeNode* node, FILE* fp) {
    stemBox* sbox = node->sBox;
    loopBox* lbox = node->lBox;
    double c[2] = { lbox->c[0], lbox->c[1] };
    double p1[2] = { sbox->c[0] + sbox->e[0] * sbox->a[0] + sbox->e[1] * sbox->b[0], sbox->c[1] + sbox->e[0] * sbox->a[1] + sbox->e[1] * sbox->b[1]};
    double p2[2] = { sbox->c[0] + sbox->e[0] * sbox->a[0] - sbox->e[1] * sbox->b[0], sbox->c[1] + sbox->e[0] * sbox->a[1] - sbox->e[1] * sbox->b[1]};
    double p3[2] = { sbox->c[0] - sbox->e[0] * sbox->a[0] + sbox->e[1] * sbox->b[0], sbox->c[1] - sbox->e[0] * sbox->a[1] + sbox->e[1] * sbox->b[1]};
    double p4[2] = { sbox->c[0] - sbox->e[0] * sbox->a[0] - sbox->e[1] * sbox->b[0], sbox->c[1] - sbox->e[0] * sbox->a[1] - sbox->e[1] * sbox->b[1]};
    double r     = lbox->r;
    fprintf(fp
            ,"[[%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f]%3.2f]\n"
            ,c[0],c[1],p1[0],p1[1],p3[0],p3[1],p4[0],p4[1],p2[0],p2[1],r
           );
}

void printBulges(treeNode* node, FILE* fp) {
    stemBox* stem = node->sBox;
    if (stem->bulgeCount > 0) {
        int i;
        for (i = 0; i < stem->bulgeCount; i++) {
            double A[2], B[2], C[2];
            getBulgeCoordinates(stem, i, A, B, C);
            fprintf(fp, ""
                        "[[%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f]]\n"
                        "", A[0], A[1], B[0], B[1], C[0], C[1]);
        }
    }
}

void printBulgesDebugToFile(treeNode* tree, FILE* fp) {

    printBulges(tree, fp);

    int i;
    for (i = 0; i < tree->childCount; i++) {
        treeNode* child = getChild(tree, i);
        printBulgesDebugToFile(child, fp);
    }
}

/*
 * function printBoxesDebug2ToFile
 */
void printBoxesDebug2ToFile(treeNode* tree, FILE* fp) {

    printNode(tree, fp);

    int i;
    for (i = 0; i < tree->childCount; i++) {
        treeNode* child = getChild(tree, i);
        printBoxesDebug2ToFile(child, fp);
    }
}

void head(FILE* fp, const int changeNumber) {
    // ------- Head -------

    fprintf(fp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(fp, "%%%%Creator: ViennaRNA-2.3.5\n");
    fprintf(fp, "%%%%CreationDate: Thu Sep  7 15:59:36 2017\n");
    fprintf(fp, "%%%%Title: RNA Secondary Structure Plot\n");
    fprintf(fp, "%%%%BoundingBox: 0 0 700 700\n");
    fprintf(fp, "%%%%DocumentFonts: Helvetica\n");
    fprintf(fp, "%%%%Pages: 1\n");
    fprintf(fp, "%%%%EndComments\n");

    fprintf(fp, ""
                "\n"
                "%% Bounding Boxes Step %d\n"
                "\n", changeNumber);

}

void doSetBounds(FILE* fp) {
    // ------- Set Bounds -------

    fprintf(fp, ""
                "\n"
                "1 setlinejoin\n"
                "1 setlinecap\n"
                "0.8 setlinewidth\n"
                "72 216 translate\n"
                "\n"
                "%% find the coordinate range\n"
                "/xmax -1000000 def /xmin 1000000 def\n"
                "/ymax -1000000 def /ymin 1000000 def\n"
                "/i 0 def\n"
                "boxes {\n"
                "    /box boxes i get def\n"
                "    /a box 0 get def\n"
                "    /b box 1 get def\n"
                "    /c box 2 get def\n"
                "    /d box 3 get def\n"
                "    /e box 4 get def\n"
                "    /r box 5 get def\n"
                "    \n"
                "    /x a 0 get r sub def\n"
                "    /y a 1 get r sub def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "    /x a 0 get r add def\n"
                "    /y a 1 get r add def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "    /x b 0 get def\n"
                "    /y b 1 get def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "    /x c 0 get def\n"
                "    /y c 1 get def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "    /x d 0 get def\n"
                "    /y d 1 get def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "    /x e 0 get def\n"
                "    /y e 1 get def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "    \n"
                "      /i i 1 add def\n"
                "} forall\n"
                "/size {xmax xmin sub ymax ymin sub max} bind def\n"
                "72 6 mul size div dup scale\n"
                "size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div\n"
                "translate\n"
                "\n");
}

void doPrintGrid(FILE* fp) {
    // ------- Print Grid -------
    fprintf(fp, ""
                "%% draw vertical line from 0 to xmin\n"
                "newpath\n"
                "/i 0 def\n"
                "i {\n"
                "   i xmin lt {\n"
                "       exit\n"
                "   } {\n"
                "       i ymin moveto i 0 lineto i ymax lineto\n"
                "       /i i 100 sub def\n"
                "   } ifelse\n"
                "   "
                "} loop\n"
                "\n"
                "%% draw vertical line from 0 to xmax\n"
                "/i 0 def\n"
                "i {\n"
                "   i xmax gt {\n"
                "       exit\n"
                "   } {\n"
                "       i ymin moveto i 0 lineto i ymax lineto\n"
                "       /i i 100 add def\n"
                "   } ifelse\n"
                "   "
                "} loop\n"
                "\n"
                "%% draw horizontal line from 0 to ymin\n"
                "/i 0 def\n"
                "i {\n"
                "   i ymin lt {\n"
                "       exit\n"
                "   } {\n"
                "       xmin i moveto 0 i lineto xmax i lineto\n"
                "       /i i 100 sub def\n"
                "   } ifelse\n"
                "   "
                "} loop\n"
                "\n"
                "%% draw horizontal line from 0 to ymax\n"
                "/i 0 def\n"
                "i {\n"
                "   i ymax gt {\n"
                "       exit\n"
                "   } {\n"
                "       xmin i moveto 0 i lineto xmax i lineto\n"
                "       /i i 100 add def\n"
                "   } ifelse\n"
                "   "
                "} loop\n"
                "1 setlinewidth\n"
                "0.5 0.5 0.5 setrgbcolor\n"
                "stroke\n"
                "\n"
                "%% draw coordinate axes\n"
                "newpath\n"
                "2 setlinewidth\n"
                "xmin 0 moveto 0 0 lineto xmax 0 lineto\n"
                "0 ymin moveto 0 0 lineto 0 ymax lineto\n"
//                "0.1 0.1 0.1 setrgbcolor\n"
                "0 0 1 setrgbcolor\n"
                "stroke\n"
                "\n");
}

void doPrintBoxesAllBlack(FILE* fp) {
    // ------- Print Boxes -------

    fprintf(fp, ""
                "/i 0 def\n"
                "boxes {\n"
                "    newpath\n"
                "    /box boxes i get def\n"
                "    /a box 0 get def\n"
                "    /b box 1 get def\n"
                "    /c box 2 get def\n"
                "    /d box 3 get def\n"
                "    /e box 4 get def\n"
                "    /r box 5 get def\n"
                "    a aload pop moveto\n"
                "    b aload pop lineto\n"
                "    c aload pop lineto\n"
                "    d aload pop lineto\n"
                "    e aload pop lineto\n"
                "    a aload pop lineto\n"
                "    /xr a 0 get r add def\n"
                "    /yr a 1 get def\n"
                "    xr yr moveto\n"
                "    a aload pop r 0.000001 0 arc\n"
                "    a aload pop moveto\n"
                "    i 0 eq { 0 0 1 setrgbcolor }\n"
                "           { 0 0 0 setrgbcolor } ifelse\n"
                "    2 setlinewidth\n"
                "    stroke\n"
                "    /i i 1 add def\n"
                "} forall\n"
                "\n");
}

void doPrintBoxesFirstBlueLastRed(FILE* fp) {
    // ------- Print Boxes -------

    fprintf(fp, ""
                "newpath\n"
                "/i 0 def\n"
                "boxes {\n"
                "    stroke\n"
                "    newpath\n"
                "    2 setlinewidth\n"
                "    i 0 eq {\n"
                "        0 0 1 setrgbcolor\n"
                "    }{\n"
                "        0 0 0 setrgbcolor\n"
                "    } ifelse\n"
                "    /box boxes i get def\n"
                "    /a box 0 get def\n"
                "    /b box 1 get def\n"
                "    /c box 2 get def\n"
                "    /d box 3 get def\n"
                "    /e box 4 get def\n"
                "    /r box 5 get def\n"
                "    a aload pop moveto\n"
                "    b aload pop lineto\n"
                "    c aload pop lineto\n"
                "    d aload pop lineto\n"
                "    e aload pop lineto\n"
                "    a aload pop lineto\n"
                "    /xr a 0 get r add def\n"
                "    /yr a 1 get def\n"
                "    xr yr moveto\n"
                "    a aload pop r 0.000001 0 arc\n"
                "    a aload pop moveto\n"
                "    /i i 1 add def\n"
                "} forall\n"
                "1 0 0 setrgbcolor\n"
                "stroke\n"
                "\n");
}

void doPrintBulges(FILE* fp) {
    // ------- Print Bulges -------

    fprintf(fp, ""
                "/i 0 def\n"
                "bulges {\n"
                "        newpath\n"
                "        /box bulges i get def\n"
                "        /j 0 def\n"
                "        /first 1 def\n"
                "        box {\n"
                "                /p box j get def\n"
                "                j 0 eq { p aload pop moveto }\n"
                "                       { p aload pop lineto } ifelse\n"
                "                /j j 1 add def\n"
                "        } forall\n"
                "        0 0 0 setrgbcolor\n"
                "        2 setlinewidth\n"
                "        stroke\n"
                "    /i i 1 add def\n"
                "} forall\n"
                "\n");
}

void doPrintExteriorIntersectionLine(FILE* fp) {
    // ------- Tail -------

    fprintf(fp, ""
                "\n"
                "newpath\n"
                "ymin %5.1f lt { 1 0 0 setrgbcolor }\n"
                "              { 0 1 0 setrgbcolor } ifelse\n"
                "xmin ymin moveto\n"
                "xmax ymin lineto\n"
                "stroke\n"
                "\n"
            , EXTERIOR_Y);
}

void doShowPage(FILE* fp) {
    fprintf(fp, ""
                "showpage\n");
}

void tailTree(FILE* fp) {
    doSetBounds(fp);
    //doPrintGrid(fp);
    doPrintBoxesAllBlack(fp);
    doPrintBulges(fp);
    doPrintExteriorIntersectionLine(fp);
    doShowPage(fp);
}

void tailPath(FILE* fp) {
    doSetBounds(fp);
    //doPrintGrid(fp);
    doPrintBoxesFirstBlueLastRed(fp);
    doPrintBulges(fp);
    doPrintExteriorIntersectionLine(fp);
    doShowPage(fp);
}

/*
 * function printBoxesDebugToFile
 */
void printBoxesDebugToFile(
        treeNode* tree,
        const puzzlerOptions* puzzler,
        const int changeNumber
) {
    if (indexNR == changeNumber) {
        return;
    } else {
        indexNR = changeNumber;
    }

    char filename[1000];
    if (puzzler->filename != NULL) {
        snprintf(filename, sizeof(filename), "%s_DEBUG_%05d.ps", puzzler->filename, changeNumber);
    } else {
        snprintf(filename, sizeof(filename), "BoundingBoxes_%05d.ps", changeNumber);
    }

    FILE* fp;
    fp = fopen(filename, "w");

    head(fp, changeNumber);

    // ------- Boxes -------

    fprintf(fp,
            "/boxes [\n");
    if (isExterior(tree)) {
        // root has no bounding box... so skip this one in printing
        int i;
        int childCount = tree->childCount;
        for (i = 0; i < childCount; i++) {
            treeNode* child = getChild(tree, i);
            printBoxesDebug2ToFile(child, fp);
        }
    } else {
        printBoxesDebug2ToFile(tree, fp);
    }
    fprintf(fp, ""
                "] def \n"
                "");

    // ------- Bulges -------

    fprintf(fp,
            "/bulges [\n");
    if (isExterior(tree)) {
        // root has no bounding box... so skip this one in printing
        int i;
        int childCount = tree->childCount;
        for (i = 0; i < childCount; i++) {
            treeNode* child = getChild(tree, i);
            printBulgesDebugToFile(child, fp);
        }
    } else {
        printBulgesDebugToFile(tree, fp);
    }
    fprintf(fp, ""
                "] def \n"
                "");

    tailTree(fp);
    fclose(fp);
}

/*
 * function printBoxesDebug2ToFile
 */
void printBoxesDebug2(treeNode* tree) {
    stemBox* sbox = tree->sBox;
    loopBox* lbox = tree->lBox;
    double c[2] = { lbox->c[0], lbox->c[1] };
    double p1[2] = { sbox->c[0] + sbox->e[0] * sbox->a[0] + sbox->e[1] * sbox->b[0], sbox->c[1] + sbox->e[0] * sbox->a[1] + sbox->e[1] * sbox->b[1]};
    double p2[2] = { sbox->c[0] + sbox->e[0] * sbox->a[0] - sbox->e[1] * sbox->b[0], sbox->c[1] + sbox->e[0] * sbox->a[1] - sbox->e[1] * sbox->b[1]};
    double p3[2] = { sbox->c[0] - sbox->e[0] * sbox->a[0] + sbox->e[1] * sbox->b[0], sbox->c[1] - sbox->e[0] * sbox->a[1] + sbox->e[1] * sbox->b[1]};
    double p4[2] = { sbox->c[0] - sbox->e[0] * sbox->a[0] - sbox->e[1] * sbox->b[0], sbox->c[1] - sbox->e[0] * sbox->a[1] - sbox->e[1] * sbox->b[1]};
    double r     = lbox->r;
    printf("[[%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f]%3.2f]\n"
           ,c[0],c[1],p1[0],p1[1],p3[0],p3[1],p4[0],p4[1],p2[0],p2[1],r
           );

    int childCount = tree->childCount;
    int i;
    for (i = 0; i < childCount; i++) {
        printBoxesDebug2(getChild(tree, i));
    }
}

/**
 * @brief printBoxesDebug
 *      - Debug method that allows to show bounding boxes in PS output.
 *        Copy the resulting terminal output to the end of the PS-file to see bounding boxes.
 * @param tree
 *      - root node
 */
void printBoxesDebug(treeNode* tree) {
    printf("\n");
    printf("[DEBUG] printBoxesDebug - start\n");
    printf("/boxes [\n");

    if (tree->parent == NULL) {
        // root has no bounding box... so skip this one in printing
        int i;
        int childCount = tree->childCount;
        for (i = 0; i < childCount; i++) {
            printBoxesDebug2(getChild(tree, i));
        }
    } else {
        printBoxesDebug2(tree);
    }

    printf(""
            "] def \n"
            "newpath\n"
            "/i 0 def\n"
            "boxes {\n"
            "    /box boxes i get def\n"
            "    /a box 0 get def\n"
            "    /b box 1 get def\n"
            "    /c box 2 get def\n"
            "    /d box 3 get def\n"
            "    /e box 4 get def\n"
            "    /r box 5 get def\n"
            "    a aload pop moveto\n"
            "    b aload pop lineto\n"
            "    c aload pop lineto\n"
            "    d aload pop lineto\n"
            "    e aload pop lineto\n"
            "    a aload pop lineto\n"
            "    a aload pop r 0.000001 0 arc\n"
            "    /i i 1 add def\n"
            "} forall\n"
            "2 setlinewidth\n"
            "stroke\n");
    printf("\n");
    printf("[DEBUG] printBoxesDebug - end\n");
}

void PS_printTree(
        treeNode* root,
        const puzzlerOptions* puzzler
) {

    char filename[1000];
    if (puzzler->filename != NULL) {
        snprintf(filename, sizeof(filename), "%s_DEBUG_TREE_%05d_%04d.ps", puzzler->filename, puzzler->numberOfChangesAppliedToConfig, getNodeID(root));
    } else {
        snprintf(filename, sizeof(filename), "BoundingBoxes_Tree_%05d.ps", puzzler->numberOfChangesAppliedToConfig);
    }

    FILE* fp;
    fp = fopen(filename, "w");

    head(fp, puzzler->numberOfChangesAppliedToConfig);

    // ------- Boxes -------

    fprintf(fp,
            "/boxes [\n");
    if (isExterior(root)) {
        // root has no bounding box... so skip this one in printing
        int i;
        int childCount = root->childCount;
        for (i = 0; i < childCount; i++) {
            treeNode* child = getChild(root, i);
            printBoxesDebug2ToFile(child, fp);
        }
    } else {
        printBoxesDebug2ToFile(root, fp);
    }
    fprintf(fp, ""
                "] def \n"
                "");

    // ------- Bulges -------

    fprintf(fp,
            "/bulges [\n");
    if (isExterior(root)) {
        // root has no bounding box... so skip this one in printing
        int i;
        int childCount = root->childCount;
        for (i = 0; i < childCount; i++) {
            treeNode* child = getChild(root, i);
            printBulgesDebugToFile(child, fp);
        }
    } else {
        printBulgesDebugToFile(root, fp);
    }
    fprintf(fp, ""
                "] def \n"
                "");

    tailTree(fp);
    fclose(fp);
}

void PS_printPath(
        treeNode* ancestor,
        treeNode* intersector,
        const puzzlerOptions* puzzler
) {

    int ancestorID = getNodeID(ancestor);
    int intersectorID = getNodeID(intersector);

    char filename[1000];
    if (puzzler->filename != NULL) {
        snprintf(filename, sizeof(filename), "%s_DEBUG_PATH_%05d_%04d_vs_%04d.ps", puzzler->filename, puzzler->numberOfChangesAppliedToConfig, ancestorID, intersectorID);
    } else {
        snprintf(filename, sizeof(filename), "BoundingBoxes_Path_%05d_%04d_vs_%04d.ps", puzzler->numberOfChangesAppliedToConfig, ancestorID, intersectorID);
    }

    FILE* fp;
    fp = fopen(filename, "w");

    head(fp, puzzler->numberOfChangesAppliedToConfig);

    // ------- Head -------
    fprintf(fp, ""
                "\n"
                "%% Bounding Boxes for Path %c(%d) to %c(%d) at Step %d\n"
                "\n"
            , getNodeName(ancestor), getNodeID(ancestor)
            , getNodeName(intersector), getNodeID(intersector)
            , puzzler->numberOfChangesAppliedToConfig);

    int pathLength = 1;
    treeNode* node = intersector;
    while (node != ancestor) {
        node = getParent(node);
        ++pathLength;
    }
    treeNode** path = (treeNode**) vrna_alloc(pathLength * sizeof(treeNode*));
    node = intersector;
    for (int i = pathLength - 1; i >= 0; i--) {
        path[i] = node;
        node = getParent(node);
    }
    int* childIndex = (int*) vrna_alloc((pathLength-1) * sizeof(int));
    for (int i = 0; i < pathLength - 1; i++) {
        childIndex[i] = getChildIndex(path[i], getNodeID(path[i+1]));
    }

    fprintf(fp, "%%[PATH]");
    for (int i = 0; i < pathLength; i++) {
        fprintf(fp, " %c(%d/%d)", getNodeName(path[i]), (i < pathLength - 1 ? childIndex[i] + 1 : -1), path[i]->childCount);
    }
    fprintf(fp, "\n\n");

    // ------- Boxes -------
    fprintf(fp,
            "/boxes [\n");
    // print nodes on path
    for (int i = 0; i < pathLength; i++) {
        printNode(path[i], fp);
    }
    fprintf(fp, ""
                "] def \n"
                "");

    // ------- Bulges -------
    fprintf(fp,
            "/bulges [\n");
    // print bulges on path
    for (int i = 0; i < pathLength; i++) {
        printBulges(path[i], fp);
    }
    fprintf(fp, ""
                "] def \n"
                "");

    tailPath(fp);
    fclose(fp);

    free(path);
    free(childIndex);
}

void PS_printColoredBox(
        FILE* fp,
        treeNode* node,
        const double rgb[3]
) {
    fprintf(fp, "[");

    // color triple
    fprintf(fp, " [%f %f %f]", rgb[0], rgb[1], rgb[2]);

    loopBox* loop = node->lBox;
    // circle
    fprintf(fp, " [[%f %f]%f]", loop->c[0], loop->c[1], loop->r);

    stemBox* stem = node->sBox;
    // stem
    double p1[2] = { stem->c[0] + stem->e[0] * stem->a[0] + stem->e[1] * stem->b[0], stem->c[1] + stem->e[0] * stem->a[1] + stem->e[1] * stem->b[1]};
    double p2[2] = { stem->c[0] + stem->e[0] * stem->a[0] - stem->e[1] * stem->b[0], stem->c[1] + stem->e[0] * stem->a[1] - stem->e[1] * stem->b[1]};
    double p3[2] = { stem->c[0] - stem->e[0] * stem->a[0] + stem->e[1] * stem->b[0], stem->c[1] - stem->e[0] * stem->a[1] + stem->e[1] * stem->b[1]};
    double p4[2] = { stem->c[0] - stem->e[0] * stem->a[0] - stem->e[1] * stem->b[0], stem->c[1] - stem->e[0] * stem->a[1] - stem->e[1] * stem->b[1]};
    fprintf(fp
            ," [[%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f]]"
            ,p1[0],p1[1],p3[0],p3[1],p4[0],p4[1],p2[0],p2[1]);

    // bulges
    fprintf(fp, " [");
    if (stem->bulgeCount > 0) {
        int i;
        for (i = 0; i < stem->bulgeCount; i++) {
            double A[2], B[2], C[2];
            getBulgeCoordinates(stem, i, A, B, C);
            fprintf(fp, ""
                        "[[%3.2f %3.2f][%3.2f %3.2f][%3.2f %3.2f]]"
                        "", A[0], A[1], B[0], B[1], C[0], C[1]);
        }
    }
    fprintf(fp, "]");

    fprintf(fp, " ]\n");
}

void PS_printFancyTreeRec(
        FILE* fp,
        treeNode* currentNode,
        treeNode* featuredNode,
        const short nodeType_in
) {

    const short typeAncestor = 0;
    const short typeFeatured = 1;
    const short typeBackground = 2;
    const short typeFeaturedSubtree = 3;

    short nodeType = nodeType_in;
    if (currentNode == featuredNode) {
        nodeType = typeFeatured;
    }

    if (currentNode->sBox && currentNode->lBox) {
        double rgb[3];

        if (nodeType == typeAncestor) {
            setRGB(rgb, RGB_ORANGE);
        } else
        if (nodeType == typeFeatured) {
            setRGB(rgb, RGB_PURPLE);
        } else
        if (nodeType == typeBackground) {
            setRGB(rgb, RGB_GREY);
        } else
        if (nodeType == typeFeaturedSubtree) {
            setRGB(rgb, RGB_GREEN);
        } else {
            setRGB(rgb, RGB_BLACK);
        }

        rgb[0] /= 255.0;
        rgb[1] /= 255.0;
        rgb[2] /= 255.0;

        PS_printColoredBox(fp, currentNode, rgb);
    }

    int indexFeatured = -1;
    if (nodeType == typeAncestor) {
        indexFeatured = getChildIndex(currentNode, getNodeID(featuredNode));
    }

    // recursion
    for (int i = 0; i < currentNode->childCount; i++) {
        treeNode* child = getChild(currentNode, i);

        // determine childType
        int childType = typeBackground;
        if (nodeType == typeAncestor) {
            if (i == indexFeatured) {
                if (child == featuredNode) {
                    childType = typeFeatured;
                } else {
                    childType = typeAncestor;
                }
            } else {
                childType = typeBackground;
            }
        }
        if (nodeType == typeFeatured || nodeType == typeFeaturedSubtree) {
            childType = typeFeaturedSubtree;
        }
        if (nodeType == typeBackground) {
            childType = typeBackground;
        }

        // recursive call
        PS_printFancyTreeRec(fp, child, featuredNode, childType);
    }

}

void PS_printFancySiblingsRec(
        FILE* fp,
        treeNode* currentNode,
        treeNode* featuredNode,
        treeNode* left,
        treeNode* right,
        const short nodeType_in
) {

    const short typeAncestor = 0;
    const short typeFeatured = 1;
    const short typeBackground = 2;
    const short typeFeaturedSubtree = 3;
    const short typeLeftSubtree = 4;
    const short typeRightSubtree = 5;

    short nodeType = nodeType_in;
    if (currentNode == featuredNode) {
        nodeType = typeFeatured;
    }
    if (currentNode == left) {
        nodeType = typeLeftSubtree;
    }
    if (currentNode == right) {
        nodeType = typeRightSubtree;
    }

    if (currentNode->sBox && currentNode->lBox) {
        double rgb[3];

        if (nodeType == typeAncestor) {
            setRGB(rgb, RGB_ORANGE);
        } else
        if (nodeType == typeFeatured) {
            setRGB(rgb, RGB_PURPLE);
        } else
        if (nodeType == typeBackground) {
            setRGB(rgb, RGB_GREY);
        } else
        if (nodeType == typeFeaturedSubtree) {
            setRGB(rgb, RGB_GREEN);
        } else
        if (nodeType == typeRightSubtree){
            setRGB(rgb, RGB_RED);
        } else
        if (nodeType == typeLeftSubtree){
            setRGB(rgb, RGB_BLUE);
        } else {
            setRGB(rgb, RGB_BLACK);
        }

        rgb[0] /= 255.0;
        rgb[1] /= 255.0;
        rgb[2] /= 255.0;

        PS_printColoredBox(fp, currentNode, rgb);
    }

    int indexFeatured = -1;
    if (nodeType == typeAncestor) {
        indexFeatured = getChildIndex(currentNode, getNodeID(featuredNode));
    }

    // recursion
    for (int i = 0; i < currentNode->childCount; i++) {
        treeNode* child = getChild(currentNode, i);

        // determine childType
        int childType = typeBackground;
        if (nodeType == typeAncestor) {
            if (i == indexFeatured) {
                if (child == featuredNode) {
                    childType = typeFeatured;
                } else {
                    childType = typeAncestor;
                }
            } else {
                childType = typeBackground;
            }
        }
        if (nodeType == typeFeatured || nodeType == typeFeaturedSubtree) {
            childType = typeFeaturedSubtree;
        }
        if (nodeType == typeBackground) {
            childType = typeBackground;
        }
        if (nodeType == typeLeftSubtree) {
            childType = typeLeftSubtree;
        }
        if (nodeType == typeRightSubtree) {
            childType = typeRightSubtree;
        }

        // recursive call
        PS_printFancySiblingsRec(fp, child, featuredNode, left, right, childType);
    }

}

void tailColoredBoxes(
        FILE *fp
) {
    fprintf(fp, ""
                "1 setlinejoin\n"
                "1 setlinecap\n"
                "2.5 setlinewidth\n"
                "72 216 translate\n"
                "\n"
                "%% find the coordinate range\n"
                "/xmax -1000000 def /xmin 1000000 def\n"
                "/ymax -1000000 def /ymin 1000000 def\n"
                "/iBox 0 def\n"
                "boxes {\n"
                "    /box boxes iBox get def\n"
                "\n"
                "    %% eval circle\n"
                "    /loop box 1 get def\n"
                "    /loopC loop 0 get def\n"
                "    /loopR loop 1 get def\n"
                "    /x loopC 0 get loopR sub def\n"
                "    /y loopC 1 get loopR sub def\n"
                "    y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "\n"
                "    %% eval stem\n"
                "    /stem box 2 get def\n"
                "    /iStem 0 def\n"
                "    stem {\n"
                "      /p stem iStem get def\n"
                "      /x p 0 get def\n"
                "      /y p 1 get def\n"
                "      y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "\n"
                "      /iStem iStem 1 add def\n"
                "    } forall\n"
                "\n"
                "    %% eval bulges\n"
                "    /bulges box 3 get def\n"
                "    /iBulge 0 def\n"
                "    bulges {\n"
                "      /bulge bulges iBulge get def\n"
                "\n"
                "      /j 0 def\n"
                "      bulge {\n"
                "        /p bulge j get def\n"
                "        /x p 0 get def\n"
                "        /y p 1 get def\n"
                "        y ymin lt { /ymin y def } if y ymax gt { /ymax y def } if x xmin lt { /xmin x def } if x xmax gt { /xmax x def } if\n"
                "\n"
                "        /j j 1 add def\n"
                "      } forall\n"
                "\n"
                "      /iBulge iBulge 1 add\n"
                "    } forall\n"
                "\n"
                "    /iBox iBox 1 add def\n"
                "} forall\n"
                "/size {xmax xmin sub ymax ymin sub max} bind def\n"
                "72 6 mul size div dup scale\n"
                "size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div\n"
                "translate\n"
                "\n"
                "%% ------------------\n"
                "%% --- draw boxes ---\n"
                "%% ------------------\n"
                "/iBox 0 def\n"
                "boxes {\n"
                "    /box boxes iBox get def\n"
                "\n"
                "    newpath\n"
                "\n"
                "    %% eval color\n"
                "    /color box 0 get def\n"
                "    color aload pop setrgbcolor\n"
                "\n"
                "    %% eval circle\n"
                "    /loop box 1 get def\n"
                "    /loopC loop 0 get def\n"
                "    /loopR loop 1 get def\n"
                "    loopC aload pop loopR 0.000001 0 arc\n"
                "\n"
                "    %% eval stem\n"
                "    /stem box 2 get def\n"
                "    /iStem 0 def\n"
                "    stem {\n"
                "      /p stem iStem get def\n"
                "      iStem 0 eq { p aload pop moveto }\n"
                "                 { p aload pop lineto } ifelse\n"
                "\n"
                "      /iStem iStem 1 add def\n"
                "    } forall\n"
                "\n"
                "    %% eval bulges\n"
                "    /bulges box 3 get def\n"
                "    /iBulge 0 def\n"
                "    bulges {\n"
                "      /bulge bulges iBulge get def\n"
                "\n"
                "      /iP 0 def\n"
                "      bulge {\n"
                "        /p bulge iP get def\n"
                "        iP 0 eq { p aload pop moveto }\n"
                "                { p aload pop lineto } ifelse\n"
                "\n"
                "        /iP iP 1 add def\n"
                "      } forall\n"
                "\n"
                "      /iBulge iBulge 1 add\n"
                "    } forall\n"
                "\n"
                "    stroke\n"
                "\n"
                "    /iBox iBox 1 add def\n"
                "} forall\n"
                "\n"
                "showpage\n"
            );

}

void PS_printFancyPathRec(
        FILE *fp,
        treeNode *current,
        treeNode *start,
        treeNode *end,
        treeNode *feature,
        short nodeType_in
) {
    const short typeBeforeStart = 0;
    const short typeStart = 1;
    const short typePath = 2;
    const short typeEnd = 3;
    const short typeFeature = 4;
    const short typeBackground = 5;

    short nodeType = nodeType_in;
    if (current == start) {
        nodeType = typeStart;
    }
    if (current == feature) {
        nodeType = typeFeature;
    }
    if (current == end) {
        nodeType = typeEnd;
    }

    if (current->sBox && current->lBox) {
        double rgb[3];

        if (nodeType == typeBeforeStart) {
            setRGB(rgb, RGB_ORANGE);
        } else
        if (nodeType == typeBackground) {
            setRGB(rgb, RGB_GREY);
        } else
        if (nodeType == typeStart) {
            setRGB(rgb, RGB_BLUE);
        } else
        if (nodeType == typePath) {
            setRGB(rgb, RGB_ORANGE);
        } else
        if (nodeType == typeEnd) {
            setRGB(rgb, RGB_RED);
        } else
        if (nodeType == typeFeature) {
            setRGB(rgb, RGB_PURPLE);
        } else {
            setRGB(rgb, RGB_BLACK);
        }

        rgb[0] /= 255.0;
        rgb[1] /= 255.0;
        rgb[2] /= 255.0;

        PS_printColoredBox(fp, current, rgb);
    }

    int indexEnd = -1;
    if (nodeType == typeStart || nodeType == typePath || nodeType == typeFeature || nodeType == typeBeforeStart) {
        indexEnd = getChildIndex(current, getNodeID(end));
    }

    // recursion
    for (int i = 0; i < current->childCount; i++) {
        treeNode* child = getChild(current, i);

        // determine childType
        short childType = -1;
        if (nodeType == typeBeforeStart) {
            childType = typeBeforeStart;
        }
        if (nodeType == typeBackground) {
            childType = typeBackground;
        }
        if (nodeType == typeStart || nodeType == typePath || nodeType == typeFeature || nodeType == typeBeforeStart) {
            if (i == indexEnd) {
                if (nodeType == typeBeforeStart) {
                    childType = typeBeforeStart;
                }
                if (nodeType == typeStart || nodeType == typeFeature || nodeType == typePath) {
                    childType = typePath;
                }
            } else {
                childType = typeBackground;
            }
        }
        if (nodeType == typeEnd || nodeType == typeBackground) {
            childType = typeBackground;
        }

        // recursive call
        PS_printFancyPathRec(fp, child, start, end, feature, childType);
    }

}

void PS_printFancyPath(
        treeNode* start,
        treeNode* end,
        treeNode* feature,
        puzzlerOptions* puzzler
) {
    (puzzler->psNumber)++;

    char filename[1000];
    snprintf(filename, sizeof(filename), "%s_%06d_FANCY_PATH_%05d_%04d_%04d_%04d.ps", puzzler->filename, puzzler->psNumber, puzzler->numberOfChangesAppliedToConfig, getNodeID(start), getNodeID(end), getNodeID(feature));

    FILE* fp;
    fp = fopen(filename, "w");

    head(fp, puzzler->numberOfChangesAppliedToConfig);

    // ------- Head -------
//    fprintf(fp, ""
//                "\n"
//                "%% Bounding Boxes for Path %c(%d) to %c(%d) at Step %d\n"
//                "\n"
//            , getNodeName(start), getNodeID(start)
//            , getNodeName(end), getNodeID(end)
//            , puzzler->numberOfChangesAppliedToConfig);

//    int pathLength = 1;
//    treeNode* node = end;
//    while (node != start) {
//        node = getParent(node);
//        ++pathLength;
//    }
//    treeNode** path = (treeNode**) vrna_alloc(pathLength * sizeof(treeNode*));
//    node = end;
//    for (int i = pathLength - 1; i >= 0; i--) {
//        path[i] = node;
//        node = getParent(node);
//    }
//    int* childIndex = (int*) vrna_alloc((pathLength-1) * sizeof(int));
//    for (int i = 0; i < pathLength - 1; i++) {
//        childIndex[i] = getChildIndex(path[i], getNodeID(path[i+1]));
//    }

//    fprintf(fp, "%%[PATH]");
//    for (int i = 0; i < pathLength; i++) {
//        fprintf(fp, " %c(%d/%d)", getNodeName(path[i]), (i < pathLength - 1 ? childIndex[i] + 1 : -1), path[i]->childCount);
//    }

    // ------- Boxes -------

    treeNode* exterior = getExterior(start);
    fprintf(fp,
            "/boxes [\n");
    PS_printFancyPathRec(fp, exterior, start, end, feature, 0);
    fprintf(fp, ""
                "] def \n");


    // ------- Tail -------
    tailColoredBoxes(fp);

    fclose(fp);
}

void PS_printFancySiblings(
        treeNode *node,
        treeNode *left,
        treeNode *right,
        puzzlerOptions* puzzler
) {
    (puzzler->psNumber)++;

    char filename[1000];
    if (puzzler->filename != NULL) {
        snprintf(filename, sizeof(filename), "%s_%06d_FANCY_TREE_%05d_%04d.ps", puzzler->filename, puzzler->psNumber, puzzler->numberOfChangesAppliedToConfig, getNodeID(node));
    } else {
        snprintf(filename, sizeof(filename), "FancyTree_%05d.ps", puzzler->numberOfChangesAppliedToConfig);
    }

    FILE* fp;
    fp = fopen(filename, "w");

    head(fp, puzzler->numberOfChangesAppliedToConfig);

    // ------- Boxes -------

    treeNode* exterior = getExterior(node);
    fprintf(fp,
            "/boxes [\n");
    PS_printFancySiblingsRec(fp, exterior, node, left, right, 0);
    fprintf(fp, ""
                "] def \n");


    // ------- Tail -------
    tailColoredBoxes(fp);

    fclose(fp);
}

void PS_printFancyTree(
        treeNode* node,
        puzzlerOptions* puzzler
) {
    (puzzler->psNumber)++;

    char filename[1000];
    if (puzzler->filename != NULL) {
        snprintf(filename, sizeof(filename), "%s_%06d_FANCY_TREE_%05d_%04d.ps", puzzler->filename, puzzler->psNumber, puzzler->numberOfChangesAppliedToConfig, getNodeID(node));
    } else {
        snprintf(filename, sizeof(filename), "FancyTree_%05d.ps", puzzler->numberOfChangesAppliedToConfig);
    }

    FILE* fp;
    fp = fopen(filename, "w");

    head(fp, puzzler->numberOfChangesAppliedToConfig);

    // ------- Boxes -------

    treeNode* exterior = getExterior(node);
    fprintf(fp,
            "/boxes [\n");
    PS_printFancyTreeRec(fp, exterior, node, 0);
    fprintf(fp, ""
                "] def \n");


    // ------- Tail -------
    tailColoredBoxes(fp);

    fclose(fp);
}
