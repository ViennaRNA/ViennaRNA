#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/output/output.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

void updateAABB(
        AABB *aabb,
        stemBox *sBox,
        loopBox *lBox
) {
    const double stem_ea[2] = { sBox->e[0] * sBox->a[0],
                                sBox->e[0] * sBox->a[1] };
    const double stem_eb[2] = { sBox->e[1] * sBox->b[0],
                                sBox->e[1] * sBox->b[1] };

    int numPoints = 6 + sBox->bulgeCount;

    /// array of relevant points
    double** p = (double**) vrna_alloc(numPoints * sizeof(double*));
    for (int i = 0; i < numPoints; i++) {
        p[i] = (double*) vrna_alloc(2 * sizeof(double));
    }

    /// corners of stem
    p[0][0] = sBox->c[0] - stem_ea[0] + stem_eb[0];
    p[0][1] = sBox->c[1] - stem_ea[1] + stem_eb[1];
    p[1][0] = sBox->c[0] + stem_ea[0] + stem_eb[0];
    p[1][1] = sBox->c[1] + stem_ea[1] + stem_eb[1];
    p[2][0] = sBox->c[0] + stem_ea[0] - stem_eb[0];
    p[2][1] = sBox->c[1] + stem_ea[1] - stem_eb[1];
    p[3][0] = sBox->c[0] - stem_ea[0] - stem_eb[0];
    p[3][1] = sBox->c[1] - stem_ea[1] - stem_eb[1];

    /// lower left of loop AABB
    p[4][0] = lBox->c[0] - lBox->r;
    p[4][1] = lBox->c[1] - lBox->r;
    /// upper right of loop AABB
    p[5][0] = lBox->c[0] + lBox->r;
    p[5][1] = lBox->c[1] + lBox->r;

    /// bulge points
    double pPrev[2], pNext[2];
    for (int i = 0; i < sBox->bulgeCount; i++) {
        getBulgeCoordinates(sBox, i, pPrev, p[6+i],pNext);
    }

    /// set aabb
    aabb->min[0] = p[0][0];
    aabb->min[1] = p[0][1];
    aabb->max[0] = p[0][0];
    aabb->max[1] = p[0][1];
    for (int i = 1; i < numPoints; i++) {
        if (aabb->min[0] > p[i][0]) {
            aabb->min[0] = p[i][0];
        }
        if (aabb->min[1] > p[i][1]) {
            aabb->min[1] = p[i][1];
        }
        if (aabb->max[0] < p[i][0]) {
            aabb->max[0] = p[i][0];
        }
        if (aabb->max[1] < p[i][1]) {
            aabb->max[1] = p[i][1];
        }
    }

    /// free allocated memory
    for (int i = 0; i < numPoints; i++) {
        free(p[i]);
    }
    free(p);
}

void updateBoundingBoxes(
        treeNode* node,
        const puzzlerOptions* puzzler
) {
    /// fix this node's loop
    /// then for each child fix the stem and bulges
    /// and call recursively

    if (!isExterior(node)) {
        long numStemBackBones = lround((2.0 * node->sBox->e[0]) / puzzler->unpaired);
        double stemLength = puzzler->unpaired * numStemBackBones;
        double distanceStemEndToLoopCenter = sqrt(  node->cfg->radius *   node->cfg->radius - 0.25 * puzzler->paired * puzzler->paired);
        double distanceStemCenterToLoopCenter = 0.5 * stemLength + distanceStemEndToLoopCenter;
        node->lBox->c[0] = node->sBox->c[0] + distanceStemCenterToLoopCenter * node->sBox->a[0];
        node->lBox->c[1] = node->sBox->c[1] + distanceStemCenterToLoopCenter * node->sBox->a[1];
        node->lBox->r = node->cfg->radius;

        updateAABB(&(node->aabb), node->sBox, node->lBox);
    }

    double childAngleRad = 0.0;

    for (int i = 0; i < node->childCount; i++) {

        treeNode* child = getChild(node, i);
        stemBox* sBox = child->sBox;
        loopBox* lBox = child->lBox;

        double parentLoopCenter[2];
        if (isExterior(node)) {
            parentLoopCenter[0] = lBox->c[0];
            parentLoopCenter[1] = EXTERIOR_Y;
        } else {
            getLoopCenter(node, parentLoopCenter);
        }

        /// ... fix the stem's extensions ...
        long numStemBackBones = lround((2.0 * sBox->e[0]) / puzzler->unpaired);
        double stemLength = puzzler->unpaired * numStemBackBones;
        sBox->e[0] = 0.5 * stemLength;
        sBox->e[1] = 0.5 * puzzler->paired;

        /// ... fix the stem's directions ...
        if (isExterior(node)) {
            childAngleRad = MATH_PI;
        } else {
            childAngleRad += getArcAngle(node->cfg, i);
        }
        double aFixed[2];
        if (isExterior(node)) {
            aFixed[0] = 0.0;
            aFixed[1] = 1.0;
        } else {
            double gamma = childAngleRad - MATH_PI;
            rotateVectorByAngle(node->sBox->a, gamma, aFixed);
        }
        sBox->a[0] = aFixed[0];
        sBox->a[1] = aFixed[1];

        double bFixed[2];
        normal(aFixed, bFixed);
        bFixed[0] *= -1;
        bFixed[1] *= -1;
        sBox->b[0] = bFixed[0];
        sBox->b[1] = bFixed[1];

        double s0 = 0;
        if (!isExterior(node)) {
            s0 = sqrt(node->cfg->radius * node->cfg->radius - 0.25 * puzzler->paired * puzzler->paired);
        }
        double distanceStemCenter = s0 + 0.5 * stemLength;

        /// ... fix the stem's position.
        sBox->c[0] = parentLoopCenter[0] + distanceStemCenter * aFixed[0];
        sBox->c[1] = parentLoopCenter[1] + distanceStemCenter * aFixed[1];

        if (stemLength == 0) {
            sBox->e[0] = epsilon7;
        }
    }

    for (int i = 0; i < node->childCount; i++) {
        updateBoundingBoxes(getChild(node, i), puzzler);
    }
}

void applyChangesToConfigAndBoundingBoxes(
    treeNode*       tree,
    const double*   deltaCfg,
    const double    radiusNew,
    const puzzlerOptions* puzzler
) {
    /// Apply all changes to config
    config *cfg = tree->cfg;

    /// - start with adjusting config radius and angles
    cfgApplyChanges(cfg, getNodeName(tree), deltaCfg, radiusNew, puzzler);

    /// - apply changes of config to bounding boxes
    updateBoundingBoxes(tree, puzzler);

    //if (GEOGEBRA_FLAG) {
    //    GEOGEBRA_generateTree(global_root2, puzzler->numberOfChangesAppliedToConfig);
    //}
}

void freeBulges(stemBox* sBox) {
    if (sBox->bulges != NULL) {
        for (int currentBulge = 0; currentBulge < sBox->bulgeCount; currentBulge++) {
            free(sBox->bulges[currentBulge]);
        }
        free(sBox->bulges);
    }
}

void freeTree(treeNode* node) {
    for (int currentChild = 0; currentChild < node->childCount; currentChild++) {
        freeTree(getChild(node, currentChild));
    }

    if (node->cfg) {
        cfgFreeConfig(node->cfg);
    }
    if (node->children) {
        free(node->children);
    }
    if (node->lBox) {
        free(node->lBox);
    }
    if (node->sBox) {
        freeBulges(node->sBox);
        free(node->sBox);
    }
    free(node);
}

int countSubtreeNodes(const treeNode* node) {
    int count = 1; // count this node

    for (int currentChild = 0; currentChild < node->childCount; currentChild++) {
        // count children and add child count
        count += countSubtreeNodes(getChild(node, currentChild));
    }

    return count;
}

int countAncestorNodes(const treeNode* node) {
    int count = 0;

    treeNode *ancestor = getParent(node);
    while (ancestor != NULL) {
        ++count;
        ancestor = getParent(ancestor);
    }

    return count;
}

int collectSubtreeNodes(
        treeNode* node,
        treeNode** const allNodes,
        const int currentIndex
) {
    allNodes[currentIndex] = node;
    int nextIndex = currentIndex + 1; // increase index as this one was just taken

    for (int currentChild = 0; currentChild < node->childCount; currentChild++) {
        nextIndex = collectSubtreeNodes(getChild(node, currentChild), allNodes, nextIndex);
    }

    return nextIndex;
}

void collectAncestorNodes(
        const treeNode* node,
        treeNode** const ancestorList
) {
    int currentIndex = 0;
    treeNode *ancestor = getParent(node);
    while (ancestor != NULL) {
        ancestorList[currentIndex] = ancestor;
        ++currentIndex;
        ancestor = getParent(ancestor);
    }
}

double getPairedAngle(const treeNode* node) {
    /// get the current node's stem's bounding wedge
    stemBox* sBox = node->sBox;
    double pStemTopCorner[2];
    pStemTopCorner[0] = sBox->c[0] + sBox->e[0] * sBox->a[0] + sBox->e[1] * sBox->b[0];
    pStemTopCorner[1] = sBox->c[1] + sBox->e[0] * sBox->a[1] + sBox->e[1] * sBox->b[1];
    double pLoopCenter[2];
    getLoopCenter(node, pLoopCenter);
    double vLoopCenterToStemTopCorner[2];
    vector(pLoopCenter, pStemTopCorner, vLoopCenterToStemTopCorner);
    double vLoopCenterToStemCenter[2] = { (-1) * sBox->a[0], (-1) * sBox->a[1] };
    double minOuterAngle = angleBetweenVectors2D(vLoopCenterToStemCenter, vLoopCenterToStemTopCorner);

    /// all arcs share the same stem angle
    double stemAngle = 2 * minOuterAngle;
    return stemAngle;
}

short isExterior(const treeNode* node) {
    return getNodeID(node) == 0;
}

treeNode* getExterior(treeNode* node) {
    treeNode* exteriorCandidate = node;
    while (!isExterior(exteriorCandidate)) {
        exteriorCandidate = getParent(exteriorCandidate);
    }
    return exteriorCandidate;
}

short isInteriorLoop(const treeNode* node) {
    return
      (!isExterior(node)
       && node->childCount == 1);
}

short isMultiLoop(const treeNode* node) {
    return
      (!isExterior(node)
       && node->childCount > 1);
}

int getNodeID(const treeNode* node) {
    if (node != NULL) {
        return node->id;
    } else {
        return -1;
    }
}

char getNodeName(const treeNode* node) {
    /**
     * @brief cfgMotivBlank
     *      - name of exterior loops and small bulge loops
     *        initialized at cfgGenerateMotivs
     */
    const char cfgMotivBlank = '_';

    int id = getNodeID(node);
    if (id == -1) {
        return cfgMotivBlank;
    }

    int motivId = (id + 33) % 128;
    while (motivId < 33) {
        motivId = (motivId + 33) % 128;
    }
    char motiv = (char) motivId;
    //printf("[CONVERT] %3d -> %3d -> %c\n", id, motivId, motiv);
    return motiv;
}

void setChild(
        treeNode* parent,
        const int index,
        treeNode* child
) {
    if (0 <= index && index < parent->childCount) {
        parent->children[index] = child;
    }
}

/**
 * @brief treeGetChildCount
 *      - counts the number of children this loop will have in the configtree
 * @param loopStart
 *      - index of the loops first base
 * @param pair_table
 *      - the RNA's pairing information
 * @return
 *      - number of child nodes
 */
int treeGetChildCount(const int loopStart, const short* const pair_table) {
    int childCount = 0;

    int end = pair_table[loopStart];

    for (int i = loopStart + 1; i < end; ++i ) {
        if (pair_table[i] > i) {
            // found new stem
            childCount++;
            i = pair_table[i];
        }
    }
    return childCount;
}

/**
 * @brief createTreeNode
 *      - this method can be referred to as a constructor method for configtree nodes.
 * @param parent
 *      - parent node (the prior loop), NULL for the root node
 * @param loopStart
 *      - index of the loops first node, 1 for root node
 * @param stemStart
 *      - index of the prior stems first node, -1 for root node
 * @param pair_table
 *      - the RNA's pairing information
 * @param cfg
 *      - the configuration found in baseInformation for that loop.
 *        NULL for exterior loop (root node)
 * @return
 *      - an initialized configtree tBaseInformation with set parent, loopStart, cfg, childCount and initialized children array
 */
treeNode* createTreeNode(
    const int id,
    treeNode* parent,
    const int loopStart,
    const int stemStart,
    const short* const pair_table,
    config* cfg
) {
    // allocate children array
    int childCount;
    if (cfg == NULL) {
        childCount = treeGetChildCount(0, pair_table);
    } else {
        childCount = treeGetChildCount(loopStart, pair_table);
    }
    treeNode** children = (childCount > 0) ? (treeNode**) vrna_alloc(childCount*sizeof(treeNode*))
                                           : NULL;

    treeNode* node = (treeNode*) vrna_alloc(1 * sizeof(treeNode));

    node->id = id;

    node->parent = parent;
    node->children = children;
    node->childCount = childCount;

    node->cfg = cfg;
    node->loop_start = loopStart;
    node->stem_start = stemStart;

    node->lBox = NULL;
    node->sBox = NULL;

    return node;
}

// prototype for stem - loop recursion
treeNode* treeHandleStem(
    treeNode* parent,
    int *nodeID,
    const int stemStart,
    const short* const pair_table,
    const tBaseInformation* baseInformation
);

/**
 * @brief treeHandleLoop
 *      - method for configtree construction.
 *        uses recursive calls alternating with treeHandleStem method to get the whole RNA
 * @param parent
 *      - parent node of the current loop in configtree
 * @param loopStart
 *      - index of the loop's first base
 * @param stemStart
 *      - index of the prior stem's first base
 * @param pairTable
 *      - the RNA's pairing information
 * @param baseInformation
 *      - array of tBaseInformation annotations (grants config)
 * @return
 *      - pointer to the subtree (configtree) that has this loop as root
 */
treeNode* treeHandleLoop(
    treeNode* parent,
    int *nodeID,
    const int loopStart,
    const int stemStart,
    const short* const pair_table,
    const tBaseInformation* baseInformation
) {
    int addedChildren = 0;

    treeNode* subtree = createTreeNode(
        *nodeID,
        parent,
        loopStart,
        stemStart,
        pair_table,
        baseInformation[loopStart].config);

    int end = pair_table[loopStart];

    for (int i = loopStart + 1; i < end; ++i) {
        if (pair_table[i] > i) {
            // found new stem
            treeNode* child = treeHandleStem(subtree, nodeID, i, pair_table, baseInformation);
            child->parent = subtree;
            setChild(subtree, addedChildren, child);
            addedChildren++;

            i = pair_table[i];
        }
    }

    return subtree;
}

/**
 * @brief treeHandleStem
 *      - method for configtree construction.
 *        uses recursive calls alternating with treeHandleLoop method to get the whole RNA
 * @param parent
 *      - parent node of the current loop in configtree
 * @param stemStart
 *      - index of the stem's first base
 * @param pair_table
 *      - the RNA's pairing information
 * @param baseInformation
 *      - array of tBaseInformation annotations (grants config)
 * @return
 *      - pointer to the cunsecutive subtree (configtree) that is created from the consecutive loop
 */
treeNode* treeHandleStem(
    treeNode* parent,
    int *nodeID,
    const int stemStart,
    const short* const pair_table,
    const tBaseInformation* baseInformation
) {
    ++(*nodeID);
    //printf("New stem: %c\n", *nodeID);
    int i = stemStart;
    while (baseInformation[i].config == NULL) {
        ++i;
    }

    treeNode* child = treeHandleLoop(parent, nodeID, i, stemStart, pair_table, baseInformation);
    return child;
}

/**
 * @brief buildBoundingBoxes
 * @param tree
 * @param pair_table
 * @param baseInformation
 * @param x
 * @param y
 * @param bulge
 *      - distance between regular stem and bulge base
 */
void buildBoundingBoxes(
    treeNode* tree,
    const short* const pair_table,
    const tBaseInformation* baseInformation,
    const double* x,
    const double* y,
    const double bulge
) {
    short isRoot = (tree->parent == NULL);

    if (!isRoot) {
        loopBox* lBox = buildLoopBox(tree->loop_start, pair_table, baseInformation, x, y);
        stemBox* sBox = buildStemBox(tree->stem_start, tree->loop_start, pair_table, x, y, bulge);

        lBox->parent = tree;
        sBox->parent = tree;

        tree->lBox = lBox;
        tree->sBox = sBox;

        updateAABB(&(tree->aabb), sBox, lBox);
    }

    for (int currentChild = 0; currentChild < tree->childCount; currentChild++) {
        treeNode* child = getChild(tree, currentChild);
        buildBoundingBoxes(child, pair_table, baseInformation, x, y, bulge);
    }
}

// documentation at header file
treeNode* buildConfigtree(
        const short* const pair_table,
        const tBaseInformation* baseInformation,
        const double* x,
        const double* y,
        const double bulge
) {
    //printf("buildConfigtree: started\n");

    // create root
    int nodeID = 0;
    treeNode* root = createTreeNode(nodeID, NULL, 1, -1, pair_table, NULL);

    int addedChildren = 0;
    int length = pair_table[0];
    for (int i = 1; i < length; ++i) {
        if (pair_table[i] > i) {
            // found stem
            treeNode* child = treeHandleStem(root, &nodeID, i, pair_table, baseInformation);
            setChild(root, addedChildren, child);
            addedChildren++;

            i = pair_table[i];
        }
    }

    buildBoundingBoxes(root, pair_table, baseInformation, x, y, bulge);

    return root;
}

/**
 * @brief translateBoundingBoxesByVector
 *      - Performs a translation of a whole branch by a given vector.
 *        Used to apply changes in config to the tree and its bounding boxes.
 * @param tree
 *      - tree that is being translated
 * @param vector
 *      - translation vector as double array[2]
 */
void translateBoundingBoxes(
        treeNode* tree,
        const double* vector
) {
    translateStemBox(tree->sBox, vector);
    translateLoopBox(tree->lBox, vector);
    updateAABB(&(tree->aabb), tree->sBox, tree->lBox);
    for (int currentChild = 0; currentChild < tree->childCount; currentChild++) {
        translateBoundingBoxes(getChild(tree, currentChild), vector);
    }
}

/**
 * @brief getChildIndex
 *      - gets the index of child node where to find the node with given name.
 * @param tree
 *      - configtree you want to search in.
 * @param childID
 *      - ID of childnode you are looking for.
 * @return
 *      - child index or -1 if tree does not contain such a childnode.
 */
int getChildIndex(
        const treeNode* tree,
        const int childID
) {
    // check if there are further nodes to check
    int childIndex = tree->childCount - 1;
    for (int currentChild = 0;
         currentChild < tree->childCount;
         ++currentChild) {
        treeNode* child = getChild(tree, currentChild);
        if (getNodeID(child) > childID) {
            childIndex = currentChild - 1;
            break;
        }
    }

    return childIndex;
}

/**
 * @brief getChildAngle
 *      - Calculates the clockwise angle of a given child node (given by name) at a loop.
 *        This child needs to be a direct child node
 *        The rotation angle of its center node will be calculated.
 * @param parentNode
 *      - tree node acting as parent loop
 * @param childNode
 *      - child node of which you want to know the angle
 * @return
 *      - angle of child node
 *        the resulting angle might be smaller than 0° or greater than 360°
 */
double getChildAngle(
    const treeNode* parentNode,
    const treeNode* childNode
) {
    double parentLoopCenter[2] = { parentNode->lBox->c[0], parentNode->lBox->c[1] };
    double parentStemCenter[2] = { parentNode->sBox->c[0], parentNode->sBox->c[1] };
    double parentLoopStemVector[2];
    vector(parentLoopCenter, parentStemCenter, parentLoopStemVector);

    double childLoopCenter[2] = { childNode->lBox->c[0], childNode->lBox->c[1] };
    double angle = anglePtPtPt2D(parentStemCenter, parentLoopCenter, childLoopCenter);
    if (!isToTheRightPointVector(parentLoopCenter, parentLoopStemVector, childLoopCenter)) {
        angle = MATH_TWO_PI - angle;
    }

    return angle;
}

/**
 * @brief getChildAngleByIndex
 *      - Calculates the clockwise angle of a given child node (given by name) at a loop.
 *        This child needs to be a direct child node
 *        The rotation angle of its center node will be calculated.
 * @param parentNode
 *      - tree node acting as parent loop
 * @param childIndex
 *      - index of child node of which you want to know the angle
 * @return
 *      - angle of child node
 *        the resulting angle might be smaller than 0° or greater than 360°
 */
double getChildAngleByIndex(
    const treeNode* parentNode,
    const int childIndex
) {
    return getChildAngle(parentNode, getChild(parentNode, childIndex));
}

/**
 * @brief getLoopCenter
 *      - Getter for the center coordinates of a tree node's loop.
 * @param node
 *      - your tree node
 * @param p
 *      - double[2] return value for the loop's center coordinates
 */
void getLoopCenter(
        const treeNode* node,
        double p[2]
) {
    getLBoxCenter(node->lBox, p);
}

/**
 * @brief getStemCenter
 *      - Getter for the center coordinates of a tree node's stem.
 * @param node
 *      - your tree node
 * @param p
 *      - double[2] return value for the stem's center coordinates
 */
void getStemCenter(
        const treeNode* node,
        double p[2]
) {
    getSBoxCenter(node->sBox, p);
}

/**
 * @brief getChildNode
 *      - searches for the node with given childname in a configtree.
 * @param tree
 *      - configtree you want to search in.
 * @param childID
 *      - ID of childnode you are looking for.
 * @return
 *      - ptr to node or NULL if tree does not contain such a childnode.
 */
treeNode* getChildNode(
    treeNode* tree,
    const int childID
) {
    short treeIsRoot = isExterior(tree);

    // check if this is the wanted node
    if (!treeIsRoot) {
        if (getNodeID(tree) == childID) {
            return tree;
        }
    }

    // get index of next child on our path
    treeNode* child = getChild(tree, getChildIndex(tree, childID));
    if (child == NULL) {
        return NULL;
    } else {
        return getChildNode(child, childID);
    }
}

treeNode* getChild(
        const treeNode* node,
        const int index
) {
    if (node == NULL) {
      return NULL;
    } else if (index < 0) {
      return NULL;
    } else if (index >= node->childCount) {
      return NULL;
    } else {
      return node->children[index];
    }
}

treeNode* getParent(const treeNode* node) {
    if (node == NULL) {
      return NULL;
    } else {
      return node->parent;
    }
}

void printTree(const treeNode* node, const int level) {
    for (int i = 0; i < level; i++) {
      printf("#");
    }

    printf(" %c(%d)\n", getNodeName(node), getNodeID(node));

    for (int i = 0; i < node->childCount; i++) {
        printTree(getChild(node, i), level + 1);
    }
}

