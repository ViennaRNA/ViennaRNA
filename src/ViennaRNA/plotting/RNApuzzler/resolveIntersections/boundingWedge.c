#include "ViennaRNA/RNApuzzler/resolveIntersections/boundingWedge.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

void getBoundingWedgeRec(
    const treeNode* root,
    const treeNode* node,
    const double parentAngle,
    double* minAngle,
    double* maxAngle
) {

    /// --- Documentation ---
    ///
    /// How to ... get the bounding wedge of root's i-th child tree?
    ///
    /// get interesting points of current node
    ///     (2) touch points of tangents from root.center to node's loop circle
    ///     (n) bulge points of node's stem
    ///     (2) corners of node's stem that coincide with root's loop [only for direct child]
    ///
    /// update min-angle...
    ///     for each interesting point
    ///         check if the point is on the left side of the min-angle axis
    ///         if so
    ///             min-angle -= diff-angle of min-axis and point-axis
    /// update max-angle...
    ///     for each interesting point
    ///         check if the point is on the right side of the min-angle axis
    ///         if so
    ///             max-angle += diff-angle of max-axis and point-axis

    const double distance = epsilonFix;
    treeNode* parent = getParent(node);

    double centerRoot[2];
    getLoopCenter(root, centerRoot);
    double centerNode[2];
    getLoopCenter(node, centerNode);
    double vRootNode[2];
    vector(centerRoot, centerNode, vRootNode);

    /// set appropriate nodeAngle
    /// this could have been done using getChildAngle function
    /// but in terms of performance we can get this for free O(1)
    /// costs of getChildAngle: O( maxDegreeOnPath * ( depth(node) - depth(root) ) ) per call
    double nodeAngle;
    if (parent == root) {
        /// this happens only for the initial call and not for the recursive calls
        /// initialize min/max with the direct child's angle
        nodeAngle = getChildAngle(root, node);
        *minAngle = nodeAngle;
        *maxAngle = nodeAngle;
    } else {
        /// compare getChildAngle function
        double centerParent[2];
        getLoopCenter(parent, centerParent);
        double vRootParent[2];
        vector(centerRoot, centerParent, vRootParent);
        double diffParent = angleBetweenVectors2D(vRootParent, vRootNode);
        if (!isToTheRightPointVector(centerRoot, vRootParent, centerNode)) {
            diffParent *= -1;
        }
        nodeAngle = parentAngle + diffParent;
    }

    /// get all bounding boxes
    loopBox* loopNode = node->lBox;
    stemBox* stemNode = node->sBox;

    /// allocate space for points of interest
    int numPoints = stemNode->bulgeCount;
    if (parent == root) { numPoints += 2; } // for bottom corners of direct child's stem
    double** points = (double**) vrna_alloc(numPoints * sizeof(double*));
    int pointIndex = 0;

    /// points of interest (part 1)
    /// bulge points
    for (int i = 0; i < stemNode->bulgeCount; i++) {
        double o[2], q[2]; // o, q are unused but necessary for function call
        double* bulgePoint = (double*) vrna_alloc(2 * sizeof(double));
        getBulgeCoordinatesExtraDistance(stemNode, i, distance, o, bulgePoint, q);

        points[pointIndex] = bulgePoint;
        pointIndex++;
    }

    /// points of interest (part 2, only for direct child of root)
    /// for the direct child of this computation's root
    /// the corners of that child that coincide with root's loop (-> stem bottom corners)
    /// have a huge effect on the size of the bounding wedge
    /// for all further descendants these points do not matter because
    /// of the greater impact of the loops
    if (parent == root) {
        double* pStemBottomCornerL = (double*) vrna_alloc(2 * sizeof(double));
        pStemBottomCornerL[0] = stemNode->c[0] - stemNode->e[0] * stemNode->a[0] + stemNode->e[1] * stemNode->b[0];
        pStemBottomCornerL[1] = stemNode->c[1] - stemNode->e[0] * stemNode->a[1] + stemNode->e[1] * stemNode->b[1];
        points[pointIndex] = pStemBottomCornerL;
        pointIndex++;

        double* pStemBottomCornerR = (double*) vrna_alloc(2 * sizeof(double));
        pStemBottomCornerR[0] = stemNode->c[0] - stemNode->e[0] * stemNode->a[0] - stemNode->e[1] * stemNode->b[0];
        pStemBottomCornerR[1] = stemNode->c[1] - stemNode->e[0] * stemNode->a[1] - stemNode->e[1] * stemNode->b[1];
        points[pointIndex] = pStemBottomCornerR;
        pointIndex++;
    }

    /// we compute the two tangents from root's center to node's circle periphery
    /// using these for our min/max calculation ensures that the whole loop
    /// is contained in the wedge
    ///
    /// we can directly compute the diffAngle for the touching points of the tangents
    /// by using pythagoras' sentence
    double radiusNode = loopNode->r + distance;
    double distanceRootNode = vectorLength2D(vRootNode);
    /// positive angle and negative angle share their size
    double angle1 = asin(radiusNode / distanceRootNode);
    /// no need for double computation, just flip sign
    double angle2 = (-1) * angle1;

    // store both angles in an array to avoid code duplication
    double angles[2] = { angle1, angle2 };
    /// update min / max for the tangents touching points (angles)
    for (int currentAngle = 0; currentAngle < 2; currentAngle++) {
        double diffAngle = angles[currentAngle];
        double pointAngle = nodeAngle + diffAngle;

        /// actual updating for min and max
        if (pointAngle < *minAngle) {
            *minAngle = pointAngle;
        }
        if (pointAngle > *maxAngle) {
            *maxAngle = pointAngle;
        }
    }

    /// update min / max for bulge points (and stem bottom points in first level of recursion)
    for (int currentPoint = 0; currentPoint < numPoints; currentPoint++) {
        double* point = points[currentPoint];

        /// vector from root.loop.center to point
        double vCenterPoint[2];
        vector(centerRoot, point, vCenterPoint);

        /// (positive) offset angle between node.loop.center and point
        double diffAngle = angleBetweenVectors2D(vRootNode, vCenterPoint);
        double sign;
        if (isToTheRightPointVector(centerRoot, vRootNode, point)) {
            sign =  1;
        } else {
            sign = -1;
        }
        /// now the offset angle has some direction information
        diffAngle *= sign;

        /// the current point's angle has to be set with regard to the node's angle
        double pointAngle = nodeAngle + diffAngle;

        /// actual updating for min and max
        if (pointAngle < *minAngle) {
            *minAngle = pointAngle;
        }
        if (pointAngle > *maxAngle) {
            *maxAngle = pointAngle;
        }
    }

    /// free allocated space
    for (int currentPoint = 0; currentPoint < numPoints; currentPoint++) {
        double* point = points[currentPoint];
        free(point);
    }
    free(points);

    /// recursive call
    for (int currentChild = 0; currentChild < node->childCount; currentChild++) {
        treeNode* child = getChild(node, currentChild);
        getBoundingWedgeRec(root, child, nodeAngle, minAngle, maxAngle);
    }
}

void getBoundingWedge(
    const treeNode* root,
    const int childIndex,
    double* minAngle,
    double* maxAngle
) {
    treeNode* child = getChild(root, childIndex);

    /// not needed, done inside getBoundingWedgeRec at first recursion level
    /// we keep this stuff here for maintainance reasons
    //double childAngle;
    //getChildAngle(root, child, &childAngle);
    //*minAngle = *maxAngle = childAngle;

    getBoundingWedgeRec(root, child, 0, minAngle, maxAngle);
}

