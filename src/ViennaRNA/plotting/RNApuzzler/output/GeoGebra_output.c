#include "ViennaRNA/RNApuzzler/output/GeoGebra_output.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/optimize.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/resolveIntersections.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"

#include <stdio.h>
#include <stdlib.h>

void GEOGEBRA_printStem(stemBox* rect) {
    char outName[10];
    snprintf(outName, sizeof(outName), "stem_{%c}", getNodeName(rect->parent));
    GEOGEBRA_printRectangle(outName, rect->a, rect->b, rect->c, rect->e);
}

void GEOGEBRA_printLoop(loopBox* circ) {
    char outName[10];
    snprintf(outName, sizeof(outName), "loop_{%c}", getNodeName(circ->parent));
    printf("[GEOGEBRA] %s : %f = (x - %f)² + (y - %f)²\n"
           , outName, (circ->r * circ->r), circ->c[0], circ->c[1]);
}

void GEOGEBRA_printNode(treeNode* node, short printStem, short printLoop) {
    char name = getNodeName(node);
    printf("[GEOGEBRA] Node %c...\n", name);

    if (printStem) {
        GEOGEBRA_printStem(node->sBox);
    }

    if (printLoop) {
        GEOGEBRA_printLoop(node->lBox);
    }
}

void GEOGEBRA_printPath(treeNode* nodeTOP, treeNode* nodeBOTTOM) {
    treeNode* node = nodeTOP;
    while (node != nodeBOTTOM) {
        GEOGEBRA_printNode(node, 1, 1);
        int indexNext = getChildIndex(node, getNodeID(nodeBOTTOM));
        node = getChild(node, indexNext);
    }
    GEOGEBRA_printNode(node, 1, 1);
}

void GEOGEBRA_printTree(treeNode* node) {
    GEOGEBRA_printNode(node, 1, 1);
    int i;
    for (i = 0; i < node->childCount; i++) {
        GEOGEBRA_printTree(getChild(node, i));
    }
}

void GEOGEBRA_generateHead(FILE* fp) {

    fprintf(fp,
    "<?xml version='1.0' encoding='utf-8'?>\n"
    "<geogebra format='4.0' id='cd099ee0-1acf-462b-8863-25eb5684fe6a'  xsi:noNamespaceSchemaLocation='http://www.geogebra.org/ggb.xsd' xmlns='' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' >\n"
    "	<gui>\n"
    "		<window width='1215' height='776' />\n"
    "		<perspectives>\n"
    "			<perspective id='tmp'>\n"
    "				<panes>\n"
    "					<pane location='' divider='0.3432098765432099' orientation='1' />\n"
    "				</panes>\n"
    "				<views>\n"
//    "					<view id='4' toolbar='0 || 2020 , 2021 , 2022 , 66 || 2001 , 2003 , 2002 , 2004 , 2005 || 2040 , 2041 , 2042 , 2044 , 2043' visible='false' inframe='false' stylebar='false' location='1,1' size='300' window='100,100,600,400' />\n"
    "					<view id='1' visible='true' inframe='false' stylebar='false' location='1' size='788' window='100,100,600,400' />\n"
//    "					<view id='2' visible='true' inframe='false' stylebar='false' location='3' size='417' window='100,100,250,400' />\n"
//    "					<view id='16' visible='false' inframe='true' stylebar='false' location='1' size='150' window='50,50,500,500' />\n"
//    "					<view id='32' visible='false' inframe='true' stylebar='true' location='1' size='150' window='50,50,500,500' />\n"
//    "					<view id='64' visible='false' inframe='true' stylebar='false' location='1' size='150' window='50,50,500,500' />\n"
    "				</views>\n"
//    "				<toolbar show='true' items='0 39 59 || 1 501 67 , 5 19 , 72 | 2 15 45 , 18 65 , 7 37 | 4 3 8 9 , 13 44 , 58 , 47 || 16 51 64 , 70 | 10 34 53 11 , 24  20 22 , 21 23 | 55 56 57 , 12 || 36 46 , 38 49 50 , 71 | 30 29 54 32 31 33 | 17 26 62 , 14 66 68 | 25 52 60 61 || 40 41 42 , 27 28 35 , 6' />\n"
//    "				<input show='true' cmd='true' top='false' />\n"
    "			</perspective>\n"
    "		</perspectives>\n"
    "		<labelingStyle  val='3'/>\n"
    "		<font  size='12'/>\n"
    "		<graphicsSettings javaLatexFonts='false'/>\n"
    "		<consProtColumns  col0='true' col1='true' col2='false' col3='true' col4='false' col5='true' col6='true' col7='false'/>\n"
    "		<consProtocol useColors='true' addIcons='false' showOnlyBreakpoints='false'/>\n"
    "	</gui>\n"
    "	<kernel>\n"
    "		<continuous val='false'/>\n"
    "		<decimals val='2'/>\n"
    "		<angleUnit val='degree'/>\n"
    "		<algebraStyle val='0'/>\n"
    "		<coordStyle val='0'/>\n"
    "		<angleFromInvTrig val='false'/>\n"
    "	</kernel>\n"
    "	<scripting blocked='false' disabled='false'/>\n"
    );

}

void getBounds(treeNode* node, double* xMin, double* xMax, double* yMin, double* yMax) {
    loopBox* lBox = node->lBox;
    // compute current node's loopbox
    double xMaxLoop = lBox->c[0] + lBox->r;
    double xMinLoop = lBox->c[0] - lBox->r;
    double yMaxLoop = lBox->c[1] + lBox->r;
    double yMinLoop = lBox->c[1] - lBox->r;
    *xMax = (xMaxLoop > *xMax) ? xMaxLoop : *xMax;
    *xMin = (xMinLoop < *xMin) ? xMinLoop : *xMin;
    *yMax = (yMaxLoop > *yMax) ? yMaxLoop : *yMax;
    *yMin = (yMinLoop < *yMin) ? yMinLoop : *yMin;

    stemBox* sBox = node->sBox;
    // compute current node's stembox
    short signs[2] = { 1, -1 };
    int k, l;
    for (k = 0; k < 2; k++) {
        for (l = 0; l < 2; l++) {
            double p[2] = { sBox->c[0] + signs[k] * sBox->a[0] + signs[l] * sBox->b[0],
                            sBox->c[1] + signs[k] * sBox->a[1] + signs[l] * sBox->b[1] };
            *xMax = (p[0] > *xMax) ? p[0] : *xMax;
            *xMin = (p[0] < *xMin) ? p[0] : *xMin;
            *yMax = (p[1] > *yMax) ? p[1] : *yMax;
            *yMin = (p[1] < *yMin) ? p[1] : *yMin;
        }
    }
    // compute current node's bulges
    for (k = 0; k < sBox->bulgeCount; k++) {
        double o[2], p[2], q[2]; // o, q are unused but necessary for function call
        getBulgeCoordinates(sBox, k, o, p ,q);
        *xMax = (p[0] > *xMax) ? p[0] : *xMax;
        *xMin = (p[0] < *xMin) ? p[0] : *xMin;
        *yMax = (p[1] > *yMax) ? p[1] : *yMax;
        *yMin = (p[1] < *yMin) ? p[1] : *yMin;
    }

    // iterate to all child nodes
    int i;
    for (i = 0; i < node->childCount; i++) {
        getBounds(getChild(node, i), xMin, xMax, yMin, yMax);
    }
}

void GEOGEBRA_adjustView(FILE* fp_xml, const treeNode* node) {
    double xMin =  1000000.0;
    double xMax = -1000000.0;
    double yMin =  1000000.0;
    double yMax = -1000000.0;

    int i;
    for (i = 0; i < node->childCount; i++) {
        getBounds(getChild(node, i), &xMin, &xMax, &yMin, &yMax);
    }

    double height = 400;//651;//788;
    double width  = 400;//788;//651;
    double scale = 1;//0.6665023101914997;

    double w = xMax - xMin;
    double h = yMax - yMin;
    double scaleH = height / h;
    double scaleW = width / w;
    double yscale = 1;

    scale  = scaleH < scaleW ? scaleH : scaleW;
    yscale = scale;

//    scale = 0.5;
//    yscale = 0.5;

//    scale = scaleW;
//    yscale = scaleH;

    height = h;
    width  = w;

    double xZero = scale * -1*xMin;
    double yZero = yscale * yMax;//height;

    printf("[BOUNDS] min(%3.2f, %3.2f) max(%3.2f, %3.2f)\n", xMin, yMin, xMax, yMax);
    printf(
    "		<size  width='%f' height='%f'/>\n"
    "		<coordSystem xZero='%f' yZero='%f' scale='%f' yscale='%f'/>\n"
    , width, height
    , xZero, yZero, scale, yscale
    );
    /*
        <size  width="788" height="651"/>
        <coordSystem xZero="153.4800220446751" yZero="719.0633934846451" scale="0.6946884293059696" yscale="0.6946884293059706"/>
    */

    fprintf(fp_xml,
    "	<euclidianView>\n"
    "		<size  height='%f' width='%f'/>\n"
    "		<coordSystem xZero='%f' yZero='%f' scale='%f' yscale='%f'/>\n"
    "		<evSettings axes='true' grid='false' gridIsBold='false' pointCapturing='3' rightAngleStyle='1' checkboxSize='13' gridType='0'/>\n"
    "		<bgColor r='255' g='255' b='255'/>\n"
    "		<axesColor r='0' g='0' b='0'/>\n"
    "		<gridColor r='192' g='192' b='192'/>\n"
    "		<lineStyle axes='1' grid='10'/>\n"
    "		<axis id='0' show='false' label='' unitLabel='' tickStyle='1' showNumbers='true'/>\n"
    "		<axis id='1' show='false' label='' unitLabel='' tickStyle='1' showNumbers='true'/>\n"
//    "		<axis id='0' show='true' label='' unitLabel='' tickStyle='1' showNumbers='true'/>\n"
//    "		<axis id='1' show='true' label='' unitLabel='' tickStyle='1' showNumbers='true'/>\n"
    "	</euclidianView>\n"
    , width, height
    , xZero, yZero, scale, yscale
    );
}

void GEOGEBRA_generateHeadEnd(FILE* fp) {
    fprintf(fp,
    "	<construction title='' author='' date=''>\n"
    );
}

void GEOGEBRA_generateTail(FILE* fp) {
    fprintf(fp,
    "	</construction>\n"
    "</geogebra>\n"
    );
}

void GEOGEBRA_generateStem(FILE* fp, const treeNode* node) {
    int id = getNodeID(node);
    stemBox* s = node->sBox;
    double A0[2] = { s->c[0] + s->e[0] * s->a[0] + s->e[1] * s->b[0],
                     s->c[1] + s->e[0] * s->a[1] + s->e[1] * s->b[1] };
    double B0[2] = { s->c[0] + s->e[0] * s->a[0] - s->e[1] * s->b[0],
                     s->c[1] + s->e[0] * s->a[1] - s->e[1] * s->b[1] };
    double C0[2] = { s->c[0] - s->e[0] * s->a[0] - s->e[1] * s->b[0],
                     s->c[1] - s->e[0] * s->a[1] - s->e[1] * s->b[1] };
    double D0[2] = { s->c[0] - s->e[0] * s->a[0] + s->e[1] * s->b[0],
                     s->c[1] - s->e[0] * s->a[1] + s->e[1] * s->b[1] };
    fprintf(fp,
    "		<command name='Polygon'>\n"
    "			<input a0='(%f, %f)' a1='(%f, %f)' a2='(%f, %f)' a3='(%f, %f)'/>\n"
    "			<output a0='stem_{%d}' a1='l%d_1' a2='l%d_2' a3='l%d_3' a4='l%d_4'/>\n"
    "		</command>\n"
    , A0[0], A0[1], B0[0], B0[1], C0[0], C0[1], D0[0], D0[1]
    , id, id, id, id, id
    );
}

void GEOGEBRA_generateLoop(FILE* fp, const treeNode* node) {
    int id = getNodeID(node);
    loopBox* l = node->lBox;
    fprintf(fp,
    "           <expression label='loop_{%d}' exp='%f² = (x - (%f))² + (y - (%f))²'/>\n"
    , id, l->r, l->c[0], l->c[1]
    );
}

void GEOGEBRA_generateNode(FILE* fp, treeNode* node, short printStem, short printLoop) {
    if (node == NULL) {
        return;
    }
    if (isExterior(node)) {
        return;
    }
    if (printStem) {
        GEOGEBRA_generateStem(fp, node);
    }
    if (printLoop) {
        GEOGEBRA_generateLoop(fp, node);
    }
}

void GEOGEBRA_generatingTree(FILE* fp, treeNode* node) {
    GEOGEBRA_generateNode(fp, node, 1, 1);
    int i;
    for (i = 0; i < node->childCount; i++) {
        GEOGEBRA_generatingTree(fp, getChild(node, i));
    }
}

void GEOGEBRA_generatingPath(FILE* fp_xml, treeNode* nodeTOP, treeNode* nodeBOTTOM) {
    treeNode* node = nodeTOP;
    while (node != nodeBOTTOM) {
        GEOGEBRA_generateNode(fp_xml, node, 1, 1);
        int indexNext = getChildIndex(node, getNodeID(nodeBOTTOM));
        node = getChild(node, indexNext);
    }
    GEOGEBRA_generateNode(fp_xml, node, 1, 1);
}

void GEOGEBRA_generating(
        treeNode* node1,
        treeNode* node2,
        int changeNumber,
        intersectionType it,
        short step
) {
    FILE* fp_xml;
    fp_xml = fopen("geogebra.xml", "w");
    GEOGEBRA_generateHead(fp_xml);
    // get bounds and adjust view
    GEOGEBRA_adjustView(fp_xml, node1);

    GEOGEBRA_generateHeadEnd(fp_xml);
    if (node2 == NULL) {
        // generateTree
        GEOGEBRA_generatingTree(fp_xml, node1);
    } else {
        // generatePath
        GEOGEBRA_generatingPath(fp_xml, node1, node2);
    }

    GEOGEBRA_generateTail(fp_xml);
    fclose(fp_xml);

    FILE* fp_js;
    fp_js = fopen("geogebra_javascript.js", "w");
    fprintf(fp_js, "function ggbOnInit() {}\n");
    fclose(fp_js);

    char script[100];
    snprintf(script, sizeof(script),
            ""
            "/bin/bash \n"
            "zip -m GEO.zip geogebra.xml geogebra_javascript.js \n"
            "mv GEO.zip GEO_%s%s%s%d%s_%s.ggb \n"
            ""
             , (changeNumber < 1000 ? "0" : "")
             , (changeNumber <  100 ? "0" : "")
             , (changeNumber <   10 ? "0" : "")
             ,  changeNumber
             , (step == 0 ? "_0" :
               (step == 1 ? "_1" : ""))
             , intersectionTypeToString(it)
             );

    system(script);
}

void GEOGEBRA_generatePath(treeNode* nodeTOP, treeNode* nodeBOTTOM, int changeNumber) {
    GEOGEBRA_generating(nodeTOP, nodeBOTTOM, changeNumber, -1, -1);
}

void GEOGEBRA_generateTree(treeNode* node, int changeNumber) {
    GEOGEBRA_generating(node, NULL, changeNumber, -1, -1);
}

void GEOGEBRA_generateTreeWithIntersectionType(
        treeNode* node,
        int changeNumber,
        intersectionType it,
        short step
) {
    GEOGEBRA_generating(node, NULL, changeNumber, it, step);
}


void GEOGEBRA_printLxL(treeNode* source, treeNode* target) {
    printf("[GEOGEBRA] LxL\n");
    treeNode *current = target;
    while (current != source) {
      GEOGEBRA_printNode(current, 1, 1);
      current = current->parent;
    }
    GEOGEBRA_printNode(source, 1, 1);
}

void GEOGEBRA_printLxS(treeNode* source, treeNode* target) {
    printf("[GEOGEBRA printLxS] implement me...\n");
}

void GEOGEBRA_printSxL(treeNode* source, treeNode* target) {
    printf("[GEOGEBRA printLxS] implement me...\n");
}

void GEOGEBRA_printSxS(treeNode* source, treeNode* target) {
    printf("[GEOGEBRA printLxS] implement me...\n");
}

void GEOGEBRA_printCircle(const char* name, const double c[2], const double r) {
    printf("[GEOGEBRA] %s : %f² = ( x - %f )² + ( y - %f )²\n"
           , name, r, c[0], c[1]);
}

void GEOGEBRA_printLinePointPoint(const char* name, const double p1[2], const double p2[2]) {
    printf("[GEOGEBRA] %s : Line[(%f, %f), (%f, %f)]\n", name, p1[0], p1[1], p2[0], p2[1]);
}

void GEOGEBRA_printLinePointDir(const char* name, const double p[2], const double v[2]) {
    printf("[GEOGEBRA] %s : Line[(%f, %f), (%f, %f)]\n", name, p[0], p[1], p[0] + v[0], p[1] + v[1]);
}

void GEOGEBRA_printRay(const char* name, const double p[2], const double v[2]) {
    printf("[GEOGEBRA] %s : Ray[(%f, %f), (%f, %f)]\n", name, p[0], p[1], p[0] + v[0], p[1] + v[1]);
}

void GEOGEBRA_printPoint(const char* name, const double p[2]) {
    printf("[GEOGEBRA] %s = (%f, %f)\n", name, p[0], p[1]);
}

void GEOGEBRA_printRectangle(
    const char* name,
    const double a[2],
    const double b[2],
    const double c[2],
    const double e[2]
) {
    double A0[2] = { c[0] + e[0] * a[0] + e[1] * b[0],
                     c[1] + e[0] * a[1] + e[1] * b[1] };
    double B0[2] = { c[0] + e[0] * a[0] - e[1] * b[0],
                     c[1] + e[0] * a[1] - e[1] * b[1] };
    double C0[2] = { c[0] - e[0] * a[0] - e[1] * b[0],
                     c[1] - e[0] * a[1] - e[1] * b[1] };
    double D0[2] = { c[0] - e[0] * a[0] + e[1] * b[0],
                     c[1] - e[0] * a[1] + e[1] * b[1] };
    printf("[GEOGEBRA] %s : Polygon[(%f, %f), (%f, %f), (%f, %f), (%f, %f)]\n"
           , name, A0[0], A0[1], B0[0], B0[1], C0[0], C0[1], D0[0], D0[1]);
}


/// ------------------------------------------------------------------------------------------------
/// ------  |                                                                              |  ------
/// ------  |                        print interactive                                     |  ------
/// ------  v                                                                              v  ------
/// ------------------------------------------------------------------------------------------------


void GEOGEBRA_printInteractiveLoop2(const treeNode* node, FILE* fp) {
    char* fnName = "INTERACTIVE LOOP";
    char loopName = getNodeName(node);

    int configSize = node->childCount + 1;
    config *cfg = node->cfg;
    double* alphas = (double*) malloc(configSize * sizeof(double));

    double pairedAngle = getPairedAngle(node);
    for (int currentArc = 0; currentArc < configSize; currentArc++) {
        alphas[currentArc] = (getArcAngleDegree(cfg, currentArc) - pairedAngle) / (cfg->cfgArcs[currentArc]).numberOfArcSegments;
    }

    for (int i = 0; i < configSize; i++) {
        printf("[%s] %c backbones:%d angle:%f° alpha:%f°\n", fnName, loopName,
               node->cfg->cfgArcs[i].numberOfArcSegments,
               getArcAngleDegree(node->cfg, i),
               alphas[i]);
    }
//    fprintf(fp,
//    "           <expression label='loop_{%c}' exp='%f² = (x - (%f))² + (y - (%f))²'/>\n"
//    , name, l->r, l->c[0], l->c[1]
//    );


//    fprintf(fp,
//    "		<command name='Polygon'>\n"
//    "			<input a0='(%f, %f)' a1='(%f, %f)' a2='(%f, %f)' a3='(%f, %f)'/>\n"
//    "			<output a0='stem_{%c}' a1='l%c_1' a2='l%c_2' a3='l%c_3' a4='l%c_4'/>\n"
//    "		</command>\n"
//    , A0[0], A0[1], B0[0], B0[1], C0[0], C0[1], D0[0], D0[1]
//    , name, name, name, name, name
//    );

    char* labelOption = "";//" label=' '";

    int stemID = getNodeID(node);
    double PStem_id_C[2];
    getStemCenter(node, PStem_id_C);
    fprintf(fp, "<expression%s exp='PStem%dC=(%f, %f)'/>\n", labelOption, stemID, PStem_id_C[0], PStem_id_C[1]);
    fprintf(fp,
            "<element type=\"point\" label=\"PStem%dC\">"
            "<show object=\"true\" label=\"false\"/>"
            "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
            "</element>\n"
            , stemID);

    int loopID = getNodeID(node);
    double PLoop_id_C[2];
    getLoopCenter(node, PLoop_id_C);
    fprintf(fp, "<expression%s exp='PLoop%dC=(%f, %f)'/>\n", labelOption, loopID, PLoop_id_C[0], PLoop_id_C[1]);
    fprintf(fp,
            "<element type=\"point\" label=\"PLoop%dC\">"
            "<show object=\"true\" label=\"false\"/>"
            "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
            "</element>\n"
            , loopID);

    double radius = cfg->radius;
    double vLoopCenterToStemCenter[2];
    vector(PLoop_id_C, PStem_id_C, vLoopCenterToStemCenter);
    double vLoopCenterToStemCenterLength = vectorLength2D(vLoopCenterToStemCenter);
    double latestPoint[2];
    latestPoint[0] = PLoop_id_C[0] + radius / vLoopCenterToStemCenterLength * vLoopCenterToStemCenter[0];
    latestPoint[1] = PLoop_id_C[1] + radius / vLoopCenterToStemCenterLength * vLoopCenterToStemCenter[1];

    int latest = 0;
    fprintf(fp, "<expression exp='HelperLoop%dP%d=(%f, %f)'/>\n", loopID, latest, latestPoint[0], latestPoint[1]);
    fprintf(fp,
            "<element type=\"point\" label=\"HelperLoop%dP%d\">"
            "<show object=\"true\" label=\"false\"/>"
            "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
            "</element>\n"
            , loopID, latest);

    /// circle
    fprintf(fp,
            "<command name=\"Circle\">"
            "<input a0=\"PLoop%dC\" a1=\"HelperLoop%dP%d\"/>"
            "<output a0=\"circleLoop%d\"/>"
            "</command>\n"
            , loopID, loopID, latest
            , loopID);

    /// current loop's paired angle (in half)
    fprintf(fp,
            "<element type=\"angle\" label=\"helperLoop%dStemAngle\">"
            "<value val=\"%f\"/>"
            "<show object=\"false\" label=\"false\"/>"
            "<objColor r=\"0\" g=\"100\" b=\"0\" alpha=\"0.1\"/>"
            "<layer val=\"0\"/>"
            "<labelMode val=\"0\"/>"
            "<lineStyle thickness=\"2\" type=\"0\" typeHidden=\"1\"/>"
            "<arcSize val=\"30\"/>"
            "<animation step=\"1\" speed=\"1\" type=\"0\" playing=\"false\"/>"
            "</element>"
            , loopID
            , toRad(0.5 * pairedAngle));

    for (int i = 0; i < configSize; i++) {
        /// <command name="Rotate">
        ///     <input a0="HelperLoop65" a1="315*°" a2="PLoop65C"/>
        ///     <output a0="HelperLoop65'"/>
        /// </command>
        /// <element type="point" label="HelperLoop65'">
        ///     <show object="true" label="false"/>
        ///     <objColor r="0" g="0" b="255" alpha="0.0"/>
        ///     <layer val="0"/>
        ///     <labelMode val="0"/>
        ///     <coords x="-168.25455158676465" y="403.61802035897176" z="1.0"/>
        ///     <pointSize val="3"/>
        ///     <pointStyle val="0"/>
        /// </element>
        /// <command name="Angle">
        ///     <input a0="HelperLoop65'" a1="PLoop65C" a2="HelperLoop65"/>
        ///     <output a0="α"/>
        /// </command>
        /// <element type="angle" label="α">
        ///     <value val="0.7853981633974485"/>
        ///     <show object="true" label="true"/>
        ///     <objColor r="0" g="100" b="0" alpha="0.1"/>
        ///     <layer val="0"/><labelMode val="2"/>
        ///     <lineStyle thickness="2" type="0" typeHidden="1"/>
        ///     <arcSize val="30"/>
        ///     <allowReflexAngle val="true"/>
        /// </element>

        /// current arc's unpaired angle
        fprintf(fp,
                "<element type=\"angle\" label=\"helperLoop%dAlpha%d\">"
                "<value val=\"%f\"/>"
                "<show object=\"false\" label=\"false\"/>"
                "<objColor r=\"0\" g=\"100\" b=\"0\" alpha=\"0.1\"/>"
                "<layer val=\"0\"/>"
                "<labelMode val=\"0\"/>"
                "<lineStyle thickness=\"2\" type=\"0\" typeHidden=\"1\"/>"
                "<arcSize val=\"30\"/>"
                "<animation step=\"1\" speed=\"1\" type=\"0\" playing=\"false\"/>"
                "</element>"
                , loopID, i, toRad(alphas[i]));

        /// ---------------------
        /// arc bases - start
        /// ---------------------

        /// open arc from latest stem's center line
        fprintf(fp,
                "<command name=\"Rotate\">"
                "<input a0=\"HelperLoop%dP%d\" a1=\"-helperLoop%dStemAngle\" a2=\"PLoop%dC\"/>"
                "<output a0=\"HelperLoop%dP%d\"/>"
                "</command>"
                , loopID, latest, loopID, loopID
                , loopID, (latest + 1));
        fprintf(fp,
                "<element type=\"point\" label=\"HelperLoop%dP%d\">"
                "<show object=\"true\" label=\"false\"/>"
                "<objColor r=\"254\" g=\"178\" b=\"76\" alpha=\"0.0\"/>"
                "</element>\n"
                , loopID, (latest + 1));
        latest++;

        /// bases between stems
        for (int arcSegment = 0; arcSegment < (cfg->cfgArcs[i]).numberOfArcSegments; arcSegment++) {
            fprintf(fp,
                "<command name=\"Rotate\">"
                "<input a0=\"HelperLoop%dP%d\" a1=\"-helperLoop%dAlpha%d\" a2=\"PLoop%dC\"/>"
                "<output a0=\"HelperLoop%dP%d\"/>"
                "</command>\n"
                , loopID, latest, loopID, i, loopID
                , loopID, (latest+1));
            fprintf(fp,
                    "<element type=\"point\" label=\"HelperLoop%dP%d\">"
                    "<show object=\"true\" label=\"false\"/>"
                    "<objColor r=\"254\" g=\"178\" b=\"76\" alpha=\"0.0\"/>"
                    "</element>\n"
                    , loopID, (latest + 1));
            latest++;
        }

        /// close arc to next stem's center line
        fprintf(fp,
                "<command name=\"Rotate\">"
                "<input a0=\"HelperLoop%dP%d\" a1=\"-helperLoop%dStemAngle\" a2=\"PLoop%dC\"/>"
                "<output a0=\"HelperLoop%dP%d\"/>"
                "</command>"
                , loopID, latest, loopID, loopID
                , loopID, (latest + 1));
        fprintf(fp,
                "<element type=\"point\" label=\"HelperLoop%dP%d\">"
                "<show object=\"true\" label=\"false\"/>"
                "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
                "</element>\n"
                , loopID, (latest + 1));
        latest++;
        /// ---------------------
        /// arc bases - end
        /// ---------------------

    }

    short printStemAndLoop = 0;
    if (printStemAndLoop) {
        GEOGEBRA_generateLoop(fp, node);
        GEOGEBRA_generateStem(fp, node);
    }

    free(alphas);
}

void GEOGEBRA_printInteractivePoint_InteractionSpec(FILE* fp, const char* const pointName, double* coords) {
    fprintf(fp,
            "<element type=\"point\" label=\"%s\">"
            "<show object=\"true\" label=\"false\"/>"
            "<objColor r=\"49\" g=\"163\" b=\"84\" alpha=\"0.0\"/>"
            , pointName);

    if (coords != NULL) {
        fprintf(fp,
            "<coords x=\"%f\" y=\"%f\" z=\"1.0\"/>"
            , coords[0], coords[1]);
    }

    fprintf(fp,
            "</element>\n");
}

void GEOGEBRA_printInteractivePoint_HelperSpec(FILE* fp, const char* const pointName) {
    fprintf(fp,
            "<element type=\"point\" label=\"%s\">"
            "<show object=\"true\" label=\"false\"/>"
            "<selectionAllowed val=\"false\"/>"
            "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
            "</element>\n"
            , pointName);
}

void GEOGEBRA_printInteractivePoint_BaseSpec(FILE* fp, const char* const pointName) {
    fprintf(fp,
            "<element type=\"point\" label=\"%s\">"
            "<show object=\"true\" label=\"false\"/>"
            "<objColor r=\"254\" g=\"178\" b=\"76\" alpha=\"0.0\"/>"
            "</element>\n"
            , pointName);
}

#define POINT_TYPE_BASE   1
#define POINT_TYPE_HELPER 2
#define POINT_TYPE_INTERACTION 3

void GEOGEBRA_printInteractivePoint_finalize(FILE* fp, const char* const pointName, const int pointType) {
    switch (pointType) {
    case POINT_TYPE_BASE:
        GEOGEBRA_printInteractivePoint_BaseSpec(fp, pointName);
        break;
    case POINT_TYPE_HELPER:
        GEOGEBRA_printInteractivePoint_HelperSpec(fp, pointName);
        break;
    case POINT_TYPE_INTERACTION:
        GEOGEBRA_printInteractivePoint_InteractionSpec(fp, pointName, NULL);
        break;
    }

}

void GEOGEBRA_printInteractivePoint_PointVector(FILE* fp, const char* const outPointName, const char* const inPointName, const char* const inVectorName, const int pointType) {
    fprintf(fp,
            "<command name=\"Point\">"
            "<input a0=\"%s\" a1=\"%s\"/>"
            "<output a0=\"%s\"/>"
            "</command>\n"
            , inPointName, inVectorName
            , outPointName);
    GEOGEBRA_printInteractivePoint_finalize(fp, outPointName, pointType);
}

void GEOGEBRA_printInteractivePoint_Rotate(FILE* fp, const char* const outPointName, const char* const inPointName, const char* const center, const char* const angle, const int pointType) {
    fprintf(fp,
            "<command name=\"Rotate\">"
            "<input a0=\"%s\" a1=\"-%s\" a2=\"%s\"/>"
            "<output a0=\"%s\"/>"
            "</command>\n"
            , inPointName, angle, center
            , outPointName);
    GEOGEBRA_printInteractivePoint_finalize(fp, outPointName, pointType);
}

void GEOGEBRA_printInteractivePoint_coordinates2(FILE* fp, const char* const pointName, const double x, const double y, const int pointType) {
    fprintf(fp,
            "<expression exp=\"%s = (%f, %f)\"/>\n"
            , pointName, x, y);
    GEOGEBRA_printInteractivePoint_finalize(fp, pointName, pointType);
}

void GEOGEBRA_printInteractivePoint_coordinates1(FILE* fp, const char* const pointName, const double coords[2], const int pointType) {
    GEOGEBRA_printInteractivePoint_coordinates2(fp, pointName, coords[0], coords[1], pointType);
}

void GEOGEBRA_printInteractiveNode(
        FILE* fp,
        const treeNode* node
) {
    int loopID = getNodeID(node);
    char loopCenter[100];
    snprintf(loopCenter, sizeof(loopCenter), "P%05dLoopCenter", loopID);
    int stemID = getNodeID(node);
    char stemCenter[100];
    snprintf(stemCenter, sizeof(stemCenter), "P%05dStemCenter", stemID);

    char loopRadius[100];
    snprintf(loopRadius, sizeof(loopRadius), "cfg%05dradius", loopID);
    double radiusValue = node->cfg->radius;
    fprintf(fp,
            "<element type=\"numeric\" label=\"%s\">"
            "<value val=\"%f\"/>"
            "<show object=\"false\" label=\"false\"/>"
            "</element>"
            , loopRadius, radiusValue);

    char stemAngle[100];
    snprintf(stemAngle, sizeof(stemAngle), "helperLoop%05dStemAngle", loopID);
    fprintf(fp,
            "<expression exp=\"%s = (180 / π * arcsin(paired / (2 * %s)))°\"/>\n"
            , stemAngle, loopRadius);
    fprintf(fp,
            "<element type=\"angle\" label=\"%s\">"
            "<show object=\"false\" label=\"false\"/>"
            "</element>"
            , stemAngle);


    double pLoopCenter[2];
    getLoopCenter(node, pLoopCenter);
    double pStemCenter[2];
    getStemCenter(node, pStemCenter);
//    GEOGEBRA_printInteractivePoint_coordinates1(fp, stemCenter, pStemCenter, POINT_TYPE_INTERACTION);

    char extA[100], extB[100];
    snprintf(extA, sizeof(extA), "stem%05dextA", stemID);
    snprintf(extB, sizeof(extB), "stem%05dextB", stemID);
    fprintf(fp,
            "<element type=\"numeric\" label=\"%s\">"
            "<value val=\"%f\"/>"
            "<show object=\"false\" label=\"false\"/>"
            "</element>"
            , extA, node->sBox->e[0]);
    fprintf(fp,
            "<element type=\"numeric\" label=\"%s\">"
            "<value val=\"%f\"/>"
            "<show object=\"false\" label=\"false\"/>"
            "</element>"
            , extB, node->sBox->e[1]);

    char initDirA[100], initDirB[100];
    snprintf(initDirA, sizeof(initDirA), "stem%05dinitDirA", stemID);
    snprintf(initDirB, sizeof(initDirB), "stem%05dinitDirB", stemID);
    fprintf(fp,
            "<element type=\"vector\" label=\"%s\">"
            "<show object=\"false\" label=\"false\"/>"
            "<coords x=\"%f\" y=\"%f\" z=\"0.0\"/>"
            "<coordStyle style=\"cartesian\"/>"
            "</element>\n"
            , initDirA
            , node->sBox->a[0], node->sBox->a[1]);
    fprintf(fp,
            "<element type=\"vector\" label=\"%s\">"
            "<show object=\"false\" label=\"false\"/>"
            "<coords x=\"%f\" y=\"%f\" z=\"0.0\"/>"
            "<coordStyle style=\"cartesian\"/>"
            "</element>\n"
            , initDirB
            , node->sBox->b[0], node->sBox->b[1]);

    treeNode* parent = getParent(node);
    int parentLoopID = isExterior(parent) ? 0 : getNodeID(parent);
    char parentLoopRadius[100];
    snprintf(parentLoopRadius, sizeof(parentLoopRadius), "cfg%05dradius", parentLoopID);
    char parentLoopCenter[100];
    snprintf(parentLoopCenter, sizeof(parentLoopCenter), "P%05dLoopCenter", parentLoopID);

    char stemStart[100];
    snprintf(stemStart, sizeof(stemStart), "P%05dStemStart", stemID);

    if (!isExterior(parent)) {

        char parentStemAngle[100];
        snprintf(parentStemAngle, sizeof(parentStemAngle), "helperLoop%05dStemAngle", parentLoopID);
        /// circle
        char radius[100];
        snprintf(radius, sizeof(radius), "radius%05dStemCenter", stemID);
        fprintf(fp,
                "<expression exp=\"%s = %s + (%s * cos(%s))\"/>\n"
                , radius, extA, parentLoopRadius, parentStemAngle);
        char circle[100];
        snprintf(circle, sizeof(circle), "circle%05dStemCenter", stemID);
        fprintf(fp,
                "<command name=\"Circle\">"
                "<input a0=\"%s\" a1=\"%s\"/>"
                "<output a0=\"%s\"/>"
                "</command>\n"
                , parentLoopCenter, radius
                , circle);
        fprintf(fp,
                "<element type=\"conic\" label=\"%s\">"
                "<show object=\"false\" label=\"false\"/>"
                //"<objColor r=\"255\" g=\"102\" b=\"204\" alpha=\"0.0\"/>"
                //"<lineStyle thickness=\"5\" type=\"0\" typeHidden=\"1\"/>"
                "</element>\n"
                , circle);

        fprintf(fp,
                "<command name=\"Point\">"
                "<input a0=\"%s\"/>"
                "<output a0=\"%s\"/>"
                "</command>\n"
                , circle
                , stemCenter);
        GEOGEBRA_printInteractivePoint_InteractionSpec(fp, stemCenter, pStemCenter);
//        fprintf(fp,
//                "<element type=\"point\" label=\"%s\">"
//                "<show object=\"true\" label=\"false\"/>"
//                "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
//                "<coords x=\"%f\" y=\"%f\" z=\"1.0\"/>"
//                "</element>\n"
//                , stemCenter, pStemCenter[0], pStemCenter[1]);

        char vParentLoopCenterToStemCenter[100];
        snprintf(vParentLoopCenterToStemCenter, sizeof(vParentLoopCenterToStemCenter), "vParentLoopCenterTo%05dStemCenter", stemID);
        fprintf(fp,
                "<command name=\"Vector\">"
                "<input a0=\"%s\" a1=\"%s\"/>"
                "<output a0=\"%s\"/>"
                "</command>"
                , parentLoopCenter, stemCenter
                , vParentLoopCenterToStemCenter);
        fprintf(fp,
                "<element type=\"vector\" label=\"%s\">"
                "<show object=\"false\" label=\"false\"/>"
                "</element>\n"
                , vParentLoopCenterToStemCenter);

        fprintf(fp,
                "<command name=\"Length\">"
                "<input a0=\"%s\"/>"
                "<output a0=\"%sLength\"/>"
                "</command>"
                , vParentLoopCenterToStemCenter
                , vParentLoopCenterToStemCenter);

        char vStemCenterToLoopCenter[100];
        snprintf(vStemCenterToLoopCenter, sizeof(vStemCenterToLoopCenter), "v%05dStemCenterToLoopCenter", loopID);
        fprintf(fp,
                "<expression exp=\"%s = %s * ((%s + %s * cos(%s)) / %sLength)\"/>"
                , vStemCenterToLoopCenter, vParentLoopCenterToStemCenter, extA, loopRadius, stemAngle, vParentLoopCenterToStemCenter);
        fprintf(fp,
                "<element type=\"vector\" label=\"%s\">"
                "<show object=\"false\" label=\"false\"/>"
                "</element>\n"
                , vStemCenterToLoopCenter);
        GEOGEBRA_printInteractivePoint_PointVector(fp, loopCenter, stemCenter, vStemCenterToLoopCenter, POINT_TYPE_HELPER);
//        fprintf(fp,
//                "<command name=\"Point\">"
//                "<input a0=\"%s\" a1=\"%s\"/>"
//                "<output a0=\"%s\"/>"
//                "</command>\n"
//                , stemCenter, vStemCenterToLoopCenter
//                , loopCenter);




    } else {
//        GEOGEBRA_printInteractivePoint_coordinates1(fp, stemCenter, pStemCenter, POINT_TYPE_INTERACTION);
        fprintf(fp,
                "<expression exp=\"exterior%05d : y = %f + %s\"/>\n"
                , stemID, EXTERIOR_Y, extA);
        fprintf(fp,
                "<element type=\"line\" label=\"exterior%05d\">"
                "<show object=\"false\" label=\"false\"/>"
                "</element>\n"
                , stemID);
        fprintf(fp,
                "<command name=\"Point\">"
                "<input a0=\"exterior%05d\"/>"
                "<output a0=\"%s\"/>"
                "</command>\n"
                , stemID
                , stemCenter);
        GEOGEBRA_printInteractivePoint_InteractionSpec(fp, stemCenter, pStemCenter);

        GEOGEBRA_printInteractivePoint_coordinates1(fp, loopCenter, pLoopCenter, POINT_TYPE_HELPER);
    }



}

void GEOGEBRA_printInteractiveLoop_Bases(FILE* fp, const treeNode* node) {
    char* fnName = "INTERACTIVE BASES";

    int loopID = getNodeID(node);
    char loopCenter[100];
    snprintf(loopCenter, sizeof(loopCenter), "P%05dLoopCenter", loopID);
    int stemID = getNodeID(node);
    char stemCenter[100];
    snprintf(stemCenter, sizeof(stemCenter), "P%05dStemCenter", stemID);

    int configSize = node->childCount + 1;
    config *cfg = node->cfg;
    double* alphas = (double*) malloc(configSize * sizeof(double));
    double pairedAngle = getPairedAngle(node);
    for (int currentArc = 0; currentArc < configSize; currentArc++) {
        alphas[currentArc] = (getArcAngleDegree(cfg, currentArc) - pairedAngle) / (cfg->cfgArcs[currentArc]).numberOfArcSegments;
    }

//    double radiusValue = node->cfg->radius;
    char radius[100];
    snprintf(radius, sizeof(radius), "cfg%05dradius", loopID);
//    fprintf(fp,
//            "<element type=\"numeric\" label=\"%s\">"
//            "<value val=\"%f\"/>"
//            "<show object=\"false\" label=\"false\"/>"
//            "</element>"
//            , radius, radiusValue);

    int latest = 0;

    /// vector von loop zu stem (centers)
    fprintf(fp,
            "<command name=\"Vector\">"
            "<input a0=\"%s\" a1=\"%s\"/>"
            "<output a0=\"v%05dLoopStemVector\"/>"
            "</command>"
            , loopCenter, stemCenter
            , loopID
            );
    fprintf(fp,
            "<element type=\"vector\" label=\"v%05dLoopStemVector\">"
            "<show object=\"false\" label=\"false\"/>"
            "</element>\n"
            , loopID
            );

    /// länge des vectors
    fprintf(fp,
            "<command name=\"Length\">"
            "<input a0=\"v%05dLoopStemVector\"/>"
            "<output a0=\"v%05dLoopStemLength\"/>"
            "</command>"
            , loopID
            , loopID);

    fprintf(fp,
            "<expression exp=\"v%05dLoopFirstHelperPoint = (%s / v%05dLoopStemLength) * v%05dLoopStemVector\"/>"
            , loopID, radius, loopID, loopID
            );
    fprintf(fp,
            "<element type=\"vector\" label=\"v%05dLoopFirstHelperPoint\">"
            "<show object=\"false\" label=\"false\"/>"
            "</element>\n"
            , loopID
            );

    char helperThis[100];
    char helperNext[100];
    snprintf(helperThis, sizeof(helperThis), "H%05dP%03d", loopID, latest);
    char vLoopStem[100];
    snprintf(vLoopStem, sizeof(vLoopStem), "v%05dLoopFirstHelperPoint", loopID);
    GEOGEBRA_printInteractivePoint_PointVector(fp, helperThis, loopCenter, vLoopStem, POINT_TYPE_HELPER);

//    fprintf(fp, "<expression exp='H%05dP%03d=(%f, %f)'/>\n", loopID, latest, latestPoint[0], latestPoint[1]);
//    fprintf(fp,
//            "<element type=\"point\" label=\"H%05dP%03d\">"
//            "<show object=\"true\" label=\"false\"/>"
//            "<objColor r=\"189\" g=\"189\" b=\"189\" alpha=\"0.0\"/>"
//            "</element>\n"
//            , loopID, latest);

//    /// circle
//    char circle[100];
//    snprintf(circle, sizeof(circle), "loop%05dcircle", loopID);
//    fprintf(fp,
//            "<command name=\"Circle\">"
//            "<input a0=\"%s\" a1=\"%s\"/>"
//            "<output a0=\"%s\"/>"
//            "</command>\n"
//            , loopCenter, helperThis
//            , circle);
//    fprintf(fp,
//            "<element type=\"conic\" label=\"%s\">"
//            "<show object=\"true\" label=\"false\"/>"
//            "<objColor r=\"255\" g=\"102\" b=\"204\" alpha=\"0.0\"/>"
////            "<lineStyle thickness=\"2\" type=\"20\" typeHidden=\"1\"/>"
//            "<lineStyle thickness=\"5\" type=\"0\" typeHidden=\"1\"/>"
//            "</element>"
//            , circle);

    /// current loop's paired angle (in half)
    char stemAngle[100];
    snprintf(stemAngle, sizeof(stemAngle), "helperLoop%05dStemAngle", loopID);
//    fprintf(fp,
//            "<expression exp=\"%s = (180 / π * arcsin(paired / (2 * %s)))°\"/>\n"
//            , stemAngle, radius);
//    fprintf(fp,
//            "<element type=\"angle\" label=\"%s\">"
//            "<show object=\"false\" label=\"false\"/>"
//            "</element>"
//            , stemAngle);
    int i;
    for (i = 0; i < configSize; i++) {

        int idLeft = (i == 0) ? getNodeID(node) : getNodeID(getChild(node, i-1));
        int idRight = (i == configSize - 1) ? getNodeID(node) : getNodeID(getChild(node, i));
        char refStemLeft[100];
        char refStemRight[100];
        char arcAngle[100];
        snprintf(refStemLeft, sizeof(refStemLeft), "P%05dStemCenter", idLeft);
        snprintf(refStemRight, sizeof(refStemRight), "P%05dStemCenter", idRight);
        snprintf(arcAngle, sizeof(arcAngle), "cfg%05dArc%03dAngle", loopID, i);
        fprintf(fp,
                "<command name=\"Angle\">"
                "<input a0=\"%s\" a1=\"%s\" a2=\"%s\"/>"
                "<output a0=\"%s\"/>"
                "</command>\n"
                , refStemRight, loopCenter, refStemLeft
                , arcAngle);
        fprintf(fp,
                "<element type=\"angle\" label=\"%s\">"
                "<show object=\"false\" label=\"false\"/>"
                "</element>"
                , arcAngle);

        /// current arc's backbone count
        char backboneCount[100];
        snprintf(backboneCount, sizeof(backboneCount), "cfg%05dArc%03dBackbones", loopID, i);
        fprintf(fp,
                "<element type=\"numeric\" label=\"%s\">"
                "<value val=\"%d\"/>"
                "<show object=\"false\" label=\"false\"/>"
                "</element>"
                , backboneCount, (cfg->cfgArcs[i]).numberOfArcSegments);

        /// current arc's unpaired angle
        char alpha[100];
        snprintf(alpha, sizeof(alpha), "h%05dAlpha%03d", loopID, i);
        fprintf(fp,
                "<command name=\"If\">"
                "<input a0=\"%s == 0°\" a1=\"(360° - 2 * %s) / %s\" a2=\"(%s - 2 * %s) / %s\"/>"
                "<output a0=\"%s\"/>"
                "</command>\n"
                , arcAngle, stemAngle, backboneCount, arcAngle, stemAngle, backboneCount
                , alpha);
        fprintf(fp,
                "<element type=\"angle\" label=\"%s\">"
                "<show object=\"false\" label=\"false\"/>"
                "</element>"
                , alpha);

        /// ---------------------
        /// arc bases - start
        /// ---------------------

        /// open arc from latest stem's center line
        snprintf(helperThis, sizeof(helperThis), "H%05dP%03d", loopID, latest);
        snprintf(helperNext, sizeof(helperNext), "H%05dP%03d", loopID, (latest+1));
        GEOGEBRA_printInteractivePoint_Rotate(fp, helperNext, helperThis, loopCenter, stemAngle, POINT_TYPE_BASE);
        latest++;

        /// bases between stems
        for (int arcSegment = 0; arcSegment < (cfg->cfgArcs[i]).numberOfArcSegments; arcSegment++) {
            snprintf(helperThis, sizeof(helperThis), "H%05dP%03d", loopID, latest);
            snprintf(helperNext, sizeof(helperNext), "H%05dP%03d", loopID, (latest+1));
            GEOGEBRA_printInteractivePoint_Rotate(fp, helperNext, helperThis, loopCenter, alpha, POINT_TYPE_BASE);

            fprintf(fp,
                    "<command name=\"CircleArc\">"
                    "<input a0=\"%s\" a1=\"%s\" a2=\"%s\"/>"
                    "<output a0=\"vis%05dL%03dArc%03d\"/>"
                    "</command>"
                    , loopCenter, helperNext, helperThis
                    , loopID, i, arcSegment);
//            fprintf(fp,
//                    "<command name=\"Segment\">"
//                    "<input a0=\"%s\" a1=\"%s\"/>"
//                    "<output a0=\"vis%05dL%03dArc%03d\"/>"
//                    "</command>"
//                    , helperThis, helperNext
//                    , loopID, i, arcSegment);

            latest++;
        }

        /// close arc to next stem's center line
        snprintf(helperThis, sizeof(helperThis), "H%05dP%03d", loopID, latest);
        snprintf(helperNext, sizeof(helperNext), "H%05dP%03d", loopID, (latest+1));
        GEOGEBRA_printInteractivePoint_Rotate(fp, helperNext, helperThis, loopCenter, stemAngle, POINT_TYPE_HELPER);
        latest++;
        /// ---------------------
        /// arc bases - end
        /// ---------------------

    }

    short printStemAndLoop = 0;
    if (printStemAndLoop) {
        GEOGEBRA_generateLoop(fp, node);
        GEOGEBRA_generateStem(fp, node);
    }

    free(alphas);
}

void GEOGEBRA_printInteractiveTree_LoopsAndStems(const treeNode* node, FILE* fp) {
    int i;
    for (i = 0; i < node->childCount; i++) {
        const treeNode* child = getChild(node, i);
        GEOGEBRA_printInteractiveNode(fp, child);

        GEOGEBRA_printInteractiveTree_LoopsAndStems(child, fp);
    }
}

void GEOGEBRA_printInteractiveTree_Bases(const treeNode* node, FILE* fp) {
    int i;
    for (i = 0; i < node->childCount; i++) {
        const treeNode* child = getChild(node, i);
        GEOGEBRA_printInteractiveLoop_Bases(fp, child);

        GEOGEBRA_printInteractiveTree_Bases(child, fp);
    }
}

void GEOGEBRA_printInteractiveTree(const treeNode* tree) {
    FILE* fp_xml;
    fp_xml = fopen("geogebra.xml", "w");

//    GEOGEBRA_generateHead(fp_xml);
//    // get bounds and adjust view
//    GEOGEBRA_adjustView(fp_xml, tree);
//    GEOGEBRA_generateHeadEnd(fp_xml);
    /// -----------------------------------------------------------------------
    fprintf(fp_xml,
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
            "<geogebra format=\"4.0\" id=\"cd099ee0-1acf-462b-8863-25eb5684fe6a\"  xsi:noNamespaceSchemaLocation=\"http://www.geogebra.org/ggb.xsd\" xmlns=\"\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" >\n"
            "<gui>\n"
            "        <window width=\"1215\" height=\"776\" />\n"
            "        <perspectives>\n"
            "<perspective id=\"tmp\">\n"
            "        <panes>\n"
            "                <pane location=\"\" divider=\"0.3390946502057613\" orientation=\"1\" />\n"
            "        </panes>\n"
            "        <views>\n"
            "                <view id=\"4\" toolbar=\"0 || 2020 , 2021 , 2022 , 66 || 2001 , 2003 , 2002 , 2004 , 2005 || 2040 , 2041 , 2042 , 2044 , 2043\" visible=\"false\" inframe=\"false\" stylebar=\"false\" location=\"1,1\" size=\"300\" window=\"100,100,600,400\" />\n"
            "                <view id=\"1\" visible=\"true\" inframe=\"false\" stylebar=\"false\" location=\"1\" size=\"793\" window=\"100,100,600,400\" />\n"
            "                <view id=\"2\" visible=\"true\" inframe=\"false\" stylebar=\"false\" location=\"3\" size=\"412\" window=\"100,100,250,400\" />\n"
            "                <view id=\"16\" visible=\"false\" inframe=\"true\" stylebar=\"false\" location=\"1\" size=\"150\" window=\"50,50,500,500\" />\n"
            "                <view id=\"32\" visible=\"false\" inframe=\"true\" stylebar=\"true\" location=\"1\" size=\"150\" window=\"50,50,500,500\" />\n"
            "                <view id=\"64\" visible=\"false\" inframe=\"true\" stylebar=\"false\" location=\"1\" size=\"150\" window=\"50,50,500,500\" />\n"
            "        </views>\n"
            "        <toolbar show=\"true\" items=\"0 39 59 || 1 501 67 , 5 19 , 72 | 2 15 45 , 18 65 , 7 37 | 4 3 8 9 , 13 44 , 58 , 47 || 16 51 64 , 70 | 10 34 53 11 , 24  20 22 , 21 23 | 55 56 57 , 12 || 36 46 , 38 49 50 , 71 | 30 29 54 32 31 33 | 17 26 62 , 14 66 68 | 25 52 60 61 || 40 41 42 , 27 28 35 , 6\" />\n"
            "        <input show=\"true\" cmd=\"true\" top=\"false\" />\n"
            "</perspective>\n"
            "        </perspectives>\n"
            "        <labelingStyle  val=\"3\"/>\n"
            "        <font  size=\"12\"/>\n"
            "        <graphicsSettings javaLatexFonts=\"false\"/>\n"
            "        <consProtColumns  col0=\"true\" col1=\"true\" col2=\"false\" col3=\"true\" col4=\"false\" col5=\"true\" col6=\"true\" col7=\"false\"/>\n"
            "        <consProtocol useColors=\"true\" addIcons=\"false\" showOnlyBreakpoints=\"false\"/>\n"
            "</gui>\n"
            "<euclidianView>\n"
            "        <size  width=\"793\" height=\"651\"/>\n"
            "        <coordSystem xZero=\"478.685766\" yZero=\"512.750691\" scale=\"0.072757\" yscale=\"0.07275699999999999\"/>\n"
            "        <evSettings axes=\"false\" grid=\"false\" gridIsBold=\"false\" pointCapturing=\"3\" rightAngleStyle=\"1\" checkboxSize=\"13\" gridType=\"0\"/>\n"
            "        <bgColor r=\"255\" g=\"255\" b=\"255\"/>\n"
            "        <axesColor r=\"0\" g=\"0\" b=\"0\"/>\n"
            "        <gridColor r=\"192\" g=\"192\" b=\"192\"/>\n"
            "        <lineStyle axes=\"1\" grid=\"10\"/>\n"
            "        <axis id=\"0\" show=\"false\" label=\"\" unitLabel=\"\" tickStyle=\"1\" showNumbers=\"true\"/>\n"
            "        <axis id=\"1\" show=\"false\" label=\"\" unitLabel=\"\" tickStyle=\"1\" showNumbers=\"true\"/>\n"
            "</euclidianView>\n"
            "<kernel>\n"
            "        <continuous val=\"false\"/>\n"
            "        <decimals val=\"2\"/>\n"
            "        <angleUnit val=\"degree\"/>\n"
            "        <algebraStyle val=\"0\"/>\n"
            "        <coordStyle val=\"0\"/>\n"
            "        <angleFromInvTrig val=\"false\"/>\n"
            "</kernel>\n"
            "<scripting blocked=\"false\" disabled=\"false\"/>\n"
            "<construction title=\"\" author=\"\" date=\"\">\n"
            );
    /// -----------------------------------------------------------------------

    fprintf(fp_xml, "<expression exp=\"paired = 35\"/>\n");
    fprintf(fp_xml, "<expression exp=\"unpaired = 25\"/>\n");
    fprintf(fp_xml, "<expression exp=\"exterior : y = %f\"/>\n", EXTERIOR_Y);

    GEOGEBRA_printInteractiveTree_LoopsAndStems(tree, fp_xml);
    GEOGEBRA_printInteractiveTree_Bases(tree, fp_xml);

    GEOGEBRA_generateTail(fp_xml);
    fclose(fp_xml);

    FILE* fp_js;
    fp_js = fopen("geogebra_javascript.js", "w");
    fprintf(fp_js, "function ggbOnInit() {}\n");
    fclose(fp_js);

    char script[1000];
    snprintf(script, sizeof(script),
            ""
            "/bin/bash \n"
            "zip -m GeoGebra_live_example.zip geogebra.xml geogebra_javascript.js \n"
            "mv GeoGebra_live_example.zip GeoGebra_live_example.ggb \n"
            ""
             );

    system(script);
}
