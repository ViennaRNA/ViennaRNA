#ifndef RNAPUZZLER_INTERSECTLEVELBOUNDINGBOXES_H
#define RNAPUZZLER_INTERSECTLEVELBOUNDINGBOXES_H

/*
 *      RNApuzzler intersect bounding boxes
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

#include <stdlib.h>
#include <math.h>

#include "intersectLevelLines.inc"
#include "boundingBoxes.inc"
#include "configtree.inc"
#include "vector_math.inc"
#include "definitions.inc"

/**
 * @brief intersectStemStem
 * Checks for intersection of two given stems.
 * @param stem1
 * @param stem2
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short
intersectStemStem(const stemBox *stem1,
                  const stemBox *stem2);


/**
 * @brief intersectLoopLoop
 * Checks for intersection of two given loops.
 * @param loop1
 * @param loop2
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short
intersectLoopLoop(const loopBox *loop1,
                  const loopBox *loop2);


/**
 * @brief intersectStemLoop
 * Check for intersection of given stem and given loop.
 * @param stem
 * @param loop
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short
intersectStemLoop(const stemBox *stem,
                  const loopBox *loop);


/**
 * @brief intersectLoopBulges
 * Check for intersection of given loop and given stem's bulges.
 * @param loop
 * @param stem
 * @param bulge
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short
intersectLoopBulges(const loopBox *loop,
                    const stemBox *stem,
                    int           *bulge);


/**
 * @brief intersectBulgesBulges
 * @return
 */
PRIVATE short
intersectBulgesBulges(const stemBox *stem1,
                      const stemBox *stem2,
                      int           *bulge1,
                      int           *bulge2);


/**
 * @brief intersectStemBulges
 * @param stem1
 * @param stem2
 * @param bulge2
 * @return
 */
PRIVATE short
intersectStemBulges(const stemBox *stem1,
                    const stemBox *stem2,
                    int           *bulge2);


PRIVATE short
intersectLoopLoop(const loopBox *loop1,
                  const loopBox *loop2)
{
  double  c1[2] = {
    loop1->c[0], loop1->c[1]
  };
  double  r1    = loop1->r + 0.5 * epsilonRecognize;
  double  c2[2] = {
    loop2->c[0], loop2->c[1]
  };
  double  r2 = loop2->r + 0.5 * epsilonRecognize;

  return intersectCircleCircle(c1, r1, c2, r2);
}


PRIVATE void
projectPointOntoLine(const double a[2],
                     const double b[2],
                     const double p[2],
                     double       ret_p[2])
{
  double  u[2];

  vector(a, p, u);
  double  v[2];
  vector(a, b, v);
  double  w[2] = {
    -v[1], v[0]
  };

  /* compute r which grants all needed information */
  double  r = (u[1] - u[0] * w[1] / w[0]) / (v[1] - v[0] * w[1] / w[0]);

  if (r < 0.0) {
    ret_p[0]  = a[0];
    ret_p[1]  = a[1];
  } else if (r > 1.0) {
    ret_p[0]  = b[0];
    ret_p[1]  = b[1];
  } else {
    ret_p[0]  = a[0] + r * v[0];
    ret_p[1]  = a[1] + r * v[1];
  }

  return;
}


PRIVATE void
ClosestPtPointBulge(const double  p[2],
                    const double  a[2],
                    const double  b[2],
                    const double  c[2],
                    double        ret_p[2])
{
  char *fnName = "CLOSEST ON BULGE";

  /*
   * Note:
   *
   * In contrast to ClosestPtPointTriangle (taken from book Real-Time Collision Detection)
   * this function does not work for general triangles.
   * We implicitely make use of the fact that our bulges are equiliteral triangles
   * or to be more precise there are no angles greater than 90° in our triangles.
   * This allows for only checking for the first occurance of a side where point p
   * is on the outer side of the triangle and not checking any further.
   *
   * In fact applying this function to a triangle with some angle being greater than 90°
   * may lead to wrong results as follows as described in the mentioned book
   * (at ClosestPtPointTriangle) as well.
   */


  /* check if p on outer side of AB */
  short orientABC = isToTheRightPointPoint(a, b, c);
  short orientABP = isToTheRightPointPoint(a, b, p);

  if (orientABC != orientABP) {
    /* p is on outer side of AB */
    projectPointOntoLine(a, b, p, ret_p);
    return;
  }

  /* check if p on outer side of BC */
  short orientBCA = isToTheRightPointPoint(b, c, a);
  short orientBCP = isToTheRightPointPoint(b, c, p);

  if (orientBCA != orientBCP) {
    /* p is on outer side of BC */
    projectPointOntoLine(b, c, p, ret_p);
    return;
  }

  /* check if p on outer side of CA */
  short orientCAB = isToTheRightPointPoint(c, a, b);
  short orientCAP = isToTheRightPointPoint(c, a, p);

  if (orientCAB != orientCAP) {
    /* p is on outer side of CA */
    projectPointOntoLine(c, a, p, ret_p);
    return;
  }

  /* p is inside ABC */
  ret_p[0]  = p[0];
  ret_p[1]  = p[1];
  return;
}


PRIVATE void
ClosestPtPointTriangle(const double p[2],
                       const double a[2],
                       const double b[2],
                       const double c[2],
                       double       ret_p[2])
{
  /* Check if P in vertex region outside A */
  double  ab[2];
  double  ac[2];
  double  ap[2];

  vector(a, b, ab);
  vector(a, c, ac);
  vector(a, p, ap);
  double  d1  = scalarProduct2D(ab, ap);
  double  d2  = scalarProduct2D(ac, ap);

  if (d1 <= 0.0 && d2 <= 0.0) {
    /* barycentric coordinates (1,0,0) */
    ret_p[0]  = a[0];
    ret_p[1]  = a[1];
    return;
  }

  /* Check if P in vertex region outside B */
  double  bp[2];
  vector(b, p, bp);
  double  d3  = scalarProduct2D(ab, bp);
  double  d4  = scalarProduct2D(ac, bp);

  if (d3 >= 0.0 && d4 <= 0.0) {
    /* barycentric coordinates (0,1,0) */
    ret_p[0]  = b[0];
    ret_p[1]  = b[1];
    return;
  }

  /* Check if P in edge region of AB, if so return projection of P onto AB */
  double vc = d1 * d4 - d3 * d2;

  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    /* barycentric coordinates (1-v,v,0) */
    double v = d1 / (d1 - d3);
    ret_p[0]  = a[0] + v * ab[0];
    ret_p[1]  = a[1] + v * ab[1];
    return;
  }

  /* Check if P in vertex region outside C */
  double  cp[2];
  vector(c, p, cp);
  double  d5  = scalarProduct2D(ab, cp);
  double  d6  = scalarProduct2D(ac, cp);

  if (d6 >= 0.0 && d5 <= d6) {
    /* barycentric coordinates (0,0,1) */
    ret_p[0]  = c[0];
    ret_p[1]  = c[1];
    return;
  }

  /* Check if P in edge region of AC, if so return projection of P onto AC */
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    /* barycentric coordinates (1-w,0,w) */
    double w = d2 / (d2 - d6);
    ret_p[0]  = a[0] + w * ac[0];
    ret_p[1]  = a[1] + w * ac[1];
    return;
  }

  /* Check if P in edge region of BC, if so return projection of P onto BC */
  double va = d3 * d6 - d5 * d4;

  if (va <= 0.0) {
    if ((d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
      /* barycentric coordinates (0,1-w,w) */
      double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      ret_p[0]  = b[0] + w * (c[0] - b[0]);
      ret_p[1]  = b[1] + w * (c[1] - b[1]);
      return;
    } else if ((d4 - d3) <= 0.0 && (d5 - d6) >= 0.0) {
      /* barycentric coordinates (0,1-w,w) */
      double w = (d3 - d4) / ((d3 - d4) + (d5 - d6));
      ret_p[0]  = b[0] + w * (c[0] - b[0]);
      ret_p[1]  = b[1] + w * (c[1] - b[1]);
      return;
    }
  }

  /* P inside face region. Compute Q through its barycentric coordinates (u,v,w) */
  double  denom = 1.0 / (va + vb + vc);
  double  v     = vb * denom;
  double  w     = vc * denom;
  ret_p[0]  = a[0] + ab[0] * v + ac[0] * w;
  ret_p[1]  = a[1] + ab[1] * v + ac[1] * w;
  return;
}


/**
 * @brief ClosestPtPointOBB
 * Implementation of ClostestPtPointOBB from the book "Real Time Collision Detection".
 * Calculates the point on a rectangle (stem) that is closest to a given point p.
 * @param stem      - rectangle / stem bounding box
 * @param p_x       - x value of point p
 * @param p_y       - y value of point p
 * @param ret_p_x   - pointer to closest point's x value (used as return)
 * @param ret_p_y   - pointer to closest point's y value (used as return)
 */
PRIVATE void
ClosestPtPointOBB(const stemBox stem,
                  const double  p[2],
                  double        ret_p[2])
{
  double  u0[2] = {
    stem.a[0], stem.a[1]
  };
  double  u1[2] = {
    stem.b[0], stem.b[1]
  };
  double  dv[2];

  vector(stem.c, p, dv);

  double  dist_0  = scalarProduct2D(dv, u0);
  double  dist_1  = scalarProduct2D(dv, u1);

  short   sign_d0 = dist_0 < 0 ? -1 : 1;
  short   sign_d1 = dist_1 < 0 ? -1 : 1;
  short   sign_e0 = stem.e[0] < 0 ? -1 : 1;
  short   sign_e1 = stem.e[1] < 0 ? -1 : 1;
  double  abs_d0  = sign_d0 * dist_0;
  double  abs_d1  = sign_d1 * dist_1;
  double  abs_e0  = sign_e0 * stem.e[0];
  double  abs_e1  = sign_e1 * stem.e[1];

  /* clamp dist_0/1 to the extents of the OBB */
  double  clamped_0 = abs_d0 > abs_e0 ? sign_d0 * abs_e0 : sign_d0 * abs_d0;
  double  clamped_1 = abs_d1 > abs_e1 ? sign_d1 * abs_e1 : sign_d1 * abs_d1;

  ret_p[0]  = stem.c[0] + clamped_0 * stem.a[0] + clamped_1 * stem.b[0];
  ret_p[1]  = stem.c[1] + clamped_0 * stem.a[1] + clamped_1 * stem.b[1];
}


PRIVATE short
intersectStemLoop(const stemBox *stem,
                  const loopBox *loop)
{
  short   intersect = 0;

  double  p[2];

  ClosestPtPointOBB(*stem, loop->c, p);

  double  v_c_to_p[2];
  vector(loop->c, p, v_c_to_p);

  double  distanceSquared = vectorLength2DSquared(v_c_to_p);
  intersect = (distanceSquared <
               ((loop->r + epsilonRecognize) *
                (loop->r + epsilonRecognize)));
  /* add epsilon to classify nearby objects as intersecting */
  return intersect;
}


PRIVATE short
intersectStemStem(const stemBox *stem1,
                  const stemBox *stem2)
{
  /* brute force approach for intersecting two rectangles */
  double  stem1_ea[2] = {
    stem1->e[0] * stem1->a[0],
    stem1->e[0] * stem1->a[1]
  };
  double  stem1_eb[2] = {
    stem1->e[1] * stem1->b[0],
    stem1->e[1] * stem1->b[1]
  };
  double  B1[2] = {
    stem1->c[0] + stem1_ea[0] + stem1_eb[0],
    stem1->c[1] + stem1_ea[1] + stem1_eb[1]
  };
  double  C1[2] = {
    stem1->c[0] + stem1_ea[0] - stem1_eb[0],
    stem1->c[1] + stem1_ea[1] - stem1_eb[1]
  };
  double  D1[2] = {
    stem1->c[0] - stem1_ea[0] - stem1_eb[0],
    stem1->c[1] - stem1_ea[1] - stem1_eb[1]
  };
  double  A1[2] = {
    stem1->c[0] - stem1_ea[0] + stem1_eb[0],
    stem1->c[1] - stem1_ea[1] + stem1_eb[1]
  };

  double  stem2_ea[2] = {
    stem2->e[0] * stem2->a[0],
    stem2->e[0] * stem2->a[1]
  };
  double  stem2_eb[2] = {
    stem2->e[1] * stem2->b[0],
    stem2->e[1] * stem2->b[1]
  };
  double  B2[2] = {
    stem2->c[0] + stem2_ea[0] + stem2_eb[0],
    stem2->c[1] + stem2_ea[1] + stem2_eb[1]
  };
  double  C2[2] = {
    stem2->c[0] + stem2_ea[0] - stem2_eb[0],
    stem2->c[1] + stem2_ea[1] - stem2_eb[1]
  };
  double  D2[2] = {
    stem2->c[0] - stem2_ea[0] - stem2_eb[0],
    stem2->c[1] - stem2_ea[1] - stem2_eb[1]
  };
  double  A2[2] = {
    stem2->c[0] - stem2_ea[0] + stem2_eb[0],
    stem2->c[1] - stem2_ea[1] + stem2_eb[1]
  };

  /*
   * Only the sides of the stems (AB, CD) need to be intersected against
   * each other
   */
  if (intersectLineSegments(A1, B1, A2, B2, NULL)
      || intersectLineSegments(A1, B1, C2, D2, NULL)
      || intersectLineSegments(C1, D1, A2, B2, NULL)
      || intersectLineSegments(C1, D1, C2, D2, NULL)
      )
    return 1;
  else
    return 0;
}


/* Returns true if circle circ intersects triangle ABC, false otherwise. */
PRIVATE short
TestCircleTriangle(const double circ_c[2],
                   const double circ_r,
                   const double A[2],
                   const double B[2],
                   const double C[2],
                   double       p[2])
{
  /* Find point P on triangle ABC closest to circle center */
  ClosestPtPointBulge(circ_c, A, B, C, p);

  /*
   * circle and triangle intersect if the (squared) distance from circle
   * center to point is less than the (squared) circle radius
   */
  double  v[2];
  vector(p, circ_c, v);

  short   ret = scalarProduct2D(v, v) <= (circ_r * circ_r);

  return ret;
}


PRIVATE short
TestLoopBulge(const double  c[2],
              const double  r,
              const double  pPrev[2],
              const double  pThis[2],
              const double  pNext[2])
{
  /* check if bulge point is inside loop */
  double vCenterBulge[2];

  vector(c, pThis, vCenterBulge);
  /* compare r and length of vCenterBulge; using squared distances is less expensive that sqrt */
  if (r * r > vectorLength2DSquared(vCenterBulge))

    return 1;

  double  vPrevThis[2];
  vector(pPrev, pThis, vPrevThis);
  double  vThisNext[2];
  vector(pThis, pNext, vThisNext);

  double  cut1[2], cut2[2];
  short   numCut;

  /* evaluate line pPrev->pThis */
  numCut = getCutPointsOfCircleAndLine(c, r, pPrev, vPrevThis, cut1, cut2);

  if (numCut > 0) {
    if (matchLinePoint(pPrev, vPrevThis, cut1))

      return 1;
  }

  if (numCut > 1)
    if (matchLinePoint(pPrev, vPrevThis, cut2))
      return 1;

  /* evaluate line pThis->pNext */
  numCut = getCutPointsOfCircleAndLine(c, r, pThis, vThisNext, cut1, cut2);
  if (numCut > 0)
    if (matchLinePoint(pThis, vThisNext, cut1))
      return 1;

  if (numCut > 1) {
    if (matchLinePoint(pThis, vThisNext, cut2))

      return 1;
  }

  return 0;
}


PRIVATE short
intersectLoopBulges(const loopBox *loop,
                    const stemBox *stem,
                    int           *bulge)
{
  *bulge = -1;

  double  c[2] = {
    loop->c[0], loop->c[1]
  };
  double  r = loop->r + epsilonRecognize;

  for (int currentBulge = 0; currentBulge < stem->bulgeCount; currentBulge++) {
    double  A[2], B[2], C[2];
    getBulgeCoordinates(stem, currentBulge, A, B, C);

    double  p[2];

    if (TestCircleTriangle(c, r, A, B, C, p)) {
      *bulge = currentBulge;
      return 1;
    }
  }

  return 0;
}


PRIVATE short
intersectBulgesBulges(const stemBox *stem1,
                      const stemBox *stem2,
                      int           *bulge1,
                      int           *bulge2)
{
  *bulge1 = -1;
  *bulge2 = -1;

  double distance = 0.5 * epsilonRecognize;

  for (int currentBulge1 = 0; currentBulge1 < stem1->bulgeCount; currentBulge1++) {
    double piPrev[2], piThis[2], piNext[2];
    getBulgeCoordinatesExtraDistance(stem1, currentBulge1, distance, piPrev, piThis, piNext);

    for (int currentBulge2 = 0; currentBulge2 < stem2->bulgeCount; currentBulge2++) {
      double pjPrev[2], pjThis[2], pjNext[2];
      getBulgeCoordinatesExtraDistance(stem2, currentBulge2, distance, pjPrev, pjThis, pjNext);

      if (intersectLineSegments(piPrev, piThis, pjPrev, pjThis, NULL)
          || intersectLineSegments(piPrev, piThis, pjThis, pjNext, NULL)
          || intersectLineSegments(piThis, piNext, pjPrev, pjThis, NULL)
          || intersectLineSegments(piThis, piNext, pjThis, pjNext, NULL)
          ) {
        *bulge1 = currentBulge1;
        *bulge2 = currentBulge2;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE short
intersectStemBulges(const stemBox *stem1,
                    const stemBox *stem2,
                    int           *bulge2)
{
  *bulge2 = -1;

  if (stem2->bulgeCount == 0)

    return 0;

  /*
   * simplify to only check bulge lines against left and right stem lines
   *
   * if the bulge is surrounded by the stem then there is a Stem vs. Stem intersection
   * if the bulge intersects the stem's bottom or top line then there is an intersection with the adjacent loop
   * -> no need to checks those cases
   */

  /*
   * N - North, E - East, S - South, W - West
   * north is direction to loop
   */
  double  pNW[2];
  pNW[0]  = stem1->c[0] + stem1->e[0] * stem1->a[0] - stem1->e[1] * stem1->b[0];
  pNW[1]  = stem1->c[1] + stem1->e[0] * stem1->a[1] - stem1->e[1] * stem1->b[1];

  double  pSW[2];
  pSW[0]  = stem1->c[0] - stem1->e[0] * stem1->a[0] - stem1->e[1] * stem1->b[0];
  pSW[1]  = stem1->c[1] - stem1->e[0] * stem1->a[1] - stem1->e[1] * stem1->b[1];

  double  pNE[2];
  pNE[0]  = stem1->c[0] + stem1->e[0] * stem1->a[0] + stem1->e[1] * stem1->b[0];
  pNE[1]  = stem1->c[1] + stem1->e[0] * stem1->a[1] + stem1->e[1] * stem1->b[1];

  double  pSE[2];
  pSE[0]  = stem1->c[0] - stem1->e[0] * stem1->a[0] + stem1->e[1] * stem1->b[0];
  pSE[1]  = stem1->c[1] - stem1->e[0] * stem1->a[1] + stem1->e[1] * stem1->b[1];

  double  distance = epsilonRecognize;

  for (int currentBulge2 = 0; currentBulge2 < stem2->bulgeCount; currentBulge2++) {
    double pPrev[2], pThis[2], pNext[2];
    getBulgeCoordinatesExtraDistance(stem2, currentBulge2, distance, pPrev, pThis, pNext);

    if (intersectLineSegments(pNW, pSW, pPrev, pThis, NULL)
        || intersectLineSegments(pNW, pSW, pThis, pNext, NULL)
        || intersectLineSegments(pNE, pSE, pPrev, pThis, NULL)
        || intersectLineSegments(pNE, pSE, pThis, pNext, NULL)
        ) {
      *bulge2 = currentBulge2;
      return 1;
    }
  }

  return 0;
}


#endif
