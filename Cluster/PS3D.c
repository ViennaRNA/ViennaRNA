#include <stdio.h>
#include <math.h>
#include "utils.h"

#define PUBLIC
#define PRIVATE    static

PUBLIC void ps3d_Preambel(FILE *fp, float view[3], float axis[3], char *projtype);
PUBLIC void PS_DrawSimplifiedBox(float X,float Y, float Z, float T, float P, FILE *fp);


PUBLIC void ps3d_Preambel(FILE *fp, float view[3], float axis[3], char *projtype)
{
   fprintf(fp,"%%!PS\n");
   fprintf(fp,"%%%%Title: RNA DotPlot\n");
   fprintf(fp,"%%%%Creator: RNAfold V.01.002c  - by Oymolon\n");
   fprintf(fp,"%%%%CreationDate: %s", time_stamp());
/*   fprintf(fp,"%%%%BoundingBox: 66 211 518 662\n"); */
   fprintf(fp,"%%%%Pages: 1\n");
   fprintf(fp,"%%%%EndComments: Geil eyh !?\n");
   fprintf(fp,"%%Viewing vector\n");
   fprintf(fp,"/v1 { %f } def\n", view[0]);
   fprintf(fp,"/v2 { %f } def\n", view[1]);
   fprintf(fp,"/v3 { %f } def\n", view[2]);
   fprintf(fp,"/a1 { %f } def\n", axis[0]);
   fprintf(fp,"/a2 { %f } def\n", axis[1]);
   fprintf(fp,"/a3 { %f } def\n", axis[2]);
   fprintf(fp,"%%Define some coefficients for the projections.\n");
   fprintf(fp,"/y1a { v2 a3 mul v3 a2 mul sub } def\n");
   fprintf(fp,"/y2a { v3 a1 mul v1 a3 mul sub } def\n");
   fprintf(fp,"/y3a { v1 a2 mul v2 a1 mul sub } def\n");
   fprintf(fp,"/ya_len { y1a y1a mul y2a y2a mul y3a y3a mul add add sqrt } def\n");
   fprintf(fp,"/x1a { v2 y3a mul v3 y2a mul sub } def\n");
   fprintf(fp,"/x2a { v3 y1a mul v1 y3a mul sub } def\n");
   fprintf(fp,"/x3a { v1 y2a mul v2 y1a mul sub } def\n");
   fprintf(fp,"/xa_len { x1a x1a mul x2a x2a mul x3a x3a mul add add sqrt } def\n");
   fprintf(fp,"/x1 { x1a xa_len div } def\n");
   fprintf(fp,"/x2 { x2a xa_len div } def\n");
   fprintf(fp,"/x3 { x3a xa_len div } def\n");
   fprintf(fp,"/y1 { y1a ya_len div } def\n");
   fprintf(fp,"/y2 { y2a ya_len div } def\n");
   fprintf(fp,"/y3 { y3a ya_len div } def\n");
   fprintf(fp,"/sx { v1 v3 div } def\n");
   fprintf(fp,"/sy { v2 v3 div } def\n");
   fprintf(fp,"/v_len { v1 v1 mul v2 v2 mul v3 v3 mul add add sqrt } def\n");
   fprintf(fp,"/u_len { v1 v1 mul v2 v2 mul add sqrt } def\n");
   fprintf(fp,"/u1 { v2     u_len div } def\n");
   fprintf(fp,"/u2 { v1 neg u_len div } def\n");
   fprintf(fp,"/w1 { v1 v3 mul neg u_len div v_len div } def\n");
   fprintf(fp,"/w2 { v2 v3 mul neg u_len div v_len div } def\n");
   fprintf(fp,"/w3 { v1 v1 mul v2 v2 mul add u_len div v_len div } def\n");
   fprintf(fp,"%%Projection Operators\n");
   fprintf(fp,"%%  Projection onto x-y plane in direction v\n");
   fprintf(fp,"/ProjXY { dup sx mul 4 -1 roll exch sub \n");
   fprintf(fp,"          3 1 roll sy mul sub \n");
   fprintf(fp,"    } def \n");
   fprintf(fp,"%%  Normal Projection onto a plane normal to (v1,v2,v3!=0)\n");
   fprintf(fp,"/Projnn  { w3 mul 3 1 roll dup w2 mul 3 1 roll u2 mul exch \n");
   fprintf(fp,"          dup w1 mul 3 1 roll u1 mul add 4 1 roll add add \n");
   fprintf(fp,"    } def \n");
   fprintf(fp,"%%  General Normal Projection\n");
   fprintf(fp,"/ProjN  { dup  y3 mul 4 1 roll x3 mul  3 1 roll dup  \n");
   fprintf(fp,"          y2 mul  5 1 roll x2 mul 3 1 roll dup  y1 mul \n");
   fprintf(fp,"          6 1 roll x1 mul add add 4 1 roll add add  \n");
   fprintf(fp,"          } def\n");
   fprintf(fp,"/Proj { Proj%s } def\n", projtype);
   fprintf(fp,"/L3   { Proj lineto  } def\n");
   fprintf(fp,"/RL3  { Proj rlineto } def\n");
   fprintf(fp,"/M3   { Proj moveto  } def\n");
   fprintf(fp,"/RM3  { Proj rmoveto } def\n");
   fprintf(fp,"%% end 3D macros\n\n");
}

PUBLIC void PS_DrawSimplifiedBox(float X,float Y, float Z, float T, float P, FILE *fp)
{
   float t1, p1, t, lw;
   t1 = T/sqrt(2.);
   p1 = P/sqrt(3.);
   t = X+t1+2*p1;    /* approximate size */
   lw = 0.001*t;
   fprintf(fp,"/scaling_factor { %f } def \n", 400./t);
   fprintf(fp,"%f %f %f Proj \n", X+t1, Y+t1, Z+t1); 
   fprintf(fp,"scaling_factor  mul 2 div 300 sub neg exch \n" );
   fprintf(fp,"scaling_factor  mul 2 div 400 sub neg translate \n");
   fprintf(fp,"scaling_factor dup scale\n");
   fprintf(fp,"newpath\n");
   fprintf(fp,"%f setlinewidth\n",lw);
   fprintf(fp,"0  0  0   M3\n");
   fprintf(fp,"%f %f %f  L3\n", X    , 0.   , 0.   );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , t1   , 0.   );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , Y+t1 , 0.   );
   fprintf(fp,"%f %f %f  L3\n", t1   , Y+t1 , 0.   );
   fprintf(fp,"%f %f %f  L3\n", 0.   , Y    , 0.   );
   fprintf(fp,"%f %f %f  L3\n", 0.   , 0.   , 0.   );
   fprintf(fp,"%f %f %f  L3\n", 0.   , 0.   , Z    );
   fprintf(fp,"%f %f %f  L3\n", 0.   , t1   , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", 0.   , Y+t1 , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", 0.   , Y+t1 , t1   );
   fprintf(fp,"%f %f %f  L3\n", 0.   , Y    , 0.   );
   fprintf(fp,"%f %f %f  L3\n", 0.   , 0.   , 0.   );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n", t1   , Y+t1 , 0.   );
   fprintf(fp,"%f %f %f  L3\n", 0.   , Y+t1 , t1   );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n", X+t1 , t1   , 0.   );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , 0.   , t1   );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , 0.   , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , Y    , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , Y+t1 , Z    );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , Y+t1 , 0.   );
   fprintf(fp,"%f %f %f  L3\n", X+t1 , t1   , 0.   );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n", X+t1 , 0.   , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", t1   , 0.   , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", 0.   , t1   , Z+t1 );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n", X+t1 , Y    , Z+t1 ); 
   fprintf(fp,"%f %f %f  L3\n", X    , Y+t1 , Z+t1 );
   fprintf(fp,"%f %f %f  L3\n", 0.   , Y+t1 , Z+t1 );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  X    , 0.   , 0.   );
   fprintf(fp,"%f %f %f  L3\n",  X+t1 , 0.   , t1   );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  0.   , 0.   , Z    );
   fprintf(fp,"%f %f %f  L3\n",  t1   , 0.   , Z+t1 );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  X+t1 , Y+t1 , Z    );
   fprintf(fp,"%f %f %f  L3\n",  X    , Y+t1 , Z+t1 );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  0.   , 0.   , 0.   );
   fprintf(fp,"%f %f %f  RL3\n", -p1  , -p1  , -p1  );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  X+t1 , Y+t1 , 0.   );
   fprintf(fp,"%f %f %f  RL3\n",  p1  ,  p1  , -p1  );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  X+t1 , 0.   , Z+t1 );
   fprintf(fp,"%f %f %f  RL3\n",  p1  , -p1  ,  p1  );
   fprintf(fp,"stroke\n");
   fprintf(fp,"%f %f %f  M3\n",  0.  , Y+t1  , Z+t1 );
   fprintf(fp,"%f %f %f  RL3\n", -p1  ,  p1  ,  p1  );
   fprintf(fp,"stroke\n");
   fprintf(fp,"showpage\n");
}
