/*****************************************************************************
**  Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
**  This file is part of the g2 library
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/
#ifndef _G2_FUNIX_H
#define _G2_FUNIX_H

/*
 *
 * FUNction IndeX enumeration
 *
 */

#define G2_N_FUNIX  31	  /* WARNING !! */
			  /* number of funix excl. g2_FUNIX_NULL !!! */

typedef enum g2_funix {
    g2_FUNIX_NULL=-1,	  /* null funix (_don't_ _count_ in G2_N_FUNIX !!!!) */
    
    g2_DoNothing=0,	  /* do nothing */
			  /* ... =  */

    g2_Delete,		  /* delete physical device */
			  /* ... =  */

    g2_Ink,		  /* set ink */
			  /* ... = (int)red, (int)green, (int)blue */
			  /* 0 < (red,green,blue) < 65535 (=0xFFFF) */
			  /* ret: color id(>=0), <0 if error */
    
    g2_Pen,		  /* set pen */
			  /* ... = (int)pen */
    
    g2_SetBackground,	  /* set background */
			  /* ... = (int)pen */
    
    g2_ClearPalette,	  /* reset color pallete to default values */
			  /* ... = */
			  /* for list of default colors (0,31) see later */

    g2_SetLineWidth,	  /* set line width (line, rectangle, ...) */
			  /* ... = (double)line width */
    
    g2_SetDash,		  /* set line dash */
			  /* ... = (int)number of descr., (double*) bw list */
    
    g2_SetFontSize,	  /* set font size */
			  /* ... = (double)font size */
    
    g2_Clear,		  /* clear screen(X11), print page(PostScript), ... */
			  /* ... = */
    
    g2_Flush,		  /* update output */
			  /* X11: Flush, PostScript: fflush, ... */

    g2_Save,		  /* save output to device (e.g. file) */
			  /* X11: Flush, etc. */

    g2_Plot,		  /* plot point */
			  /* ... = (double)x, (double)y */

    g2_Line,		  /* draw line */
			  /* ... = (double)x1, (double)y1,
			           (double)x2, (double)y2  */

    g2_PolyLine,	  /* draw poly line */
			  /* ... = (int)N,
			           (double*)dp
			    dp[0]=x1, dp[1]=y1,
			    ...
			    dp[2*N-2]=xN, dp[2*N-1]=yN   */

    g2_Polygon,		  /* draw polygon */
			  /* ... = (int)N,
			           (double*)dp
			    dp[0]=x1, dp[1]=y1,
			    ...
			    dp[2*N-2]=xN, dp[2*N-1]=yN   */

    g2_FilledPolygon,	  /* draw filled polygon */
			  /* ... = (int)N,
			           (double*)dp
			    dp[0]=x1, dp[1]=y1,
			    ...
			    dp[2*N-2]=xN, dp[2*N-1]=yN   */

    g2_Rectangle,	  /* draw rectangle */
			  /* ... = (double)x1, (double)y1,
			           (double)x2, (double)y2  */

    g2_FilledRectangle,	  /* draw filled rectangle */
			  /* ... = (double)x1, (double)x1,
			           (double)x2, (double)y2  */

    g2_Triangle,	  /* draw triangle (should be faster as lines) */
			  /* ... = (double)x1, (double)y1,
			           (double)x2, (double)y2,
				   (double)x3, (double)y3 */
    
    g2_FilledTriangle,	  /* draw filled triangle */
			  /* ... = (double)x1, (double)y1,
			           (double)x2, (double)y2,
				   (double)x3, (double)y3 */
    
    g2_Arc,		  /* draw arc */
			  /* ... = (double)x,  (double)y,
			           (double)r1, (double)r2,
				   (double)a1, (double)a2  */

    g2_FilledArc,	  /* draw filled arc */
			  /* ... = (double)x,  (double)y,
			           (double)r1, (double)r2,
				   (double)a1, (double)a2  */
    
    g2_Ellipse,		  /* draw ellipse */
			  /* ... = (double)x,  (double)y,
			           (double)r1, (double)r2 */

    g2_FilledEllipse,	  /* draw filled ellipse */
			  /* ... = (double)x,  (double)y,
			           (double)r1, (double)r2 */

    g2_Circle,		  /* draw circle */
			  /* ... = (double)x,  (double)y,
			           (double)r */

    g2_FilledCircle,	  /* draw filled circle */
			  /* ... = (double)x,  (double)y,
			           (double)r */

    g2_String,	          /* draw string */
			  /* ... = (double)x, (double)y, (const char*)string */

    g2_Image,		  /* draw (pen) image */
			  /* ... = (double)x, (double)y,
			           (int)x_size, (int)y_size, (int*)pen_array */

    g2_QueryPointer,	  /* query pointer position (mouse) */
			  /* ... =  (double)*x, (double)*y,
			            (unsigned int)*button */
    g2_GetPDHandles       /* get pointers to low level handles */
                          /* ... = void *handles[G2_MAX_NUMBER_PD_HANDLES] */
} g2_funix;


#endif /* _G2_FUNIX_H */
