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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2.h"
#include "g2_device.h"
#include "g2_util.h"
#include "g2_config.h"

#include "g2_gd_P.h"
#include "g2_gd.h"

#define PDP ((struct g2_gd_STRUCT *)pdp)

#include "g2_gd_funix.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */


/**
 * \ingroup physdev
 * \defgroup GD GD
 */


/**
 *
 * Create a GD (bitmap image) device.
 *
 * \param filename output file name
 * \param width width
 * \param height height
 * \param gd_type file type, see ::g2_gd_type
 *
 * \return physical device id
 *
 * \ingroup GD
 */
int  g2_open_gd(const char *filename, int width, int height, enum g2_gd_type gd_type)
{
    int pid=-1;
    g2_gd_STRUCT *pdp;
    
    pdp = (g2_gd_STRUCT *)malloc(sizeof(g2_gd_STRUCT));

    pdp->width = width;
    pdp->height = height;
    pdp->gd_type = gd_type;
    pdp->im = gdImageCreate(width,height);
    pdp->f = fopen(filename,"wb");
    pdp->NoOfInks = 0;
    pdp->BackCol = 0;
	
    pid = g2_register_physical_device(pid, pdp,
				      g2_IntCoor, g2_gd_funix,
				      1.0, -1.0,
				      0.0, height-1);

    g2_gd_Clear(pid,pdp);
    g2_set_line_width(pid, 0.0);
    g2_set_font_size(pid, 12.0);
    g2_allocate_basic_colors(pid);
    g2_pen(pid, 1);
    	
    return pid;
}



int g2_gd_Alloc_Basic(int pid, void *pdp)
	{
	int icol;
	for (icol=0;icol<32;icol++)
		gdImageColorAllocate(PDP->im,g2_Basic_Colors[icol][0]/256,g2_Basic_Colors[icol][1]/256,g2_Basic_Colors[icol][2]/256);
	PDP->NoOfInks = 32;
	return 0;
	}


int g2_gd_Clear(int pid, void *pdp)
	{
    gdImageFilledRectangle(PDP->im, 0, 0, PDP->width, PDP->height, PDP->BackCol);
	return 0;
    }

int g2_gd_Save(int pid, void *pdp)
	{
	if (PDP->gd_type == g2_gd_png)
		gdImagePng(PDP->im,PDP->f);
	else if (PDP->gd_type == g2_gd_jpeg)	
		gdImageJpeg(PDP->im,PDP->f,-1);
#ifdef DO_GIF
	else if (PDP->gd_type == g2_gd_gif)	
		gdImageGif(PDP->im,PDP->f);
#endif 
	fflush(PDP->f);
	rewind(PDP->f);
	return 0;
	}

int g2_gd_Delete(int pid, void *pdp)
	{
	g2_gd_Save(pid,pdp);
	fclose(PDP->f);
	gdImageDestroy(PDP->im);
	free(PDP);
	return 0;
    }

int g2_gd_Flush(int pid, void *pdp)
	{
	return 0;
    }

int g2_gd_Pen(int pid, void *pdp, int color)
	{
	PDP->CurCol = color;
	return 0;
    }

int g2_gd_Ink(int pid, void *pdp, double red, double green, double blue)
	{
	if(PDP->NoOfInks == 256)
		return -1;
	else
		PDP->NoOfInks++;
	return gdImageColorAllocate(PDP->im,(int)(255*red),(int)(255*green),(int)(255*blue));
	}

int g2_gd_ClearPalette(int pid, void *pdp)
	{
	int i;
	for (i=0;i<PDP->NoOfInks;i++)
		gdImageColorDeallocate(PDP->im,i);
	PDP->NoOfInks = 0;
	return 0;
    }

int g2_gd_ResetPalette(int pid, void *pdp)
	{
	g2_gd_ClearPalette(pid,pdp);
	g2_gd_Alloc_Basic(pid,pdp);
    return 0;
	}

int g2_gd_SetBackground(int pid, void *pdp, int color)
	{
	PDP->BackCol = color;
	return 0;
    }

int g2_gd_SetLineWidth(int pid, void *pdp, int w)
	{
	PDP->LineWidth = w;
	return 0;
        }
/*	  {
	  if (PDP->brush != NULL)
	    {
            gdImageDestroy(PDP->brush);
	    }
	  PDP->brush = gdImageCreate(w,w);
	  gdImageColorTransparent(PDP->brush,0);
	  gdImageColorAllocate(PDP->brush,0,0,0);
	  gdImageColorAllocate(PDP->brush,
			 gdImageRed(PDP->im,PDP->CurCol),
			 gdImageGreen(PDP->im,PDP->CurCol),
			 gdImageBlue(PDP->im,PDP->CurCol));
	  gdImageArc(PDP->brush, w/2, w/2, w/2,w/2, 0, 360, 1);
	  gdImageFill(PDP->brush,w/2,w/2,1);
          gdImageSetBrush(PDP->im, PDP->brush);
	  PDP->OldCol = PDP->CurCol;
	  PDP->CurCol = gdBrushed;
	  }
	else
	  {
	  PDP->CurCol = PDP->OldCol;
	  }
	return 0;
    }
*/
     
int g2_gd_SetDash(int pid, void *pdp, int n, char *data)
	{
	return 0;
    }

int g2_gd_SetDashX(int pid, void *pdp, int N, double *dashes)
	{
	return 0;
    }

int g2_gd_SetFontSize(int pid, void *pdp, int size)
	{
	if (size <=10)
		PDP->FontSize = gdFontTiny;
	else if (size <=12) 
		PDP->FontSize = gdFontSmall;
	else if (size <=13) 
		PDP->FontSize = gdFontMediumBold;
	else if (size <=15) 
		PDP->FontSize = gdFontLarge;
	else  
		PDP->FontSize = gdFontGiant;
	return 0;
    }

int g2_gd_Plot(int pid, void *pdp, int x, int y)
	{
	gdImageSetPixel(PDP->im, x, y, PDP->CurCol);
	return 0;
    }

int g2_gd_Line(int pid, void *pdp, int x1, int y1, int x2, int y2)
	{
	if (PDP->LineWidth <= 1)
	  gdImageLine(PDP->im, x1, y1, x2, y2, PDP->CurCol);
	else
	  {
	  float dx,dy,l;
	  gdPoint points[4];
	  dx = -(float)(y2-y1);
	  dy =  (float)(x2-x1);
	  l  =  (float)(PDP->LineWidth/sqrt(dy*dy+dx*dx)/2.);
	  dx =  dx*l;
	  dy =  dy*l;
	  points[0].x = (int)(x1+dx);
      points[0].y = (int)(y1+dy);
      points[1].x = (int)(x1-dx);
      points[1].y = (int)(y1-dy);
      points[2].x = (int)(x2-dx);
      points[2].y = (int)(y2-dy);
      points[3].x = (int)(x2+dx);
      points[3].y = (int)(y2+dy);
	  gdImageFilledPolygon(PDP->im,points,4,PDP->CurCol);
	  }
	return 0;
    }

int g2_gd_PolyLine(int pid, void *pdp, int N, int *points)
	{
	return 0;
    }

int g2_gd_Triangle(int pid, void *pdp, int x1, int y1,
		 int x2, int y2,
		 int x3, int y3)
	{
	return 0;
    }

int g2_gd_FilledTriangle(int pid, void *pdp, int x1, int y1,
		       int x2, int y2,
		       int x3, int y3)
	{
	return 0;
    }

int g2_gd_Rectangle(int pid, void *pdp, int x, int y, int x2, int y2)
	{
    gdImageRectangle(PDP->im, x, y, x2, y2, PDP->CurCol);
	return 0;
    }

int g2_gd_FilledRectangle(int pid, void *pdp, int x, int y, int x2, int y2)
	{
    gdImageFilledRectangle(PDP->im, x, y, x2, y2, PDP->CurCol);
	return 0;
	}

int g2_gd_Polygon(int pid, void *pdp, int N, int *points)
	{
	return 0;
    }

int g2_gd_FilledPolygon(int pid, void *pdp, int N, int *points)
	{
	gdPoint *GIFPolygon;
	int i;
	GIFPolygon = (gdPoint *)malloc(N*sizeof(gdPoint));
	for (i=0;i<N;i++)
		{
		GIFPolygon[i].x = points[2*i];
		GIFPolygon[i].y = points[2*i+1];
		}
	gdImageFilledPolygon(PDP->im,GIFPolygon,N,PDP->CurCol);
	free(GIFPolygon);
	return 0;
    }

int g2_gd_Circle(int pid, void *pdp, int x, int y, int r)
	{
	gdImageArc(PDP->im, (int)x, (int)y, (int)r, (int)r, 0, 360, PDP->CurCol);
	return 0;
    }

int g2_gd_FilledCircle(int pid, void *pdp, int x, int y, int r)
	{
	return 0;
    }

int g2_gd_Ellipse(int pid, void *pdp, int x, int y, int r1, int r2)
	{
	return 0;
    }

int g2_gd_FilledEllipse(int pid, void *pdp, int x, int y, int r1, int r2)
	{
	return 0;
	}

int g2_gd_Arc(int pid, void *pdp, int x, int y, int r1, int r2, double a1, double a2)
	{
	gdImageArc(PDP->im,x,y,2*r1,2*r2,dtoi(a1),dtoi(a2),PDP->CurCol);
	return 0;
	}

int g2_gd_FilledArc(int pid, void *pdp, int x, int y,
		  int r1, int r2,
		  double a1, double a2)
	{
	return 0;
	}

int g2_gd_DrawString(int pid, void *pdp, int x, int y, const char *text)
	{
	
	gdImageString(PDP->im,PDP->FontSize,x,y+2-PDP->FontSize->h,(unsigned char *)text,PDP->CurCol);
	return 0;
	}

