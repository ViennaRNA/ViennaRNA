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
#include <string.h>
#include "g2.h"
#include "g2_util.h"

/*
 *
 * g2 Fortran Interface
 *
 */
#ifdef LINUX
#define FIF(funame) funame ## __
#else
#define FIF(funame) funame ## _
#endif

#define F_REAL         float    /* everything is float (real) !!!!!!!!!!!!!! */
#define F_CHAR         char     /* only char is char */
#define F_CHAR_LENGTH  int      /* and char length is integer (it is automatically supplied */


/**********************************************************/

#ifdef DO_PS

#include "PS/g2_PS.h"





F_REAL FIF(g2_open_ps)(F_CHAR *text, F_REAL *paper, F_REAL *orientation,
		       F_CHAR_LENGTH length)
{
    char *str;
    int rv;

    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';
    rv=g2_open_PS(str, dtoi(*paper), dtoi(*orientation));
    g2_free(str);
    
    return (F_REAL)rv;
}

F_REAL FIF(g2_open_epsf)(F_CHAR *text, F_CHAR_LENGTH length)
{
    char *str;
    int rv;

    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';
    rv=g2_open_EPSF(str);
    g2_free(str);
    
    return (F_REAL)rv;
}

F_REAL FIF(g2_open_epsf_clip)(F_CHAR *text, F_REAL *width, F_REAL *height,
		       F_CHAR_LENGTH length)
{
    char *str;
    int rv;

    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';
    rv=g2_open_PS(str, dtoi(*width), dtoi(*height));
    g2_free(str);
    
    return (F_REAL)rv;
}


#endif /* DO_PS */

/**********************************************************/

#ifdef DO_X11

#include "X11/g2_X11.h"

F_REAL FIF(g2_open_x11)(F_REAL *width, F_REAL *height)
{
    return (F_REAL)g2_open_X11(*width, *height);
}

F_REAL FIF(g2_open_x11x)( F_REAL *width, F_REAL *height, F_REAL *x, F_REAL *y,
			  F_CHAR *text1, F_CHAR *text2, F_CHAR *text3,
			  F_REAL *icon_width, F_REAL *icon_height,
			  F_CHAR_LENGTH length1, F_CHAR_LENGTH length2, F_CHAR_LENGTH length3)
{
    char *str1,*str2,*str3;
    int rv;

    str1=g2_malloc((length1+1)*sizeof(char));
    str2=g2_malloc((length2+1)*sizeof(char));
    str3=g2_malloc((length3+1)*sizeof(char));
    strncpy(str1, text1, length1);
    strncpy(str2, text2, length2);
    strncpy(str3, text3, length3);
    str1[length1]='\0';
    str2[length2]='\0';
    str3[length3]='\0';
    
    rv=g2_open_X11X(dtoi(*width), dtoi(*height),
                 dtoi(*x), dtoi(*y),
		 str1, str2, str3,
                 dtoi (*icon_width), dtoi(*icon_height));
    g2_free(str1);
    g2_free(str2);
    g2_free(str3);
    
    return (F_REAL)rv;
}

#endif /* DO_X11 */

/**********************************************************/

#ifdef DO_WIN32

#include "Win32/g2_Win32.h"

F_REAL FIF(g2_open_win32)( F_REAL *width, F_REAL *height, F_CHAR *text, F_REAL *type,
		       F_CHAR_LENGTH length)
{
    char *str;
    int rv;

    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';
    rv=g2_open_win32(dtoi(*width), dtoi(*height), str, dtoi (*type));
    g2_free(str);
    
    return (F_REAL)rv;
}

#endif /* DO_WIN32 */

/**********************************************************/

#ifdef DO_GD

#include "GD/g2_gd.h"

F_REAL FIF(g2_open_gd)(F_CHAR *text, F_REAL *width, F_REAL *height, F_REAL *gd_type,
			F_CHAR_LENGTH length)
{
    char *str;
    int rv;
    
    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';

    rv=g2_open_gd(str, *width, *height, *gd_type);
    
    g2_free(str);
    
    return (F_REAL)rv;
}

#endif /* DO_GD */

/**********************************************************/

#ifdef DO_FIG

#include "FIG/g2_FIG.h"

F_REAL FIF(g2_open_fig)(F_CHAR *text,
			F_CHAR_LENGTH length)


{
    char *str;
    int rv;
    
    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';

    rv=g2_open_FIG(str);
    
    g2_free(str);
    
    return (F_REAL)rv;
}

#endif /* DO_FIG */

/**********************************************************/

F_REAL FIF(g2_open_vd)(void)
{
    return (F_REAL)g2_open_vd();
}


void FIF(g2_attach)(F_REAL *vd_dev, F_REAL *dev)
{
    g2_attach(dtoi(*vd_dev), dtoi(*dev));
}


void FIF(g2_detach)(F_REAL *vd_dev, F_REAL *dev)
{
    g2_detach(dtoi(*vd_dev), dtoi(*dev));
}



void FIF(g2_close)(F_REAL *dev)
{
    g2_close(dtoi(*dev));
}


void FIF(g2_set_auto_flush)(F_REAL *dev, F_REAL *on_off)
{
    g2_set_auto_flush(dtoi(*dev), dtoi(*on_off));
}


void FIF(g2_set_coordinate_system)(F_REAL *dev,
				   F_REAL *x_origin, F_REAL *y_origin,
				   F_REAL *x_mul,    F_REAL *y_mul)
{
    g2_set_coordinate_system(dtoi(*dev),
			     *x_origin, *y_origin,
			     *x_mul,    *y_mul);
}


F_REAL FIF(g2_ld)(void)
{
    return (F_REAL)g2_ld();
}


void FIF(g2_set_ld)(F_REAL *dev)
{
    g2_set_ld(dtoi(*dev));
}



void FIF(g2_flush)(F_REAL *dev)
{
    g2_flush(dtoi(*dev));
}


void FIF(g2_save)(F_REAL *dev)
{
    g2_save(dtoi(*dev));
}




void FIF(g2_arc)(F_REAL *dev, F_REAL *x, F_REAL *y,
		 F_REAL *r1, F_REAL *r2, F_REAL *a1, F_REAL *a2)
{
    g2_arc(dtoi(*dev), *x,  *y, *r1,  *r2,  *a1,  *a2);
}


void FIF(g2_circle)(F_REAL *dev, F_REAL *x, F_REAL *y, F_REAL *r)
{
    g2_circle(dtoi(*dev), *x, *y, *r);
}


void FIF(g2_clear)(F_REAL *dev)
{
    g2_clear(dtoi(*dev));
}


void FIF(g2_clear_palette)(F_REAL *dev)
{
    g2_clear_palette(dtoi(*dev));
}


void FIF(g2_string)(F_REAL *dev, F_REAL *x, F_REAL *y, F_CHAR *text,
		    F_CHAR_LENGTH length)
{
    char *str;
    str=g2_malloc((length+1)*sizeof(char));
    strncpy(str, text, length);
    str[length]='\0';
    g2_string(dtoi(*dev), *x, *y, str);
    g2_free(str);
}


void FIF(g2_ellipse)(F_REAL *dev, F_REAL *x, F_REAL *y, F_REAL *r1, F_REAL *r2)
{
    g2_ellipse(dtoi(*dev), *x, *y, *r1, *r2);
}


void FIF(g2_filled_arc)(F_REAL *dev, F_REAL *x, F_REAL *y,
			F_REAL *r1, F_REAL *r2,
			F_REAL *a1, F_REAL *a2)
{
    g2_filled_arc(dtoi(*dev), *x, *y, *r1, *r2, *a1, *a2);
}


void FIF(g2_filled_circle)(F_REAL *dev, F_REAL *x, F_REAL *y, F_REAL *r)
{
    g2_filled_circle(dtoi(*dev), *x, *y, *r);
}


void FIF(g2_filled_ellipse)(F_REAL *dev, F_REAL *x, F_REAL *y, F_REAL *r1, F_REAL *r2)
{
    g2_filled_ellipse(dtoi(*dev), *x, *y, *r1, *r2);
}


void FIF(g2_filled_triangle)(F_REAL *dev, F_REAL *x1, F_REAL *y1,
			     F_REAL *x2, F_REAL *y2,
			     F_REAL *x3, F_REAL *y3)
{
    g2_filled_triangle(dtoi(*dev), *x1, *y1, *x2, *y2, *x3, *y3);
}


F_REAL  FIF(g2_ink)(F_REAL *dev, F_REAL *red, F_REAL *green, F_REAL *blue)
{
    return (F_REAL)g2_ink(dtoi(*dev), *red, *green, *blue);
}


void FIF(g2_line)(F_REAL *dev, F_REAL *x1, F_REAL *y1, F_REAL *x2, F_REAL *y2)
{
    g2_line(dtoi(*dev), *x1, *y1, *x2, *y2);
}


void FIF(g2_poly_line)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_poly_line(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}


void FIF(g2_polygon)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_polygon(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}


void FIF(g2_filled_polygon)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_filled_polygon(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}


void FIF(g2_line_r)(F_REAL *dev, F_REAL *dx, F_REAL *dy)
{
    g2_line_r(dtoi(*dev), *dx, *dy);
}


void FIF(g2_line_to)(F_REAL *dev, F_REAL *x, F_REAL *y)
{
    g2_line_to(dtoi(*dev), *x, *y);
}


void FIF(g2_move)(F_REAL *dev, F_REAL *x, F_REAL *y)
{
    g2_move(dtoi(*dev), *x, *y);
}


void FIF(g2_move_r)(F_REAL *dev, F_REAL *dx, F_REAL *dy)
{
    g2_move_r(dtoi(*dev), *dx, *dy);
}


void FIF(g2_pen)(F_REAL *dev, F_REAL *color)
{
    g2_pen(dtoi(*dev), dtoi(*color));
}


void FIF(g2_plot)(F_REAL *dev, F_REAL *x, F_REAL *y)
{
    g2_plot(dtoi(*dev), *x, *y);
}


void FIF(g2_plot_r)(F_REAL *dev, F_REAL *dx, F_REAL *dy)
{
    g2_plot_r(dtoi(*dev), *dx, *dy);
}


void FIF(g2_rectangle)(F_REAL *dev,
		       F_REAL *x1, F_REAL *y1,
		       F_REAL *x2, F_REAL *y2)
{
    g2_rectangle(dtoi(*dev), *x1, *y1, *x2, *y2);
}


void FIF(g2_filled_rectangle)(F_REAL *dev,
			      F_REAL *x1, F_REAL *y1,
			      F_REAL *x2, F_REAL *y2)
{
    g2_filled_rectangle(dtoi(*dev), *x1, *y1, *x2, *y2);
}


void FIF(g2_reset_palette)(F_REAL *dev)
{
    g2_reset_palette(dtoi(*dev));
}


void FIF(g2_set_background)(F_REAL *dev, F_REAL *color)
{
    g2_set_background(dtoi(*dev), dtoi(*color));
}


void FIF(g2_set_dash)(F_REAL *dev, F_REAL *N, F_REAL *dashes)
{
    double *d;
    d=g2_floatp2doublep(dashes, dtoi(*N));
    g2_set_dash(dtoi(*dev), dtoi(*N), d);
    g2_free(d);
}


void FIF(g2_set_font_size)(F_REAL *dev, F_REAL *size)
{
    g2_set_font_size(dtoi(*dev), *size);
}


void FIF(g2_set_line_width)(F_REAL *dev, F_REAL *w)
{
    g2_set_line_width(dtoi(*dev), *w);
}


void FIF(g2_triangle)(F_REAL *dev, F_REAL *x1, F_REAL *y1,
		      F_REAL *x2, F_REAL *y2,
		      F_REAL *x3, F_REAL *y3)
{
    g2_triangle(dtoi(*dev), *x1, *y1, *x2, *y2, *x3, *y3);
}


void FIF(g2_set_qp)(F_REAL *dev, F_REAL *d, F_REAL *shape)
{
    g2_set_QP(dtoi(*dev), *d, dtoi(*shape));
}


void FIF(g2_plot_qp)(F_REAL *dev, F_REAL *x, F_REAL *y)
{
    g2_plot_QP(dtoi(*dev), *x, *y);
}


/* thanks to Yuri Sbitnev for contributing the g2_image code for FORTRAN */
void FIF(g2_image)(F_REAL *dev, F_REAL *x, F_REAL *y, F_REAL *x_size, F_REAL *y_size,
		   F_REAL *pens)
{
    int i, j, xs, ys;
    int *mypens;
    xs=dtoi(*x_size);
    ys=dtoi(*y_size);
    mypens=(int *) g2_malloc(xs*ys*sizeof(int));
    for(j=0;j<ys;j++) 
      for(i=0;i<xs;i++) 
        mypens[j*xs+i]=dtoi(pens[j*xs+i]);         /* pens[dtoi(*y_size)][dtoi(*x_size)] */
    g2_image(dtoi(*dev), *x, *y, xs, ys, mypens);
    g2_free(mypens);
}






void FIF(g2_spline)(F_REAL *dev, F_REAL *N_pt, F_REAL *points, F_REAL *o)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_spline(dtoi(*dev), dtoi(*N_pt), d, dtoi(*o));
    g2_free(d);
}

void FIF(g2_b_spline)(F_REAL *dev, F_REAL *N_pt, F_REAL *points, F_REAL *o)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_b_spline(dtoi(*dev), dtoi(*N_pt), d, dtoi(*o));
    g2_free(d);
}

void FIF(g2_raspln)(F_REAL *dev, F_REAL *N_pt, F_REAL *points, F_REAL *tn)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_raspln(dtoi(*dev), dtoi(*N_pt), d, *tn);
    g2_free(d);
}

void FIF(g2_para_3)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_para_3(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}

void FIF(g2_para_5)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_para_5(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}

void FIF(g2_filled_spline)(F_REAL *dev, F_REAL *N_pt, F_REAL *points, F_REAL *o)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_filled_spline(dtoi(*dev), dtoi(*N_pt), d, dtoi(*o));
    g2_free(d);
}

void FIF(g2_filled_b_spline)(F_REAL *dev, F_REAL *N_pt, F_REAL *points, F_REAL *o)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_filled_b_spline(dtoi(*dev), dtoi(*N_pt), d, dtoi(*o));
    g2_free(d);
}

void FIF(g2_filled_raspln)(F_REAL *dev, F_REAL *N_pt, F_REAL *points, F_REAL *tn)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_filled_raspln(dtoi(*dev), dtoi(*N_pt), d, *tn);
    g2_free(d);
}

void FIF(g2_filled_para_3)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_filled_para_3(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}

void FIF(g2_filled_para_5)(F_REAL *dev, F_REAL *N_pt, F_REAL *points)
{
    double *d;
    d=g2_floatp2doublep(points, dtoi(*N_pt)*2);
    g2_filled_para_5(dtoi(*dev), dtoi(*N_pt), d);
    g2_free(d);
}

