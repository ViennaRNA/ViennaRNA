#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#ifdef __cplusplus
}
#endif

#include <stdlib.h>
#include <g2.h>

#ifdef DO_PS
#include <g2_PS.h>
#endif /* DO_PS */

#ifdef DO_FIG
#include <g2_FIG.h>
#endif /* DO_FIG */

#ifdef DO_X11
#include <g2_X11.h>
#endif /* DO_X11 */

#ifdef DO_GD
#include <g2_gd.h>
#endif /* DO_GD */

#ifdef DO_WIN32
#include <g2_Win32.h>
#endif /* DO_WIN32 */


static int
not_here(s)
char *s;
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(name, arg)
char *name;
int arg;
{
    errno = 0;
    switch (*name) {
    case 'A':
	break;
    case 'B':
	break;
    case 'C':
	break;
    case 'D':
	break;
    case 'E':
	break;
    case 'F':
	break;
    case 'G':
	if (strEQ(name, "G2LD"))
#ifdef G2LD
	    return G2LD;
#else
	    goto not_there;
#endif
	if (strEQ(name, "G2_H"))
#ifdef G2_H
	    return G2_H;
#else
	    goto not_there;
#endif
	if (strEQ(name, "G2_VERSION"))
#ifdef G2_VERSION
	    return atof(G2_VERSION);
#else
	    goto not_there;
#endif
	break;
    case 'H':
	break;
    case 'I':
	break;
    case 'J':
	break;
    case 'K':
	break;
    case 'L':
	break;
    case 'M':
	break;
    case 'N':
	break;
    case 'O':
	break;
    case 'P':
	break;
    case 'Q':
	break;
    case 'R':
	break;
    case 'S':
	break;
    case 'T':
	break;
    case 'U':
	break;
    case 'V':
	break;
    case 'W':
	break;
    case 'X':
	break;
    case 'Y':
	break;
    case 'Z':
	break;
    }
    errno = EINVAL;
    return 0;

not_there:
    errno = ENOENT;
    return 0;
}

typedef int* DevType;
typedef	DevType	G2__Device;


MODULE = G2		PACKAGE = G2			


double
constant(name,arg)
	char *		name
	int		arg

MODULE = G2		PACKAGE = G2::Device		PREFIX = g2_

#ifdef DO_X11

G2::Device
g2_newX11(packname="G2::Device", width=100,height=100)
	char * packname
	int width
	int height
        PROTOTYPE: $;$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_X11(width, height);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL

#endif /* DO_X11 */
	

#ifdef DO_WIN32

G2::Device
g2_newWin32(packname="G2::Device", width=100,height=100,filename="Win32 window",type=0)
	char * packname
	int width
	int height
	char * filename	
        int type
        PROTOTYPE: $;$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_win32(width, height, filename, type);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL

#endif /* DO_WIN32 */

#ifdef DO_GD

G2::Device
g2_newGD(packname="G2::Device", filename="g2.png", width=100, height=100, type=1)
	char * packname
	char * filename
	int width
	int height
	int type
        PROTOTYPE: $;$$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_gd(filename, width, height, type);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL

#endif /* DO_GD */

#ifdef DO_PS

G2::Device
g2_newPS(packname="G2::Device", filename="g2.ps", paper=1,orientation=1)
	char * packname
	char * filename
	int paper
	int orientation
        PROTOTYPE: $;$$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_PS(filename, paper, orientation);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL
	
G2::Device
g2_newEPSF(packname="G2::Device", filename="g2.eps")
	char * packname
	char * filename
        PROTOTYPE: $;$$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_EPSF(filename);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL
	
G2::Device
g2_newEPSF_CLIP(packname="G2::Device", filename="g2.eps",width=100,height=100)
	char * packname
	char * filename
	long width
	long height
        PROTOTYPE: $;$$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_EPSF_CLIP(filename,width,height);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL
	
#endif /* DO_PS */

#ifdef DO_FIG

G2::Device
g2_newFIG(packname="G2::Device", filename="g2.fig")
	char * packname
	char * filename
        PROTOTYPE: $;$$$
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_FIG(filename);
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL
	
#endif /* DO_FIG */

	
G2::Device
g2_newvd(packname="G2::Device")
	char * packname
        PROTOTYPE: $;
	CODE:
	{
		DevType theDevice;
		theDevice = (int *)malloc(sizeof(int));
		*theDevice = g2_open_vd();
		RETVAL = theDevice;
	}
	OUTPUT:
	RETVAL


void
g2_DESTROY(dev)
	G2::Device	dev
	PROTOTYPE: $
	CODE:
	{
		if(g2_device_exist(*dev)) {
			g2_close(*dev);
		}
		free(dev);
	}


void
g2_attach(vd_dev, dev)
	G2::Device	vd_dev
	G2::Device	dev
	PROTOTYPE: $
	CODE:
	{
	g2_attach(*vd_dev, *dev);
	}

void
g2_detach(vd_dev, dev)
	G2::Device 	vd_dev
	G2::Device	dev
	PROTOTYPE: $
	CODE:
	{
	g2_detach(*vd_dev, *dev);
	}

void
g2_close(dev)
	G2::Device	dev
	PROTOTYPE: 
	CODE:
	{
	g2_close(*dev);
	}

void
g2_set_auto_flush(dev, on_off)
	G2::Device	dev
	int	on_off
	PROTOTYPE: $
	CODE:
	{
		g2_set_auto_flush(*dev, on_off);
	}

void
g2_flush(dev)
	G2::Device	dev
	PROTOTYPE: 
	CODE:
	{
		g2_flush(*dev);
	}

void
g2_save(dev)
	G2::Device	dev
	PROTOTYPE: 
	CODE:
	{
		g2_save(*dev);
	}

void
g2_set_coordinate_system(dev, x_origin, y_origin, x_mul, y_mul)
	G2::Device	dev
	double	x_origin
	double	y_origin
	double	x_mul
	double	y_mul
	PROTOTYPE: $$$$
	CODE:
	{
		g2_set_coordinate_system(*dev, x_origin, y_origin, x_mul, y_mul);
	}


int
g2_ink(pd_dev, red, green, blue)
	G2::Device	pd_dev
	double	red
	double	green
	double	blue
	PROTOTYPE: $$$
	CODE:
	{
		RETVAL = g2_ink(*pd_dev, red, green, blue);
	}
	OUTPUT:
	RETVAL

void
g2_pen(dev, color)
	G2::Device	dev
	int	color
	PROTOTYPE: $
	CODE:
	{
		g2_pen(*dev, color);
	}


void
g2_set_dash(dev, N, dashes=NULL)
	G2::Device	dev
	int	N
	double *	dashes
	PROTOTYPE: $
	CODE:
	{
		g2_set_dash(*dev, N, dashes);
		free(dashes);
	}

void
g2_set_font_size(dev, size)
	G2::Device	dev
	double	size
	PROTOTYPE: $
	CODE:
	{
		g2_set_font_size(*dev, size);
	}

void
g2_set_line_width(dev, w)
	G2::Device	dev
	double	w
	PROTOTYPE: $
	CODE:
	{
		g2_set_line_width(*dev, w);
	}

void
g2_clear_palette(dev)
	G2::Device	dev
	PROTOTYPE: 
	CODE:
	{
		g2_clear_palette(*dev);
	}

void
g2_reset_palette(dev)
	G2::Device	dev
	PROTOTYPE: 
	CODE:
	{
		g2_reset_palette(*dev);
	}

void
g2_allocate_basic_colors(dev)
	G2::Device	dev
	PROTOTYPE: 
	CODE:
	{
		g2_allocate_basic_colors(*dev);
	}

void
g2_clear(dev)
	G2::Device	dev
	PROTOTYPE:
	CODE:
	{
		g2_clear(*dev);
	}

void
g2_set_background(dev, color)
	G2::Device	dev
	int	color
	PROTOTYPE: $
	CODE:
	{
		g2_set_background(*dev, color);
	}

void
g2_move(dev, x, y)
	G2::Device	dev
	double	x
	double	y
	PROTOTYPE: $$
	CODE:
	{
		g2_move(*dev, x, y);
	}

void
g2_move_r(dev, dx, dy)
	G2::Device	dev
	double	dx
	double	dy
	PROTOTYPE: $$
	CODE:
	{
		g2_move_r(*dev, dx, dy);
	}

void
g2_plot(dev, x, y)
	G2::Device	dev
	double	x
	double	y
	PROTOTYPE: $$
	CODE:
	{
		g2_plot(*dev, x, y);
	}

void
g2_plot_r(dev, dx, dy)
	G2::Device	dev
	double	dx
	double	dy
	PROTOTYPE: $$
	CODE:
	{
		g2_plot_r(*dev, dx, dy);
	}



void
g2_line(dev, x1, y1, x2, y2)
	G2::Device	dev
	double	x1
	double	y1
	double	x2
	double	y2

	PROTOTYPE: $$$$
	CODE:
	{
		g2_line(*dev, x1, y1, x2, y2);
	}

void
g2_line_r(dev, dx, dy)
	G2::Device	dev
	double	dx
	double	dy
	PROTOTYPE: $$
	CODE:
	{
		g2_line_r(*dev, dx, dy);
	}

void
g2_line_to(dev, x, y)
	G2::Device	dev
	double	x
	double	y
	PROTOTYPE: $$
	CODE:
	{
		g2_line_to(*dev, x, y);
	}

void
g2_poly_line(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_poly_line(*dev, N_pt, points);
		free(points);
	}

void
g2_triangle(dev, x1, y1, x2, y2, x3, y3)
	G2::Device	dev
	double	x1
	double	y1
	double	x2
	double	y2
	double	x3
	double	y3
	PROTOTYPE: $$$$$$$
	CODE:
	{
		g2_triangle(*dev, x1, y1, x2, y2, x3, y3);
	}


void
g2_filled_triangle(dev, x1, y1, x2, y2, x3, y3)
	G2::Device	dev
	double	x1
	double	y1
	double	x2
	double	y2
	double	x3
	double	y3
	PROTOTYPE: $$$$$$$
	CODE:
	{
		g2_filled_triangle(*dev, x1, y1, x2, y2, x3, y3);
	}

void
g2_rectangle(dev, x1, y1, x2, y2)
	G2::Device	dev
	double	x1
	double	y1
	double	x2
	double	y2
	PROTOTYPE: $$$$
	CODE:
	{
		g2_rectangle(*dev, x1, y1, x2, y2);
	}


void
g2_filled_rectangle(dev, x1, y1, x2, y2)
	G2::Device	dev
	double	x1
	double	y1
	double	x2
	double	y2
	PROTOTYPE: $$$$
	CODE:
	{
		g2_filled_rectangle(*dev, x1, y1, x2, y2);
	}

void
g2_polygon(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_polygon(*dev, N_pt, points);
		free(points);
	}

void
g2_filled_polygon(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_filled_polygon(*dev, N_pt, points);
		free(points);
	}

void
g2_circle(dev, x, y, r)
	G2::Device	dev
	double	x
	double	y
	double	r
	PROTOTYPE: $$$
	CODE:
	{
		g2_circle(*dev, x, y, r);
	}

void
g2_filled_circle(dev, x, y, r)
	G2::Device	dev
	double	x
	double	y
	double	r
	PROTOTYPE: $$$
	CODE:
	{
		g2_filled_circle(*dev, x, y, r);
	}

void
g2_ellipse(dev, x, y, r1, r2)
	G2::Device	dev
	double	x
	double	y
	double	r1
	double	r2
	PROTOTYPE: $$$$
	CODE:
	{
		g2_ellipse(*dev, x, y, r1, r2);
	}

void
g2_filled_ellipse(dev, x, y, r1, r2)
	G2::Device	dev
	double	x
	double	y
	double	r1
	double	r2
	PROTOTYPE: $$$$
	CODE:
	{
		g2_filled_ellipse(*dev, x, y, r1, r2);
	}

void
g2_arc(dev, x, y, r1, r2, a1, a2)
	G2::Device	dev
	double	x
	double	y
	double	r1
	double	r2
	double	a1
	double	a2
	PROTOTYPE: $$$$$$
	CODE:
	{
		g2_arc(*dev, x, y, r1, r2, a1, a2);
	}

void
g2_filled_arc(dev, x, y, r1, r2, a1, a2)
	G2::Device	dev
	double	x
	double	y
	double	r1
	double	r2
	double	a1
	double	a2
	PROTOTYPE: $$$$$$
	CODE:
	{
		g2_filled_arc(*dev, x, y, r1, r2, a1, a2);
	}

void
g2_string(dev, x, y, text)
	G2::Device	dev
	double	x
	double	y
	char *	text
	PROTOTYPE: $$$
	CODE:
	{
		g2_string(*dev, x, y, text);
	}

void
g2_set_QP(dev, d, shape)
	G2::Device	dev
	double	d
	enum QPshape	shape
	PROTOTYPE: $$
	CODE:
	{
		g2_set_QP(*dev, d, shape);
	}

void
g2_plot_QP(dev, x, y)
	G2::Device	dev
	double	x
	double	y
	PROTOTYPE: $$
	CODE:
	{
		g2_plot_QP(*dev, x, y);
	}

void
g2_query_pointer(dev)
	G2::Device	dev
	PROTOTYPE: $$
	CODE:
	{
		double x, y;
		unsigned int button;
		g2_query_pointer(*dev, &x, &y, &button);
	}





void
g2_spline(dev, N_pt, points, o)
	G2::Device	dev
	int	N_pt
	double * points
	int o
	PROTOTYPE: $$
	CODE:
	{
		g2_spline(*dev, N_pt, points, o);
		free(points);
	}


void
g2_b_spline(dev, N_pt, points, o)
	G2::Device	dev
	int	N_pt
	double * points
	int o
	PROTOTYPE: $$
	CODE:
	{
		g2_b_spline(*dev, N_pt, points, o);
		free(points);
	}


void
g2_raspln(dev, N_pt, points, tn)
	G2::Device	dev
	int	N_pt
	double * points
	double tn
	PROTOTYPE: $$
	CODE:
	{
		g2_raspln(*dev, N_pt, points, tn);
		free(points);
	}


void
g2_para_3(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_para_3(*dev, N_pt, points);
		free(points);
	}


void
g2_para_5(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_para_5(*dev, N_pt, points);
		free(points);
	}


void
g2_filled_spline(dev, N_pt, points, o)
	G2::Device	dev
	int	N_pt
	double * points
	int o
	PROTOTYPE: $$
	CODE:
	{
		g2_filled_spline(*dev, N_pt, points, o);
		free(points);
	}


void
g2_filled_b_spline(dev, N_pt, points, o)
	G2::Device	dev
	int	N_pt
	double * points
	int o
	PROTOTYPE: $$
	CODE:
	{
		g2_filled_b_spline(*dev, N_pt, points, o);
		free(points);
	}


void
g2_filled_raspln(dev, N_pt, points, tn)
	G2::Device	dev
	int	N_pt
	double * points
	double tn
	PROTOTYPE: $$
	CODE:
	{
		g2_filled_raspln(*dev, N_pt, points, tn);
		free(points);
	}


void
g2_filled_para_3(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_filled_para_3(*dev, N_pt, points);
		free(points);
	}


void
g2_filled_para_5(dev, N_pt, points)
	G2::Device	dev
	int	N_pt
	double * points
	PROTOTYPE: $$
	CODE:
	{
		g2_filled_para_5(*dev, N_pt, points);
		free(points);
	}
