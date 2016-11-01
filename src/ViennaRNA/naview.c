/*
*   NAVIEW -- A program to make a modified radial drawing of an RNA
*   secondary structure.
*
*   Copyright (c) 1988 Robert E. Bruccoleri
*   Copying of this software, in whole or in part, is permitted
*   provided that the copies are not made for commercial purposes,
*   appropriate credit for the use of the software is given, this
*   copyright notice appears, and notice is given that the copying
*   is by permission of Robert E. Bruccoleri. Any other copying
*   requires specific permission.
*
*   See R. Bruccoleri and G. Heinrich, Computer Applications in the
*   Biosciences, 4, 167-173 (1988) for a full description.
*
*   In November 1997, Michael Zuker made a number of changes to bring
*   naview up to modern standards. All functions defined in naview are
*   now declared before main() with arguments and argument types.
*   When functions are defined, their argument types are declared
*   with the function and these definitions are removed after the '{'.
*   The 'void' declaration was used as necessary for functions.
*
*   The troublesome na_scanf function was deleted and replaced by
*   scanf. Finally, there is now no default for the minimum separation
*   of bases. A floating point number must be entered. However, as
*   before an entry < 0 will be moved up to 0 and an entry > 0.5
*   will be reduced to 0.5.
*
*   Adapted for use as a subroutine in the Vienna RNA Package
*   by Ivo Hofacker, May 1998:
*   deleted output routines, replaced main() by naview_xy_coordinates(),
*   which fills the X and Y arrays used by PS_rna_plot() etc.
*   added ansi prototypes and fixed memory leaks.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/naview.h"

typedef int LOGICAL;
#define logical LOGICAL

#define true 1
#define false 0
#define FATAL_ERROR 1
#define SUCCESS 0

#define type_alloc(type) (type *) vrna_alloc(sizeof(type))

#define struct_alloc(structure_name) type_alloc(struct structure_name)

#define add_double_list(head,tail,newp) {\
	(newp)->next = (newp)->prev = NULL; \
        if ((head) == NULL) (head) = (tail) = (newp); \
	else { \
	     (tail)->next = (newp); \
	     (newp)->prev = (tail); \
	     (tail) = (newp); \
	     } \
	}

static double pi = 3.141592653589793;
static double anum = 9999.0;



/*
*   Function data type definitions
*/

#define minf2(x1, x2) ((x1)<(x2))?(x1):(x2)
#define maxf2(x1, x2) ((x1)>(x2))?(x1):(x2)

static struct base {
  int mate;
  double x,y;
  logical extracted;
  struct region *region;
} *bases;

struct region {
  int start1,end1,start2,end2;
};

struct loop {
  int nconnection;
  struct connection **connections;
  int number;
  int depth;
  logical mark;
  double x,y,radius;
};

struct connection {
  struct loop *loop;
  struct region *region;
  int start,end;       /* Start and end form the 1st base pair of the region. */
  double xrad,yrad,angle;
  logical extruded;	  /* True if segment between this connection and
			     the next must be extruded out of the circle */
  logical broken;	  /* True if the extruded segment must be drawn long. */
};

static int nbase, nregion, loop_count;

static struct loop *root, *loops;

static struct region *regions;

static struct loop *construct_loop(int ibase);

struct radloop {
  double radius;
  int loopnumber;
  struct radloop *next, *prev;
};

static struct radloop *rlphead;

static double lencut;

static logical debug = false;

static void read_in_bases(short *pair_table);
static void find_regions(void);
static void dump_loops(void);
static void find_central_loop(void);
static void determine_depths(void);
static void traverse_loop(struct loop *lp,struct connection *anchor_connection);
static void determine_radius(struct loop *lp,double lencut);
static void generate_region(struct connection *cp);
static void construct_extruded_segment(struct connection *cp,struct connection *cpnext);
static void find_center_for_arc(int n,double b,double *hp,double *thetap);
static int depth(struct loop *lp);

static logical connected_connection(struct connection *cp, struct connection *cpnext);
static int    find_ic_middle(int icstart, int icend, struct connection *anchor_connection, struct connection *acp, struct loop *lp);



int naview_xy_coordinates(short *pair_table, float *X, float *Y) {
  int i;

  nbase = pair_table[0]; /* length */
  bases = (struct base *) vrna_alloc(sizeof(struct base)*(nbase+1));
  regions = (struct region *) vrna_alloc(sizeof(struct region)*(nbase+1));
  read_in_bases(pair_table);
  lencut = 0.5;
  rlphead = NULL;
  find_regions();
  loop_count = 0;
  loops = (struct loop *) vrna_alloc(sizeof(struct loop)*(nbase+1));
  construct_loop(0);
  find_central_loop();
  if (debug) dump_loops();

  traverse_loop(root,NULL);
  for (i=0; i<nbase; i++) {
    X[i] = 100 + 15*bases[i+1].x;
    Y[i] = 100 + 15*bases[i+1].y;
  }
  free(bases);
  free(regions);
  free(loops);
  return nbase;
}


static void read_in_bases(short *pair_table)
{
  int i,npairs;

  /* Set up an origin.  */
  bases[0].mate = 0;
  bases[0].extracted = false;
  bases[0].x = anum;
  bases[0].y = anum;

  for (npairs=0,i=1; i<=nbase; i++) {
    bases[i].extracted = false;
    bases[i].x = anum;
    bases[i].y = anum;
    bases[i].mate = pair_table[i];
    if ((int) pair_table[i]>i) npairs++;
  }
  if (npairs==0) { /* must have at least 1 pair to avoid segfault */
    bases[1].mate=nbase;
    bases[nbase].mate=1;
  }
}


static void find_regions(void)
/*
*   Identifies the regions in the structure.
*/

{
  int i,mate,nb1;
  logical *mark;

  nb1 = nbase + 1;
  mark = (int *) vrna_alloc(sizeof(int)*nb1);
  for (i = 0; i < nb1; i++) mark[i] = false;
  nregion = 0;
  for (i=0; i<=nbase; i++) {
    if ( (mate = bases[i].mate) && !mark[i]) {
      regions[nregion].start1 = i;
      regions[nregion].end2 = mate;
      mark[i] = true;
      mark[mate] = true;
      bases[i].region = bases[mate].region = &regions[nregion];
      for (i++,mate--;
	   i<mate && bases[i].mate == mate;
	   i++,mate--) {
	mark[i] = mark[mate] = true;
	bases[i].region = bases[mate].region = &regions[nregion];
      }
      regions[nregion].end1 = --i;
      regions[nregion].start2 = mate+1;
      if (debug) {
	if (nregion == 0) printf("\nRegions are:\n");
	printf("Region %d is %d-%d and %d-%d with gap of %d.\n",
	       nregion+1,regions[nregion].start1,regions[nregion].end1,
	       regions[nregion].start2,regions[nregion].end2,
	       regions[nregion].start2-regions[nregion].end1+1);
      }
      nregion++;
    }
  }
  free(mark);
}


static struct loop *construct_loop(int ibase)
/*
*   Starting at residue ibase, recursively constructs the loop containing
*   said base and all deeper bases.
*/

{
  int i,mate;
  struct loop *retloop,*lp;
  struct connection *cp;
  struct region *rp;
  struct radloop *rlp;

  retloop = &loops[loop_count++];
  retloop->nconnection = 0;
  retloop->connections = (struct connection **) vrna_alloc(sizeof(struct connection *));
  retloop->depth = 0;
  retloop->number = loop_count;
  retloop->radius = 0.0;
  for (rlp = rlphead;  rlp;  rlp = rlp->next)
    if (rlp->loopnumber == loop_count) retloop->radius = rlp->radius;
  i = ibase;
  do {
    if ((mate = bases[i].mate) != 0) {
      rp = bases[i].region;
      if (!bases[rp->start1].extracted) {
	if (i == rp->start1) {
	  bases[rp->start1].extracted = true;
	  bases[rp->end1].extracted = true;
	  bases[rp->start2].extracted = true;
	  bases[rp->end2].extracted = true;
	  lp = construct_loop(rp->end1 < nbase ? rp->end1+1 : 0);
	}
	else if (i == rp->start2){
	  bases[rp->start2].extracted = true;
	  bases[rp->end2].extracted = true;
	  bases[rp->start1].extracted = true;
	  bases[rp->end1].extracted = true;
	  lp = construct_loop(rp->end2 < nbase ? rp->end2+1 : 0);
	}
	else {
	  vrna_message_error("naview: Error detected in construct_loop. i = %d not found in region table.",i);
	  exit(FATAL_ERROR);
	}
	retloop->connections = (struct connection **)
	  realloc(retloop->connections,
		  (++retloop->nconnection+1) *
		  sizeof(struct connection *));
	retloop->connections[retloop->nconnection-1] = cp =
	  struct_alloc(connection);
	retloop->connections[retloop->nconnection] = NULL;
	cp->loop = lp;
	cp->region = rp;
	if (i == rp->start1) {
	  cp->start = rp->start1;
	  cp->end = rp->end2;
	}
	else {
	  cp->start = rp->start2;
	  cp->end = rp->end1;
	}
	cp->extruded = false;
	cp->broken = false;
	lp->connections = (struct connection **)
	  realloc(lp->connections,
		  (++lp->nconnection+1) *
		  sizeof(struct connection *));
	lp->connections[lp->nconnection-1] = cp =
	  struct_alloc(connection);
	lp->connections[lp->nconnection] = NULL;
	cp->loop = retloop;
	cp->region = rp;
	if (i == rp->start1) {
	  cp->start = rp->start2;
	  cp->end = rp->end1;
	}
	else {
	  cp->start = rp->start1;
	  cp->end = rp->end2;
	}
	cp->extruded = false;
	cp->broken = false;
      }
      i = mate;
    }
    if (++i > nbase) i = 0;
  }
  while (i != ibase);
  return retloop;
}


static void dump_loops(void)
/*
*   Displays all the loops.
*/

{
  int il,ilp,irp;
  struct loop *lp;
  struct connection *cp,**cpp;

  printf("\nRoot loop is #%ld\n",(root-loops)+1);
  for (il=0; il < loop_count; il++) {
    lp = &loops[il];
    printf("Loop %d has %d connections:\n",il+1,lp->nconnection);
    for (cpp = lp->connections; (cp = *cpp); cpp++) {
      ilp = (cp->loop - loops) + 1;
      irp = (cp->region - regions) + 1;
      printf("  Loop %d Region %d (%d-%d)\n",
	     ilp,irp,cp->start,cp->end);
    }
  }
}


static void find_central_loop(void)
/*
*   Find node of greatest branching that is deepest.
*/

{
  struct loop *lp;
  int maxconn,maxdepth,i;

  determine_depths();
  maxconn = 0;
  maxdepth = -1;

  for (i=0; i<loop_count; i++) {
    lp = &loops[i];
    if (lp->nconnection > maxconn) {
      maxdepth = lp->depth;
      maxconn = lp->nconnection;
      root = lp;
    }
    else if (lp->depth > maxdepth && lp->nconnection == maxconn) {
      maxdepth = lp->depth;
      root = lp;
    }
  }
}


static void determine_depths(void)
/*
*   Determine the depth of all loops.
*/

{
  struct loop *lp;
  int i,j;

  for (i=0; i<loop_count; i++) {
    lp = &loops[i];
    for (j=0; j<loop_count; j++) loops[j].mark = false;
    lp->depth = depth(lp);
  }
}



static int depth(struct loop *lp)
/*
*   Determines the depth of loop, lp. Depth is defined as the minimum
*   distance to a leaf loop where a leaf loop is one that has only one
*   or no connections.
*/

{
  struct connection *cp,**cpp;
  int count,ret,d;

  if (lp->nconnection <= 1) return 0;
  if (lp->mark) return -1;
  lp->mark = true;
  count = 0;
  ret = 0;
  for (cpp=lp->connections; (cp = *cpp); cpp++) {
    d = depth(cp->loop);
    if (d >= 0) {
      if (++count == 1) ret = d;
      else if (ret > d) ret = d;
    }
  }
  lp->mark = false;
  return ret+1;
}


static void traverse_loop(struct loop *lp, struct connection *anchor_connection)
/*
*   This is the workhorse of the display program. The algorithm is
*   recursive based on processing individual loops. Each base pairing
*   region is displayed using the direction given by the circle diagram,
*   and the connections between the regions is drawn by equally spaced
*   points. The radius of the loop is set to minimize the square error
*   for lengths between sequential bases in the loops. The "correct"
*   length for base links is 1. If the least squares fitting of the
*   radius results in loops being less than 1/2 unit apart, then that
*   segment is extruded.
*
*   The variable, anchor_connection, gives the connection to the loop
*   processed in an previous level of recursion.
*/

{
  double xs,ys,xe,ye,xn,yn,angleinc,r;
  double radius,xc,yc,xo,yo,astart,aend,a;
  struct connection *cp,*cpnext,**cpp,*acp,*cpprev;
  int i,j,n,ic;
  double da,maxang;
  int count,icstart,icend,icmiddle,icroot;
  logical done,done_all_connections,rooted;
  int sign;
  double midx,midy,nrx,nry,mx,my,vx,vy,dotmv,nmidx,nmidy;
  int icstart1,icup,icdown,icnext,direction;
  double dan,dx,dy,rr;
  double cpx,cpy,cpnextx,cpnexty,cnx,cny,rcn,rc,lnx,lny,rl,ac,acn,sx,sy,dcp;
  int imaxloop;

  angleinc = 2 * pi / (nbase+1);
  acp = NULL;
  icroot = -1;
  for (cpp=lp->connections, ic = 0; (cp = *cpp); cpp++, ic++) {
    /*	xs = cos(angleinc*cp->start);
	ys = sin(angleinc*cp->start);
	xe = cos(angleinc*cp->end);
	ye = sin(angleinc*cp->end); */
    xs = -sin(angleinc*cp->start);
    ys = cos(angleinc*cp->start);
    xe = -sin(angleinc*cp->end);
    ye = cos(angleinc*cp->end);
    xn = ye-ys;
    yn = xs-xe;
    r = sqrt(xn*xn + yn*yn);
    cp->xrad = xn/r;
    cp->yrad = yn/r;
    cp->angle = atan2(yn,xn);
    if (cp->angle < 0.0) cp->angle += 2*pi;
    if (anchor_connection != NULL &&
	anchor_connection->region == cp->region) {
      acp = cp;
      icroot = ic;
    }
  }

 set_radius:
  determine_radius(lp,lencut);
  radius = lp->radius;
  if (anchor_connection == NULL) xc = yc = 0.0;
  else {
    xo = (bases[acp->start].x+bases[acp->end].x) / 2.0;
    yo = (bases[acp->start].y+bases[acp->end].y) / 2.0;
    xc = xo - radius * acp->xrad;
    yc = yo - radius * acp->yrad;
  }

  /*
   *   The construction of the connectors will proceed in blocks of
   *   connected connectors, where a connected connector pairs means
   *   two connectors that are forced out of the drawn circle because they
   *   are too close together in angle.
   */

  /*
   *   First, find the start of a block of connected connectors
   */

  if (icroot == -1)
    icstart = 0;
  else icstart = icroot;
  cp = lp->connections[icstart];
  count = 0;
  if (debug) printf("Now processing loop %d\n",lp->number);
  done = false;
  do {
    j = icstart - 1;
    if (j < 0) j = lp->nconnection - 1;
    cpprev = lp->connections[j];
    if (!connected_connection(cpprev,cp)) {
      done = true;
    }
    else {
      icstart = j;
      cp = cpprev;
    }
    if (++count > lp->nconnection) {
      /*
       *  Here everything is connected. Break on maximum angular separation
       *  between connections.
       */
      maxang = -1.0;
      for (ic = 0;  ic < lp->nconnection;  ic++) {
	j = ic + 1;
	if (j >= lp->nconnection) j = 0;
	cp = lp->connections[ic];
	cpnext = lp->connections[j];
	ac = cpnext->angle - cp->angle;
	if (ac < 0.0) ac += 2*pi;
	if (ac > maxang) {
	  maxang = ac;
	  imaxloop = ic;
	}
      }
      icend = imaxloop;
      icstart = imaxloop + 1;
      if (icstart >= lp->nconnection) icstart = 0;
      cp = lp->connections[icend];
      cp->broken = true;
      done = true;
    }
  } while    (!done);
  done_all_connections = false;
  icstart1 = icstart;
  if (debug) printf("Icstart1 = %d\n",icstart1);
  while (!done_all_connections) {
    count = 0;
    done = false;
    icend = icstart;
    rooted = false;
    while (!done) {
      cp = lp->connections[icend];
      if (icend == icroot) rooted = true;
      j = icend + 1;
      if (j >= lp->nconnection) {
	j = 0;
      }
      cpnext = lp->connections[j];
      if (connected_connection(cp,cpnext)) {
	if (++count >= lp->nconnection) break;
	icend = j;
      }
      else {
	done = true;
      }
    }
    icmiddle = find_ic_middle(icstart,icend,anchor_connection,acp,lp);
    ic = icup = icdown = icmiddle;
    if (debug)
      printf("IC start = %d  middle = %d  end = %d\n",
	     icstart,icmiddle,icend);
    done = false;
    direction = 0;
    while (!done) {
      if (direction < 0) {
	ic = icup;
      }
      else if (direction == 0) {
	ic = icmiddle;
      }
      else {
	ic = icdown;
      }
      if (ic >= 0) {
	cp = lp->connections[ic];
	if (anchor_connection == NULL || acp != cp) {
	  if (direction == 0) {
	    astart = cp->angle - asin(1.0/2.0/radius);
	    aend = cp->angle + asin(1.0/2.0/radius);
	    bases[cp->start].x = xc + radius*cos(astart);
	    bases[cp->start].y = yc + radius*sin(astart);
	    bases[cp->end].x = xc + radius*cos(aend);
	    bases[cp->end].y = yc + radius*sin(aend);
	  }
	  else if (direction < 0) {
	    j = ic + 1;
	    if (j >= lp->nconnection) j = 0;
	    cp = lp->connections[ic];
	    cpnext = lp->connections[j];
	    cpx = cp->xrad;
	    cpy = cp->yrad;
	    ac = (cp->angle + cpnext->angle) / 2.0;
	    if (cp->angle > cpnext->angle) ac -= pi;
	    cnx = cos(ac);
	    cny = sin(ac);
	    lnx = cny;
	    lny = -cnx;
	    da = cpnext->angle - cp->angle;
	    if (da < 0.0) da += 2*pi;
	    if (cp->extruded) {
	      if (da <= pi/2) rl = 2.0;
	      else rl = 1.5;
	    }
	    else {
	      rl = 1.0;
	    }
	    bases[cp->end].x = bases[cpnext->start].x + rl*lnx;
	    bases[cp->end].y = bases[cpnext->start].y + rl*lny;
	    bases[cp->start].x = bases[cp->end].x + cpy;
	    bases[cp->start].y = bases[cp->end].y - cpx;
	  }
	  else {
	    j = ic - 1;
	    if (j < 0) j = lp->nconnection - 1;
	    cp = lp->connections[j];
	    cpnext = lp->connections[ic];
	    cpnextx = cpnext->xrad;
	    cpnexty = cpnext->yrad;
	    ac = (cp->angle + cpnext->angle) / 2.0;
	    if (cp->angle > cpnext->angle) ac -= pi;
	    cnx = cos(ac);
	    cny = sin(ac);
	    lnx = -cny;
	    lny = cnx;
	    da = cpnext->angle - cp->angle;
	    if (da < 0.0) da += 2*pi;
	    if (cp->extruded) {
	      if (da <= pi/2) rl = 2.0;
	      else rl = 1.5;
	    }
	    else {
	      rl = 1.0;
	    }
	    bases[cpnext->start].x = bases[cp->end].x + rl*lnx;
	    bases[cpnext->start].y = bases[cp->end].y + rl*lny;
	    bases[cpnext->end].x = bases[cpnext->start].x - cpnexty;
	    bases[cpnext->end].y = bases[cpnext->start].y + cpnextx;
	  }
	}
      }
      if (direction < 0) {
	if (icdown == icend) {
	  icdown = -1;
	}
	else if (icdown >= 0) {
	  if (++icdown >= lp->nconnection) {
	    icdown = 0;
	  }
	}
	direction = 1;
      }
      else {
	if (icup == icstart) icup = -1;
	else if (icup >= 0) {
	  if (--icup < 0) {
	    icup = lp->nconnection - 1;
	  }
	}
	direction = -1;
      }
      done = icup == -1 && icdown == -1;
    }
    icnext = icend + 1;
    if (icnext >= lp->nconnection) icnext = 0;
    if (icend != icstart && (! (icstart == icstart1 && icnext == icstart1))) {
      /*
       *	    Move the bases just constructed (or the radius) so
       *	    that the bisector of the end points is radius distance
       *	    away from the loop center.
       */
      cp = lp->connections[icstart];
      cpnext = lp->connections[icend];
      dx = bases[cpnext->end].x - bases[cp->start].x;
      dy = bases[cpnext->end].y - bases[cp->start].y;
      midx = bases[cp->start].x + dx/2.0;
      midy = bases[cp->start].y + dy/2.0;
      rr = sqrt(dx*dx + dy*dy);
      mx = dx / rr;
      my = dy / rr;
      vx = xc - midx;
      vy = yc - midy;
      rr = sqrt(dx*dx + dy*dy);
      vx /= rr;
      vy /= rr;
      dotmv = vx*mx + vy*my;
      nrx = dotmv*mx - vx;
      nry = dotmv*my - vy;
      rr = sqrt(nrx*nrx + nry*nry);
      nrx /= rr;
      nry /= rr;
      /*
       *	    Determine which side of the bisector the center should be.
       */
      dx = bases[cp->start].x - xc;
      dy = bases[cp->start].y - yc;
      ac = atan2(dy,dx);
      if (ac < 0.0) ac += 2*pi;
      dx = bases[cpnext->end].x - xc;
      dy = bases[cpnext->end].y - yc;
      acn = atan2(dy,dx);
      if (acn < 0.0) acn += 2*pi;
      if (acn < ac) acn += 2*pi;
      if (acn - ac > pi) sign = -1;
      else sign = 1;
      nmidx = xc + sign*radius*nrx;
      nmidy = yc + sign*radius*nry;
      if (rooted) {
	xc -= nmidx - midx;
	yc -= nmidy - midy;
      }
      else {
	for (ic=icstart; ; ++ic >= lp->nconnection ? (ic = 0) : 0) {
	  cp = lp->connections[ic];
	  i = cp->start;
	  bases[i].x += nmidx - midx;
	  bases[i].y += nmidy - midy;
	  i = cp->end;
	  bases[i].x += nmidx - midx;
	  bases[i].y += nmidy - midy;
	  if (ic == icend) break;
	}
      }
    }
    icstart = icnext;
    done_all_connections = icstart == icstart1;
  }
  for (ic=0; ic < lp->nconnection; ic++) {
    cp = lp->connections[ic];
    j = ic + 1;
    if (j >= lp->nconnection) j = 0;
    cpnext = lp->connections[j];
    dx = bases[cp->end].x - xc;
    dy = bases[cp->end].y - yc;
    rc = sqrt(dx*dx + dy*dy);
    ac = atan2(dy,dx);
    if (ac < 0.0) ac += 2*pi;
    dx = bases[cpnext->start].x - xc;
    dy = bases[cpnext->start].y - yc;
    rcn = sqrt(dx*dx + dy*dy);
    acn = atan2(dy,dx);
    if (acn < 0.0) acn += 2*pi;
    if (acn < ac) acn += 2*pi;
    dan = acn - ac;
    dcp = cpnext->angle - cp->angle;
    if (dcp <= 0.0) dcp += 2*pi;
    if (fabs(dan-dcp) > pi) {
      if (cp->extruded) {
        vrna_message_warning("...from traverse_loop. Loop %d has crossed regions",
                                    lp->number);
      }
      else if ((cpnext->start - cp->end) != 1) {
	cp->extruded = true;
	goto set_radius;	    /* Forever shamed */
      }
    }
    if (cp->extruded) {
      construct_extruded_segment(cp,cpnext);
    }
    else {
      n = cpnext->start - cp->end;
      if (n < 0) n += nbase + 1;
      angleinc = dan / n;
      for (j = 1;  j < n;  j++) {
	i = cp->end + j;
	if (i > nbase) i -= nbase + 1;
	a = ac + j*angleinc;
	rr = rc + (rcn-rc)*(a-ac)/dan;
	bases[i].x = xc + rr*cos(a);
	bases[i].y = yc + rr*sin(a);
      }
    }
  }
  for (ic=0; ic < lp->nconnection; ic++) {
    if (icroot != ic) {
      cp = lp->connections[ic];
      generate_region(cp);
      traverse_loop(cp->loop,cp);
    }
  }
  n = 0;
  sx = 0.0;
  sy = 0.0;
  for (ic = 0;  ic < lp->nconnection;  ic++) {
    j = ic + 1;
    if (j >= lp->nconnection) j = 0;
    cp = lp->connections[ic];
    cpnext = lp->connections[j];
    n += 2;
    sx += bases[cp->start].x + bases[cp->end].x;
    sy += bases[cp->start].y + bases[cp->end].y;
    if (!cp->extruded) {
      for (j = cp->end + 1;  j != cpnext->start;  j++) {
	if (j > nbase) j -= nbase + 1;
	n++;
	sx += bases[j].x;
	sy += bases[j].y;
      }
    }
  }
  lp->x = sx / n;
  lp->y = sy / n;

  /* free connections (ih) */
  for (ic = 0;  ic < lp->nconnection;  ic++)
    free(lp->connections[ic]);
  free(lp->connections);
}



static void determine_radius(struct loop *lp,double lencut)
/*
*   For the loop pointed to by lp, determine the radius of
*   the loop that will ensure that each base around the loop will
*   have a separation of at least lencut around the circle.
*   If a segment joining two connectors will not support this separation,
*   then the flag, extruded, will be set in the first of these
*   two indicators. The radius is set in lp.
*
*   The radius is selected by a least squares procedure where the sum of the
*   squares of the deviations of length from the ideal value of 1 is used
*   as the error function.
*/

{
  double mindit,ci,dt,sumn,sumd,radius,dit;
  int i,j,end,start,imindit;
  struct connection *cp,*cpnext;
  static double rt2_2 = 0.7071068;

  do {
    mindit = 1.0e10;
    for (sumd=0.0, sumn=0.0, i=0;
	 i < lp->nconnection;
	 i++) {
      cp = lp->connections[i];
      j = i + 1;
      if (j >= lp->nconnection) j = 0;
      cpnext = lp->connections[j];
      end =  cp->end;
      start = cpnext->start;
      if (start < end) start += nbase + 1;
      dt = cpnext->angle - cp->angle;
      if (dt <= 0.0) dt += 2*pi;
      if (!cp->extruded)
	ci = start - end;
      else {
	if (dt <= pi/2) ci = 2.0;
	else ci = 1.5;
      }
      sumn += dt * (1.0/ci + 1.0);
      sumd += dt * dt / ci;
      dit = dt/ci;
      if (dit < mindit && !cp->extruded && ci > 1.0) {
	mindit = dit;
	imindit = i;
      }
    }
    radius = sumn/sumd;
    if (radius < rt2_2) radius = rt2_2;
    if (mindit*radius < lencut) {
      lp->connections[imindit]->extruded = true;
    }
  } while (mindit*radius < lencut);
  if (lp->radius > 0.0)
    radius = lp->radius;
  else lp->radius = radius;
}


static logical    connected_connection(struct connection *cp, struct connection *cpnext)
/*
*   Determines if the connections cp and cpnext are connected
*/

{

  if (cp->extruded) {
    return true;
  }
  else if (cp->end+1 == cpnext->start) {
    return true;
  }
  else {
    return false;
  }
}


static int    find_ic_middle(int icstart, int icend, struct connection *anchor_connection, struct connection *acp, struct loop *lp)
/*
*   Finds the middle of a set of connected connectors. This is normally
*   the middle connection in the sequence except if one of the connections
*   is the anchor, in which case that connection will be used.
*/

{
  int count,ret,ic,i;
  logical done;

  count = 0;
  ret = -1;
  ic = icstart;
  done = false;
  while (!done) {
    if (count++ > lp->nconnection * 2) {
      printf("Infinite loop detected in find_ic_middle\n");
      exit(FATAL_ERROR);
    }
    if (anchor_connection != NULL && lp->connections[ic] == acp) {
      ret = ic;
    }
    done = ic == icend;
    if (++ic >= lp->nconnection) {
      ic = 0;
    }
  }
  if (ret == -1) {
    for (i=1, ic=icstart; i<(count+1)/2; i++) {
      if (++ic >= lp->nconnection) ic = 0;
    }
    ret = ic;
  }
  return ret;
}


static void generate_region(struct connection *cp)
/*
*   Generates the coordinates for the base pairing region of a connection
*   given the position of the starting base pair.
*/

{
  int l,start,end,i,mate;
  struct region *rp;

  rp = cp->region;
  l = 0;
  if (cp->start == rp->start1) {
    start = rp->start1;
    end = rp->end1;
  }
  else {
    start = rp->start2;
    end = rp->end2;
  }
  if (bases[cp->start].x > anum - 100.0 ||
      bases[cp->end].x > anum - 100.0) {
    printf("Bad region passed to generate_region. Coordinates not defined.\n");
    exit(FATAL_ERROR);
  }
  for (i=start+1; i<=end; i++) {
    l++;
    bases[i].x = bases[cp->start].x + l*cp->xrad;
    bases[i].y = bases[cp->start].y + l*cp->yrad;
    mate = bases[i].mate;
    bases[mate].x = bases[cp->end].x + l*cp->xrad;
    bases[mate].y = bases[cp->end].y + l*cp->yrad;
  }
}


static void construct_circle_segment(int start, int end)
/*
*   Draws the segment of residue between the bases numbered start
*   through end, where start and end are presumed to be part of a base
*   pairing region. They are drawn as a circle which has a chord given
*   by the ends of two base pairing regions defined by the connections.
*/

{
  double dx,dy,rr,h,angleinc,midx,midy,xn,yn,nrx,nry,mx,my,a;
  int l,j,i;

  dx = bases[end].x - bases[start].x;
  dy = bases[end].y - bases[start].y;
  rr = sqrt(dx*dx + dy*dy);
  l = end - start;
  if (l < 0) l += nbase + 1;
  if (rr >= l) {
    dx /= rr;
    dy /= rr;
    for (j = 1;  j < l;  j++) {
      i = start + j;
      if (i > nbase) i -= nbase + 1;
      bases[i].x = bases[start].x + dx*(double)j/(double)l;
      bases[i].y = bases[start].y + dy*(double)j/(double)l;
    }
  }
  else {
    find_center_for_arc(l-1,rr,&h,&angleinc);
    dx /= rr;
    dy /= rr;
    midx = bases[start].x + dx*rr/2.0;
    midy = bases[start].y + dy*rr/2.0;
    xn = dy;
    yn = -dx;
    nrx = midx + h*xn;
    nry = midy + h*yn;
    mx = bases[start].x - nrx;
    my = bases[start].y - nry;
    rr = sqrt(mx*mx + my*my);
    a = atan2(my,mx);
    for (j = 1;  j < l;  j++) {
      i = start + j;
      if (i > nbase) i -= nbase + 1;
      bases[i].x = nrx + rr*cos(a+j*angleinc);
      bases[i].y = nry + rr*sin(a+j*angleinc);
    }
  }
}


static void construct_extruded_segment(struct connection *cp, struct connection *cpnext)
/*
*   Constructs the segment between cp and cpnext as a circle if possible.
*   However, if the segment is too large, the lines are drawn between
*   the two connecting regions, and bases are placed there until the
*   connecting circle will fit.
*/

{
  double astart,aend1,aend2,aave,dx,dy,a1,a2,ac,rr,da,dac;
  int start,end,n,nstart,nend;
  logical collision;

  astart = cp->angle;
  aend2 = aend1 = cpnext->angle;
  if (aend2 < astart) aend2 += 2*pi;
  aave = (astart + aend2) / 2.0;
  start = cp->end;
  end = cpnext->start;
  n = end - start;
  if (n < 0) n += nbase + 1;
  da = cpnext->angle - cp->angle;
  if (da < 0.0) {
    da += 2*pi;
  }
  if (n == 2) construct_circle_segment(start,end);
  else {
    dx = bases[end].x - bases[start].x;
    dy = bases[end].y - bases[start].y;
    rr = sqrt(dx*dx + dy*dy);
    dx /= rr;
    dy /= rr;
    if (rr >= 1.5 && da <= pi/2) {
      nstart = start + 1;
      if (nstart > nbase) nstart -= nbase + 1;
      nend = end - 1;
      if (nend < 0) nend += nbase + 1;
      bases[nstart].x = bases[start].x + 0.5*dx;
      bases[nstart].y = bases[start].y + 0.5*dy;
      bases[nend].x = bases[end].x - 0.5*dx;
      bases[nend].y = bases[end].y - 0.5*dy;
      start = nstart;
      end = nend;
    }
    do {
      collision = false;
      construct_circle_segment(start,end);
      nstart = start + 1;
      if (nstart > nbase) nstart -= nbase + 1;
      dx = bases[nstart].x - bases[start].x;
      dy = bases[nstart].y - bases[start].y;
      a1 = atan2(dy,dx);
      if (a1 < 0.0) a1 += 2*pi;
      dac = a1 - astart;
      if (dac < 0.0) dac += 2*pi;
      if (dac > pi) collision = true;
      nend = end - 1;
      if (nend < 0) nend += nbase + 1;
      dx = bases[nend].x - bases[end].x;
      dy = bases[nend].y - bases[end].y;
      a2 = atan2(dy,dx);
      if (a2 < 0.0) a2 += 2*pi;
      dac = aend1 - a2;
      if (dac < 0.0) dac += 2*pi;
      if (dac > pi) collision = true;
      if (collision) {
	ac = minf2(aave,astart + 0.5);
	bases[nstart].x = bases[start].x + cos(ac);
	bases[nstart].y = bases[start].y + sin(ac);
	start = nstart;
	ac = maxf2(aave,aend2 - 0.5);
	bases[nend].x = bases[end].x + cos(ac);
	bases[nend].y = bases[end].y + sin(ac);
	end = nend;
	n -= 2;
      }
    } while    (collision && n > 1);
  }
}


static void find_center_for_arc(int n,double b,double *hp,double *thetap)
/*
*   Given n points to be placed equidistantly and equiangularly on a
*   polygon which has a chord of length, b, find the distance, h, from the
*   midpoint of the chord for the center of polygon. Positive values
*   mean the center is within the polygon and the chord, whereas
*   negative values mean the center is outside the chord. Also, the
*   radial angle for each polygon side is returned in theta.
*
*   The procedure uses a bisection algorithm to find the correct
*   value for the center. Two equations are solved, the angles
*   around the center must add to 2*pi, and the sides of the polygon
*   excluding the chord must have a length of 1.
*/

{
  double h,hhi,hlow,r,disc,theta,e,phi;
  int iter;
#define maxiter 500

  hhi = (n+1) / pi;
  hlow = - hhi - b/(n+1.000001-b);  /* changed to prevent div by zero if (ih) */
  if (b<1) hlow = 0;  /* otherwise we might fail below (ih) */
  iter = 0;
  do {
    h = (hhi + hlow) / 2.0;
    r = sqrt(h*h + b*b/4.0);
    /*  if (r<0.5) {r = 0.5; h = 0.5*sqrt(1-b*b);} */
    disc = 1.0 - 0.5/(r*r);
    if (fabs(disc) > 1.0) {
      vrna_message_error("Unexpected large magnitude discriminant = %g %g", disc,r);
      exit(FATAL_ERROR);
    }
    theta = acos(disc);
    /*    theta = 2*acos(sqrt(1-1/(4*r*r))); */
    phi = acos(h/r);
    e = theta * (n+1) + 2*phi - 2*pi;
    if (e > 0.0) {
      hlow = h;
    }
    else {
      hhi = h;
    }
  } while    (fabs(e) > 0.0001 && ++iter < maxiter);
  if (iter >= maxiter) {
    vrna_message_warning("Iteration failed in find_center_for_arc");
    h = 0.0;
    theta = 0.0;
  }
  *hp = h;
  *thetap = theta;
}

