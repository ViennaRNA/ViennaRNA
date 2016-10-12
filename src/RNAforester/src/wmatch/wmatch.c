/* N-cubed weighted matching */
/* Implementation of H. Gabow's Ph.D. thesis, Stanford Univ. 1973 */
/* Written by Edward Rothberg  7/85 */
/* For complete details, please refer to the original paper */

/* Send either a Euclidean graph or an adjacency vector graph. */
/* Returns an array, the ith entry being the mate of vertex 'i'. */
/* A zero entry indicates an unmatched vertex. */

/* To add new type, see readgraph.c */

#include <stdlib.h>
#include <stdio.h>
#include "match.defs"
#include "graphtypes.h"
#include "wmatch.h"

void MERGE_PAIRS(int v);
void LINK_PATH(int e);
void PAIR(int *outcome);
void INSERT_PAIR (void);

void UNPAIR(int oldbase, int oldmate);
void REMATCH (int firstmate, int e);
void UNLINK (int oldbase);

void POINTER (int u, int v, int e);
void SCAN (int x, int del);

#include "pairs.c"
#include "pointer.c"
#include "readgraph.c"
#include "term.c"
#include "unpairs.c"

void Initialize(int maximize);
void FreeUp();

int *Weighted_Match (void *gptr,int type,int maximize)
     //int gptr,type,maximize;
{   int g, j, w, outcome, loop=1;

    /* set up internal data structure */
    SetUp(gptr,type);
    Initialize(maximize);

    for(;;) {
	/* printf("Augment #%d\n",loop++); */
	DELTA = 0;
	for (v=1; v<=U; ++v)
	    if (MATE[v] == DUMMYEDGE)
		POINTER (DUMMYVERTEX, v, DUMMYEDGE);
	for(;;) {
	    i = 1;
	    for (j=2; j<=U; ++j)
		if (NEXT_D[i] > NEXT_D[j])
		    i = j;
	    DELTA = NEXT_D[i];
	    if (DELTA == LAST_D)
		goto done;
	    v = BASE[i];
	    if (LINK[v] >= 0) {
		PAIR (&outcome);
		if (outcome == 1)
		    break;
		}
	    else {
		w = BMATE (v);
		if (LINK[w] < 0) {
		    POINTER (v, w, OPPEDGE(NEXTEDGE[i]));
		    }
		else UNPAIR (v, w);
	        }
	    }

	LAST_D -=DELTA;
	SET_BOUNDS();
	g = OPPEDGE(e);
	REMATCH (BEND(e), g);
	REMATCH (BEND(g), e);
	}
    
done:
    SET_BOUNDS();
    UNPAIR_ALL();
    for (i=1; i<=U;++i) {
	MATE[i] = END[MATE[i]];
	if (MATE[i]==DUMMYVERTEX)
	    MATE[i]=0;
	}

    FreeUp();
    return(MATE);
}


void Initialize(int maximize)
{   int i, allocsize, max_wt= -MAXWT, min_wt=MAXWT;

    DUMMYVERTEX = U+1;
    DUMMYEDGE = U+2*V+1;
    END[DUMMYEDGE] = DUMMYVERTEX;

    for (i=U+2; i<=U+2*V; i+=2) {
	if (WEIGHT[i] > max_wt)
	    max_wt = WEIGHT[i];
	if (WEIGHT[i] < min_wt)
	    min_wt = WEIGHT[i];
	}
    if (!maximize) {
	if (U%2!=0) {
	    printf("Must have an even number of vertices to do a\n");
	    printf("minimum complete matching.\n");
	    exit(0);
	    }
	max_wt += 2; /* Don't want all zero weight */
	for (i=U+1; i<=U+2*V; i++)
	    WEIGHT[i] = max_wt-WEIGHT[i];
	max_wt = max_wt-min_wt;
	}
    LAST_D = max_wt/2;

    allocsize = (U+2)*sizeof(int);
    MATE     = (int *) malloc(allocsize);
    LINK     = (int *) malloc(allocsize);
    BASE     = (int *) malloc(allocsize);
    NEXTVTX  = (int *) malloc(allocsize);
    LASTVTX  = (int *) malloc(allocsize);
    Y        = (int *) malloc(allocsize);
    NEXT_D   = (int *) malloc(allocsize);
    NEXTEDGE = (int *) malloc(allocsize);
    allocsize = (U+2*V+2)*sizeof(int);
    NEXTPAIR = (int *) malloc(allocsize);

    for (i = 1; i <= U+1; ++i) {
	MATE[i] = DUMMYEDGE;
	NEXTEDGE[i] = DUMMYEDGE;
	NEXTVTX[i] = 0;
	LINK[i] = -DUMMYEDGE;
	BASE[i] = i;
	LASTVTX[i] = i;
	Y[i] = LAST_D;
	NEXT_D[i] = LAST_D;
	}
}

void FreeUp()
{
    free(LINK);
    free(BASE);
    free(NEXTVTX);
    free(LASTVTX);
    free(Y);
    free(NEXT_D);
    free(NEXTEDGE);
    free(NEXTPAIR);
    free(A);
    free(END);
    free(WEIGHT);
}

