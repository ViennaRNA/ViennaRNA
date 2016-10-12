/* Expands a blossom.  Fixes up LINK and MATE. */

void
UNPAIR (int oldbase, int oldmate)
{   int e, newbase, u;

#ifdef DEBUG
    printf("Unpair oldbase, oldmate=%d %d\n",oldbase, oldmate);
#endif

    UNLINK (oldbase);
    newbase = BMATE (oldmate);
    if (newbase != oldbase) {
	LINK [oldbase] = -DUMMYEDGE;
	REMATCH (newbase, MATE[oldbase]);
	if (f == LASTEDGE[1])
	    LINK[secondmate] = -LASTEDGE[2];
	else LINK[secondmate] = -LASTEDGE[1];
	}
    e = LINK[oldmate];
    u = BEND (OPPEDGE (e));
    if (u == newbase) {
	POINTER (newbase, oldmate, e);
	return;
	}
    LINK[BMATE (u)] = -e;
    do {
	e = -LINK[u];
	v = BMATE (u);
	POINTER (u, v, -LINK[v]);
	u = BEND (e);
    } while (u != newbase);
    e = OPPEDGE (e);
    POINTER (newbase, oldmate, e);
}


/* changes the matching along an alternating path */
/* firstmate is the first base vertex on the path */
/* edge e is the new matched edge for firstmate   */

void
REMATCH (int firstmate, int e)
{
#ifdef DEBUG
     printf("Rematch firstmate=%d e=%d-%d\n",firstmate, END[OPPEDGE(e)], END[e]);
#endif

    MATE[firstmate] = e;
    nexte = -LINK[firstmate];
    while (nexte != DUMMYEDGE) {
	e = nexte;
	f = OPPEDGE (e);
	firstmate = BEND (e);
	secondmate = BEND (f);
	nexte = -LINK[firstmate];
	LINK[firstmate] = -MATE[secondmate];
	LINK[secondmate] = -MATE[firstmate];
	MATE[firstmate] = f;
	MATE[secondmate] = e;
	}
}


/* unlinks subblossoms in a blossom.  oldbase is the base of the blossom to */
/* be unlinked. */

void
UNLINK (int oldbase)
{   int k, j=1;

#ifdef DEBUG
    printf("Unlink oldbase=%d\n",oldbase);
#endif

    i = NEXTVTX[oldbase];
    newbase = NEXTVTX[oldbase];
    nextbase = NEXTVTX[LASTVTX[newbase]];
    e = LINK[nextbase];
UL2:
    do {
	nextedge = OPPEDGE (LINK[newbase]);
	for (k=1; k <= 2; ++k) {
	    LINK[newbase] = -LINK[newbase];
	    BASE[i] = newbase;
	    i = NEXTVTX[i];
	    while (i != nextbase) {
		BASE[i] = newbase;
		i = NEXTVTX[i];
		} 
	    newbase = nextbase;
	    nextbase = NEXTVTX[LASTVTX[newbase]];
	    }
	} while (LINK[nextbase] == nextedge);
    if (j==1) {
	LASTEDGE[1] = nextedge;
	j++;
	nextedge = OPPEDGE (e);
	if (LINK[nextbase] == nextedge) {
	    goto UL2;
	    }
	}
    LASTEDGE[2] = nextedge;

    if (BASE[LASTVTX[oldbase]] == oldbase) 
	NEXTVTX[oldbase] = newbase;
    else {
	NEXTVTX[oldbase] = DUMMYVERTEX;
	LASTVTX[oldbase] = oldbase;
	}
}


