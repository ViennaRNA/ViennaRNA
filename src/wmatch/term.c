/* updates numerical bounds for linking paths. */
/* called with LAST_D set to the bound on DELTA for the next search */

void
SET_BOUNDS ()
 
{   int del;

    for (v=1; v <= U; ++v) {
	if (LINK[v] < 0 || BASE[v] != v) {
	    NEXT_D[v] = LAST_D;
	    continue;
	    }
	LINK[v] = -LINK[v];
	i = v;
	while (i != DUMMYVERTEX) {
	    Y[i] -= DELTA;
	    i = NEXTVTX[i];
	    }
	f = MATE[v];
	if (f != DUMMYEDGE) {
	    i = BEND(f);
	    del = SLACK(f);
	    while (i != DUMMYVERTEX) {
		Y[i] -= del;
		i = NEXTVTX[i];
		}
	    }
	NEXT_D[v] = LAST_D;
	}
}


/* undoes all blossoms to get the final matching */

void
UNPAIR_ALL ()

{   int u;

    for (v=1; v <= U; ++v) {
	if (BASE[v] != v || LASTVTX[v] == v)
	    continue;
	nextu = v;
	NEXTVTX[LASTVTX[nextu]] = DUMMYVERTEX;
	while (1) {
	    u = nextu;
	    nextu = NEXTVTX[nextu];
	    UNLINK (u);
	    if (LASTVTX[u] != u) {
		f = (LASTEDGE[2] == OPPEDGE(e)) ? LASTEDGE[1] : LASTEDGE[2];
		NEXTVTX[LASTVTX[BEND(f)]] = u;
		}
	    newbase = BMATE (BMATE(u));
	    if (newbase != DUMMYVERTEX && newbase != u) {
		LINK[u] = -DUMMYEDGE;
		REMATCH (newbase, MATE[u]);
		}
	    while (LASTVTX[nextu] == nextu && nextu != DUMMYVERTEX)
		    nextu = NEXTVTX[nextu];
	    if (LASTVTX[nextu] == nextu && nextu == DUMMYVERTEX)
		break;
	    }
	}
}


