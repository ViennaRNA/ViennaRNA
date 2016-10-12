/* Process an edge linking two linked vertices */
/* Note: global variable v set to the base of one end of the linking edge */

void
PAIR(int *outcome)
{   int u, w, temp;

#ifdef DEBUG
    printf("Pair v=%d\n",v);
#endif

    e = NEXTEDGE[v];
    while (SLACK(e) != 2*DELTA)
	e = NEXTPAIR[e];
    w = BEND (e);
    LINK[BMATE (w)] = -e;
    u = BMATE (v);
    while (LINK[u] != -e) {
	LINK[u] = -e;
	if (MATE[w] != DUMMYEDGE) {
	    temp = v;
	    v = w;
	    w = temp; }
	v = BLINK (v);
	u = BMATE (v);
	}
    if (u == DUMMYVERTEX && v != w) {
	*outcome = 1;
	return;
	}
    newlast = v;
    newbase = v;
    oldfirst = NEXTVTX[v];
    LINK_PATH (e);
    LINK_PATH (OPPEDGE (e));
    NEXTVTX[newlast] = oldfirst;
    if (LASTVTX[newbase] == newbase)
	LASTVTX[newbase] = newlast;
    NEXTPAIR[DUMMYEDGE] = DUMMYEDGE;
    MERGE_PAIRS (newbase);
    i = NEXTVTX[newbase];
    do {
	MERGE_PAIRS (i);
	i = NEXTVTX[LASTVTX[i]];
	SCAN (i, 2*DELTA - SLACK(MATE[i]));
	i = NEXTVTX[LASTVTX[i]];
    } while (i != oldfirst);
    *outcome = 0;
    return;
}


/* merges a subblossom's pair vector into a new blossom's pair vector */
/* v is the base of the previously unlinked subblossom */
/* Note: global variable newbase set to the base of the new blossom */
/* 	called with NEXTPAIR[DUMMYEDGE] pointing to the first edge */
/*		on newbase's pair vector */

void
MERGE_PAIRS(int v)
{
#ifdef DEBUG
    printf("Merge Pairs v=%d\n",v);
#endif

    NEXT_D[v] = LAST_D;
    pairpoint = DUMMYEDGE;
    f = NEXTEDGE[v];
    while (f != DUMMYEDGE) {
	e = f;
	neighbor = END[e];
	f = NEXTPAIR[f];
	if (BASE[neighbor] != newbase)
	    INSERT_PAIR();
	}
}


/* links the unlinked vertices in the path P(END[e],newbase) */
/* Note: global variable newbase is set to the base vertex of the new blossom */
/*		newlast is set to the last vertex in newbase's current blossom*/

void
LINK_PATH (int e)
{   int u;

#ifdef DEBUG
    printf("Link Path e=%d-%d\n", END[OPPEDGE(e)], END[e]);
#endif

    v = BEND (e);
    while (v != newbase) {
	u = BMATE (v);
	LINK[u] = OPPEDGE (e);
	NEXTVTX[newlast] = v;
	NEXTVTX[LASTVTX[v]] = u;
	newlast = LASTVTX[u];
	i = v;
	BASE[i] = newbase;
	i = NEXTVTX[i];
	while (i != DUMMYVERTEX) {
	    BASE[i] = newbase;
	    i = NEXTVTX[i];
	    }
	e = LINK[v];
	v = BEND (e);
	}
}


/* Update a blossom's pair vector. */
/* Note: called with global variable e set to the edge to be inserted. */
/*			neighbor set to the vertex at the end of e */
/*			pairpoint set to the next pair on the pair vector */

void
INSERT_PAIR (void)
{   int del_e;

#ifdef DEBUG
    printf("Insert Pair e=%d-%d\n",END[OPPEDGE(e)],END[e]);
#endif

    del_e = SLACK(e)/2;
    nextpoint = NEXTPAIR[pairpoint];

    while (END[nextpoint] < neighbor) {
	pairpoint = nextpoint;
	nextpoint = NEXTPAIR[nextpoint];
	}
    if (END[nextpoint] == neighbor) {
	if (del_e >= SLACK (nextpoint)/2)
	    return;
	nextpoint = NEXTPAIR[nextpoint];
	}
    NEXTPAIR[pairpoint] = e;
    pairpoint = e;
    NEXTPAIR[e] = nextpoint;
    if (NEXT_D[newbase] > del_e)
	NEXT_D[newbase] = del_e;
}

