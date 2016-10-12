/* set up data structures for weighted match */

/* to add a new type, add new case in SetUp() and a Set_X() routine */

SetUp (gptr,type)
int gptr,type;

{   int i,allocsize;
    Graph g;
    EuclidGraph xy;
    MatrixGraph matg;

    if (type==1) {
	g = (Graph) gptr;
	U = Degree(g,0);
	V = NumEdges(g);
	}
    else if (type==2) {
	xy = (EuclidGraph) gptr;
	U = xy[0][0];
	V = U*(U-1)/2;
	}
    else if (type==3) {
	matg = (MatrixGraph) gptr;
	U = matg[0];
	V = U*(U-1)/2;
	}

    allocsize = (U+2*V+2)*sizeof(int);
    A      = (int *) malloc(allocsize);
    END    = (int *) malloc(allocsize);
    WEIGHT = (int *) malloc(allocsize);
    for (i=0; i<U+2*V+2; i++)
	A[i]=END[i]=WEIGHT[i]=0;

    if (type == 1) SetStandard(g);
    else if (type == 2) SetEuclid(xy);
    else if (type == 3) SetMatrix(matg);
}


/* set up from Type 1 graph. */

SetStandard(graph)
Graph graph;
{   int elabel, adj_node, i, j;
    int u, v, currentedge;
    Edge edge;

    currentedge = U+2;
    for (i=1; i<=U; ++i) {
	edge = FirstEdge(graph,i);
	for (j = 1; j <= Degree(graph,i); ++j) {
	    adj_node = EndPoint(edge);
	    if (i < adj_node) {
		elabel = ELabel(edge)*2;
		WEIGHT[currentedge-1] = WEIGHT[currentedge] = 2*elabel;
		END[currentedge-1] = i;
		END[currentedge] = adj_node;
		if (A[i] == 0)
		    A[i] = currentedge;
		else {
		    u = i;
		    v = A[i];
		    while (v != 0) {
			if (END[v] > adj_node)
			    break;
			u = v;
			v = A[v];
			}
		    A[u] = currentedge;
		    A[currentedge] = v;
		    }
		u = adj_node;
		v = A[u];
		while (v != 0) {
		    u = v;
		    v = A[v];
		    }
		A[u] = currentedge - 1;
		currentedge += 2;
		}
	    edge = NextEdge(edge);
	    }
	}
}

/* set up from Euclidean graph */

SetEuclid(graph)
EuclidGraph graph;
{   int i,j,currentedge;

    currentedge = U+2;

    for (i=U; i>=1; --i) 
	for (j = i-1; j >= 1; --j) {
	    WEIGHT[currentedge-1] = WEIGHT[currentedge]
		    = 2*eucdist2(graph,i,j);
	    END[currentedge-1] = i;
	    END[currentedge] = j;
	    A[currentedge] = A[i];
	    A[i] = currentedge;
	    A[currentedge-1] = A[j];
	    A[j] = currentedge-1;
	    currentedge += 2;
	    }
}

SetMatrix(graph)
MatrixGraph graph;
{   int i,j,currentedge;

    currentedge = U+2;

    for (i=U; i>=1; --i) 
	for (j = i-1; j >= 1; --j) {
	    WEIGHT[currentedge-1] = WEIGHT[currentedge]
		    = 2*graph[j*U+i];
	    END[currentedge-1] = i;
	    END[currentedge] = j;
	    A[currentedge] = A[i];
	    A[i] = currentedge;
	    A[currentedge-1] = A[j];
	    A[j] = currentedge-1;
	    currentedge += 2;
	    }
}

