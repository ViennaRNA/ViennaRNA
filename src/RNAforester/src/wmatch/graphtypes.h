#ifndef _GRAPHTYPES_H_
#define _GRAPHTYPES_H_

#define INF	100000000

#ifndef NULL
#define NULL	0
#endif

struct node_entry {
    int degree;
    int label;
    int x;
    int y;
    struct edge_ent *adj_vector;
};
typedef struct node_entry *Graph;

struct edge_ent {
    int endpoint;
    int label;
    int label2;
    struct edge_ent *nextedge;
    struct edge_ent *prevedge;
    struct edge_ent *otheredge;
};
typedef struct edge_ent *Edge;

#define Degree(graph,n)    (graph[n].degree)
#define NLabel(graph,n)    (graph[n].label)
#define Xcoord(graph,n)    (graph[n].x)
#define Ycoord(graph,n)    (graph[n].y)
#define FirstEdge(graph,n) (graph[n].adj_vector)

#define EndPoint(e) (e->endpoint)
#define ELabel(e)   (e->label)
#define ELabel2(e)  (e->label2)
#define Other(e)    (e->otheredge)
#define NextEdge(e) (e->nextedge)


extern Graph Prim();
//extern int *EulerTraverse(),*Match(),*Weighted_Match(),*Dijkstra(),*Concomp();

/* Euclidean graph type */
typedef int (*EuclidGraph)[2];

extern Graph EuclidPrim();
extern EuclidGraph ReadEuclid(),NewEuclid();
extern int eucdist(),eucdistsq();

extern int *CvxHull();

/* Distance matrix graph type */
typedef int *MatrixGraph;

extern int *MatrixDijkstra();
extern Graph MatrixPrim();
extern Graph MatrixMaxFlow();
extern MatrixGraph ReadMatrix(), NewMatrix();

#ifdef __cplusplus
extern "C" {
#endif

    void AddEdge (Graph g,int n,int m,int label);
    Edge FindEdge(Graph graph,int i,int j);
    int RemoveEdge(Graph graph,Edge edge);
    int NumEdges(Graph g);
    Graph NewGraph(int size);
    EuclidGraph NewEuclid(int size);
    MatrixGraph NewMatrix(int size);
    Graph CopyGraph(Graph g);
    Graph ReadGraph (int *size,char *file);
    void WriteGraph (Graph graph,char *file);
    EuclidGraph ReadEuclid(int *size,char *file);
    void WriteEuclid(EuclidGraph graph,char *file);
    MatrixGraph ReadMatrix(int *size,char *file);
    void WriteMatrix(MatrixGraph graph,char *file);
    int eucdist (EuclidGraph graph,int i,int j);
    int eucdist2 (EuclidGraph graph,int i,int j);
    int eucdistsq(EuclidGraph graph,int i,int j);

#ifdef __cplusplus
}
#endif

#endif
