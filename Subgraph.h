#ifndef _SUBGRAPH_
#define _SUBGRAPH_
#include "Graph.h"
#include "GraphIsomor.h"
#include <vector>
#include <stdio.h>
#include <string>
#include <map>
#include <iostream>

class Subgraph {

public:
    int numVertex; // Number of verices in the subgraph
	int numEdge;      // Number of edges in the subgraph
	int frequency;   // Frequency of embeddings
	char *canstr;
	Graph *g;
	//list<VInt> *embed;
	std::list< std::vector<int> > *embed;
	void subgraphCensus(Graph *g);
	Subgraph(int n);   // Create a subgraph with n vertices
	static void selectionSort(Graph *h,int prid[]);
	static void setBest(Graph *h, int bestorder[], int bestneighbours[][10], int bestconditions[][10]);
	void expand(Graph *g,Graph *h,int l,int bestorder[],int bestneighbours[][10],int bestconditions[][10],int map[],int marked[]);
	~Subgraph();
};
#endif
