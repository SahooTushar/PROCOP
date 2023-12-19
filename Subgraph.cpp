#include "Subgraph.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>

//Graph* Subgraph::subg;
//vector<int>* Subgraph::embed;
extern std::vector<Graph*> graphv;

Subgraph::Subgraph(int n) {
  numVertex = n;
  numEdge = 0;
  frequency=0;
  embed=NULL;
  g = NULL;
  canstr=NULL;
}

Subgraph::~Subgraph() {
  if (canstr != NULL) delete[] canstr;
  delete g;
}
void Subgraph::subgraphCensus(Graph *g) {
        int sz=numVertex;
        //cout << canstr << endl;
        //cout << g->numNodes() << endl;
	//cout << numVertex << endl;

        Graph *h = new Graph();
	//if (dir)  GraphUtils::strToGraph(g, canstr, size, DIRECTED);
  	//else      GraphUtils::strToGraph(g, canstr, size, UNDIRECTED);
        Graph::strToGraph(h, canstr, sz);
	numEdge = h->numberEdges();
	h->sortNeighbours();
    h->makeNeighboursArray();
	 // sort and create array of neighbours
  	//h->sortNeighbours();
   	//cout << "start" << endl;
	int bestorder[sz];
	int bestneighbours[sz][10];
	int bestconditions[sz][10];
        //cout << "canstr=" << canstr << endl;
	setBest(h, bestorder, bestneighbours, bestconditions);
	/*for (int i=0; i< sz; i++) {
		cout << bestorder[i] << " : ";
		for (int j=1; j<= bestneighbours[i][0]; j++) {
			cout << bestneighbours[i][j] << " ";
		}
		cout << " : ";
		for (int j=1; j<= bestconditions[i][0]; j++) {
			cout << bestconditions[i][j] << " ";
		}
		cout << endl;
	}*/
	//cout << "end" << endl;
	int map1[sz];
	int marked[g->numberNodes()];
	for (int i=0; i < g->numberNodes() ; i++) {
		marked[i]=0;
	}
        embed = new std::list< std::vector<int> >;
	for (int i=0; i < g->numberNodes() ; i++) {
		if (g->numberNeighbours(i) >= h->numberNeighbours(bestorder[0])) {
			//cout << "i=" << i << endl;
			map1[bestorder[0]] = i;
			marked[map1[bestorder[0]]]=1;
			this->expand(g,h,0,bestorder, bestneighbours, bestconditions,map1,marked);
		}
	}
	std::cout << "frequency=" << frequency << std::endl;
}
void Subgraph::expand(Graph *g,Graph *h,int l,int bestorder[],int bestneighbours[][10],int bestconditions[][10],int map1[], int marked[]) {
	l++;
	int i, j;
	if (l == h->numberNodes()) {
                // non induced subgraphs
		//frequency++;
                std::vector<int> emb;
                for (i = 0; i<l; i++) {
                   emb.push_back(map1[i]);
                }
                embed->push_back(emb);
                /*for (i=0; i< l; i++) {
			cout << map[i] << " ";
		}
		cout << endl;*/

		// Induced subgraphs
		for (i=0; i< l; i++) {
			for (j=i+1; j< l; j++) {
				if(g->hasEdge(map1[i], map1[j])==true && h->hasEdge(i, j) == false) break;
			}
			if (j!=l) break;
		}
		if (i==l) {
            frequency++;
            //for (i=0; i< l; i++) {
            //    std::cout << map1[i] << " ";
            //}
            //std::cout << std::endl;
            Graph *grph = new Graph();
            Graph::strToGraph(grph, canstr, numVertex);
            grph->setMap(map1);
            graphv.push_back(grph);
		}
		marked[map1[bestorder[l-1]]]=0;
		return;
	}
	/*for (int i=0; i< l; i++) {
		cout << map[i] << " ";
	}
	cout << endl;*/
	int n=map1[bestneighbours[l][1]];
	int *temp = g->getNeighboursArray(n);
	//cout << "g->numberNeighbours(n)=" << n << " " << g->numberNeighbours(n) << endl;
	for (i=0; i< g->numberNeighbours(n); i++) {
		if (marked[temp[i]])
			continue;
		if (g->numberNeighbours(temp[i]) < h->numberNeighbours(bestorder[l])) {
			continue;
		}
		for (j=2; j<=bestneighbours[l][0]; j++) {
			if (g->hasEdge(temp[i], map1[bestneighbours[l][j]]) == false) {
				break;
			}
		}
		if (j <= bestneighbours[l][0]) {
			continue;
		}
		for (j=1; j<=bestconditions[l][0]; j++) {
			if (temp[i] < map1[bestconditions[l][j]]) {
				break;
			}
		}
		if (j <= bestconditions[l][0]) {
			continue;
		}

		map1[bestorder[l]] = temp[i];
		marked[map1[bestorder[l]]]=1;
		expand(g,h,l,bestorder, bestneighbours, bestconditions,map1,marked);
	}
	marked[map1[bestorder[l-1]]]=0;
	return;
}


void Subgraph::setBest(Graph *h, int bestorder[], int bestneighbours[][10], int bestconditions[][10]) {
	int sz=h->numberNodes();
	std::list<std::pair<int,int>> *cond = new std::list<std::pair<int,int>>;
  	GraphIsomor::symmetryConditions(h, cond);
	//cout << "cond->size()= " << cond->size() << endl;
	int i,j;
	int prid[sz];
	int mark[sz];
	for (i=0; i< sz; i++)
		mark[i]=0;
	selectionSort(h,prid);
	int best=0;
	for (i=1; i< sz; i++) {
		//cout << prid[i] << "  ";
		if (prid[i] < prid[best])
			best = i;
	}
	//cout << endl;
	bestorder[0] = best;
	bestneighbours[0][0]=0;
	bestconditions[0][0]=0;
	mark[best] = 1;
	int l=1;
	for (i=1; i< sz; i++) {
		int nc[sz];
		for (j=0; j< sz; j++)
			nc[j]=0;
		for (j=0; j< sz; j++) {
			if (mark[j] == 1)
				continue;
			for (int k=0; k< l; k++) {
				if (h->hasEdge(j,bestorder[k])) {
					nc[j]++;
				}
			}
		}
		best = -1;
		for (j=0; j< sz; j++) {
			if (mark[j] == 1)
				continue;
			else if (best == -1)
				best=j;
			else if (nc[j] > nc[best])
				best=j;
			else if (nc[j] == nc[best] && prid[j] < prid[best])
				best=j;
		}
		bestorder[l] = best;
		int k1=0;
		for (int k=0; k< l; k++) {
			if (h->hasEdge(best,bestorder[k])) {
				k1++;
				bestneighbours[l][k1]=bestorder[k];
			}
		}
		bestneighbours[l][0]=k1;
		std::list<std::pair<int,int>>::iterator jj;
		//cout << "cond->size()= " << cond->size() << endl;
		k1=0;
		for (jj=cond->begin(); jj!=cond->end(); jj++) {
            		//printf("%d < %d  |  ", jj->first, jj->second);
			if (jj->second == best) {
				for (int k=0; k< l; k++) {
					if (jj->first == bestorder[k]) {
						k1++;
						bestconditions[l][k1]=bestorder[k];
					}
				}
			}
       		}
		bestconditions[l][0]=k1;
		mark[best] = 1;
		l++;
	}
}

void Subgraph::selectionSort(Graph *h,int prid[]) {
	int sz= h->numberNodes();
	int i,j;
	int d[sz];
	int nd[sz];
	int mark[sz];
	for (i=0; i< sz; i++) {
		d[i]=h->numberNeighbours(i);
		int *temp = h->getNeighboursArray(i);
		nd[i] = 0;
		for (j=0; j< h->numberNeighbours(i); j++) {
			nd[i] = nd[i] + h->numberNeighbours(temp[j]);
		}
	}
	for (i=0; i< sz; i++)
		mark[i]=0;
	for (i=0; i< sz; i++) {
		int opt=-1;
		for (j=0; j< sz; j++) {
			if (mark[j] == 1)
				continue;
			else if (opt == -1) {
				opt =j;
				continue;
			}
			else if (d[j] > d[opt]) {
				opt =j;
			}
			else if(d[j] == d[opt] && nd[j]>nd[opt]) {
				opt =j;
			}
		}
		prid[opt]=i;
		mark[opt] =1;
	}
}



