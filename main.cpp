// Copyright © 2019 Tushar R. sahoo/ Sabyasachi Patra, IIIT Bhubaneswar, All rights reserved.
// Main Function
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cmath>
#include <iostream>
// using namespace std;
#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "Graph.h"
#include "Subgraph.h"
#include "GraphIsomor.h"
#include "Timer.h"
#include "SubgraphTree.h"
#include "ExpansionTree.h"
#include "DET.h"

#define MAX_BUF 256  // Maximum string buffer size

void parse_cmdline(int argc, char **argv);
void check_inputs();
void initialize();
void prepare_graph();
void compute_real();
void compute_results();
void run_modet();
void run_mdet();
void run_subgraphs();
void read_subgraphs();
void printDenseSubgraphs();
void printComplexes();
void computeScore();
void buildHeap();
bool greaterG(Graph* g1, Graph *g2);
void prepareEdgeList();
void mergeGraph();
bool checkOverlap(Graph* g1, Graph *g2);
double intadjustcd(Graph* g1, int a, Graph *g2, int b);
void club(Graph* g1, Graph *g2, Graph *g3);
void insertHeap(Graph* maxHeapv[], int heapSize, Graph *gr);
int compare_results(const void *a, const void *b);

// Variable declarations
char graph_file[MAX_BUF];
char method[MAX_BUF];
char subgraphs_file[MAX_BUF];
int motif_size;
int rand_number;
int num_exchanges;
int num_tries;
int threshold;
double overlap_threshold;
Graph *g;
SubgraphTree realSG;
SubgraphTree realDSG;
SubgraphTree sg_original;
ExpansionTree *realET;
ExpansionTree *randET;
DET *realDET;
DET *randDET;
std::vector<Subgraph*> sgv;
// graph vectors
std::vector<Graph*> graphv;
std::vector<Graph*> complexes;


int main(int argc, char **argv) {
    std::cout << "This program is developed by Sabyasachi Patra, ";
    std::cout << "IIIT Bhubaneswar, India" << std::endl;
    initialize();
    parse_cmdline(argc, argv);
    check_inputs();
    prepare_graph();
    compute_real();
    compute_results();
    GraphIsomor::finishNauty();
    return 0;
}

// Initialize everything
void initialize() {
    motif_size = -1;
    rand_number = -1;
    num_exchanges = 2;
    num_tries = 7;
    threshold = 0;
    overlap_threshold=1.5;
    snprintf(method, strlen("none")+1, "%s\n", "none");
    snprintf(graph_file, strlen("none")+1, "%s\n", "none");
    FILE *fp;
    fp = fopen("Output.txt", "w");
    fprintf(fp, "OUTPUT");
    fprintf(fp, "\n");
    fclose(fp);
    fp = fopen("Error.txt", "w");
    fprintf(fp, "ERROR");
    fprintf(fp, "\n");
    fclose(fp);
    GraphIsomor::initNauty(motif_size);
}

// Parse all command line arguments
void parse_cmdline(int argc, char **argv) {
    for (int i=1; i < argc; i++) {
        // Graph file
        if (!strcmp("-g", argv[i]) || !strcmp("--graph", argv[i])) {
            snprintf(graph_file, strlen(argv[i+1])+1, "%s\n", argv[i+1]);
            // cout << "graph_file = " << graph_file << endl;
            i++;
        } else if (!strcmp("-s", argv[i]) || !strcmp("--size", argv[i])) {
            motif_size = atoi(argv[++i]);  // Size of motifs to consider
        } else if (!strcmp("-m", argv[i]) || !strcmp("--method", argv[i])) {
            snprintf(method, strlen(argv[i+1])+1, "%s\n", argv[i+1]);
            // cout << "method = " << method << endl;
            if (strcmp(method, "subgraph") == 0) {
                    strcpy(subgraphs_file, argv[i+2]);
                    i++;
            }
            i++;
        } else if (!strcmp("-r", argv[i]) || !strcmp("--random", argv[i])) {
            rand_number = atoi(argv[++i]);  // Method for set random graphs
        } else if (!strcmp("-e", argv[i]) || !strcmp("--exchanges", argv[i])) {
            num_exchanges = atoi(argv[++i]);  // Number of exchanges per edge
        } else if (!strcmp("-t", argv[i]) || !strcmp("--tries", argv[i])) {
            num_tries = atoi(argv[++i]);  // Number of tries per node
        } else if (!strcmp("-th", argv[i]) || !strcmp("--threshold", argv[i])) {
            threshold = atoi(argv[++i]);  // Threshold value
        }
        else if (!strcmp("-o", argv[i]) || !strcmp("--overlapthreshold", argv[i])) {
            overlap_threshold = std::stof(argv[++i]);  // Threshold score
            //std::cout << overlap_threshold << "\n";
        }
    }
}

void check_inputs() {
    if (strcmp(method, "modet") != 0 && strcmp(method, "mdet") != 0 && strcmp(method, "subgraph") != 0) {
        std::cout << "invalid method" << std::endl;
        exit(1);
    }
    if (strcmp(method, "modet") == 0) {
        if (motif_size < 3 || motif_size > 10) {
            std::cout << "invalid motf size" << std::endl;
            exit(1);
        }
    }
    if (strcmp(method, "mdet") == 0) {
        if (motif_size < 3 || motif_size > 15) {
            std::cout << "invalid motf size" << std::endl;
            exit(1);
        }
    }if (strcmp(method, "subgraph") == 0) {
        if (motif_size < 3 || motif_size > 10) {
            std::cout << "invalid motf size" << std::endl;
            exit(1);
        }
    }

    if (strcmp(graph_file, "none") == 0) {
        std::cout << "no input graph file" << std::endl;
        exit(1);
    }
}

// Prepare the real graph for computation
void prepare_graph() {
    g = new Graph();
    // Read the graph file
    Graph::readGraphFile(g, graph_file);
    // sort and create neighbours array
    g->sortNeighbours();
    g->makeNeighboursArray();
    // Print chosen parameters
    printf("motif size: %d\n", motif_size);
    printf("graph file: %s\n", graph_file);
    printf("%d nodes, %d edges\n", g->numberNodes(), g->numberEdges());
    threshold = static_cast<int>(g->numberNodes()*threshold/100.0);
    std::cout << "threshold = " << threshold << std::endl;
}

// Count subgraphs on real network
void compute_real() {
    // Print method name
    if (strcmp(method, "modet") == 0)
        printf("Method: MODET on real network\n");
    else if (strcmp(method, "mdet") == 0)
        printf("Method: MDET on real network\n");
    else if (strcmp(method, "subgraph") == 0)
        printf("Method: SUBGRAPH with subgraphs read from file\n");
    else
        printf("Invalid method\n");

    // Compute frequency
    printf("\nCounting subgraph frequency on 'REAL NETWORK'\n");
    Timer::start(0);
    if (strcmp(method, "modet") == 0)
        run_modet();
    else if (strcmp(method, "mdet") == 0)
        run_mdet();
    else if (strcmp(method, "subgraph") == 0)
        run_subgraphs();
    Timer::stop(0);
    if (strcmp(method, "modet") == 0) {
        printf("%d subgraphs, ", realSG.countSubgraphs());
        printf("%f embeddings\n", realSG.countEmbeddings());
    } else if (strcmp(method, "mdet") == 0) {
        printf("%d subgraphs, ", realDSG.countSubgraphs());
        printf("%f embeddings\n", realDSG.countEmbeddings());
    }
    printf("Time elapsed: %.6fs\n", Timer::elapsed(0));
}

// Run MODET method on graph 'g' and store results on SubgraphTree 'realSG'
void run_modet() {
    Timer::start(0);
    SubgraphTree isomorSG;
    realET = new ExpansionTree();
    realET->create(motif_size, &isomorSG);
    Timer::stop(0);
    printf("Creation time: %.2f\n", Timer::elapsed(0));
    Timer::start(0);
    realET->census(motif_size, g, &realSG);
    // realET->printEmbeddings(motif_size);
    Timer::stop(0);
    printf("census time: %.2f\n", Timer::elapsed(0));
}

// Run MDET method on graph 'g' and store results on SubgraphTree 'realDSG'
void run_mdet() {
    SubgraphTree isomorSG;
    realDET = new DET();
    Timer::start(0);
    realDET->census(motif_size, g, &realDSG, &isomorSG, threshold);
    // realDET->printEmbeddings(motif_size);
    Timer::stop(0);
    printf("census time: %.2f\n", Timer::elapsed(0));
}

// Run SUBGRAPHS method on graph 'g' and store results on SubgraphTree 'realG'
void run_subgraphs() {
  Timer::start(0);
  read_subgraphs();
  Timer::stop(0);
  printf("Creation time: %.2f\n", Timer::elapsed(0));
  std::vector<Subgraph *>::iterator ii;
  std::cout << sgv.size()<< std::endl;
  for(ii=sgv.begin(); ii!=sgv.end(); ii++) {
	std::cout << (*ii)->canstr << std::endl;
	(*ii)->subgraphCensus(g);
	if ((*ii)->frequency > 0)
		sg_original.setString((*ii)->canstr, (*ii)->frequency);
  }
  //printDenseSubgraphs();
  computeScore();
  //buildHeap();
  sort(graphv.begin(), graphv.end(), greaterG);
  //printDenseSubgraphs();
  //prepareEdgeList();
  mergeGraph();
  printComplexes();
}
bool greaterG(Graph* g1, Graph *g2)
{
    return (g1->getScore() < g2->getScore());
}

void prepareEdgeList() {
  std::vector<Graph *>::iterator ii;
  for(ii=graphv.begin(); ii!=graphv.end(); ii++) {
	(*ii)->makeEdgeList();
  }
}

void printDenseSubgraphs() {
  int i;
  std::vector<Graph *>::iterator ii;
  std::cout << graphv.size()<< std::endl;
  for(ii=graphv.begin(); ii!=graphv.end(); ii++) {
	std::cout << (*ii)->numberEdges() << std::endl;
	std::cout << (*ii)->getScore() << std::endl;
	for (i=0; i<(*ii)->numberNodes(); i++) {
        std::cout << (*ii)->maps(i) << "  ";
	}
	std::cout << std::endl;
  }
}

void printComplexes() {
  int i;
  std::ofstream outfile;
  outfile.open("complexes.txt");
  std::vector<Graph *>::iterator ii;
  std::cout << complexes.size()<< std::endl;
  for(ii=complexes.begin(); ii!=complexes.end(); ii++) {
	std::cout << (*ii)->numberEdges() << std::endl;
	std::cout << (*ii)->getScore() << std::endl;
	for (i=0; i<(*ii)->numberNodes(); i++) {
        std::cout << (*ii)->maps(i) << "  ";
        outfile << (*ii)->maps(i) << " ";
	}
	std::cout << std::endl;
	outfile << std::endl;
  }
}
void computeScore() {
    int i;
    std::vector<Graph *>::iterator ii;
    for(ii=graphv.begin(); ii!=graphv.end(); ii++) {
        (*ii)->computeScore();
        //std::cout << (*ii)->getScore() << "  ";
	}
	//std::cout << std::endl;
}

void mergeGraph() {
    int i,k,j,flag1,flag;
    std::vector<Graph *>::iterator ii;
    std::vector<Graph *>::iterator jj;
    std::vector<Graph *>::iterator kk;
    int counts=0;
    int counts2=0;
    while (graphv.size()>1) {
        flag=0;
        ii=graphv.begin();
        jj=ii;
        jj++;

        j=0;
        for(; jj!=graphv.end(); jj++,j++) {
            if(checkOverlap(*ii,*jj)==true) {
                flag=1;
                break;
            }
        }
        //std::cout << "checkOverlap count = " << j << "\n";
        if (flag==1) {
            std::cout << "overlap count = " << counts << "\n";
            //std::cout << "graphv size= " << graphv.size() << "\n";
            counts++;
            Graph *g1=*ii, *g2=*jj;
            Graph *g3 = new Graph();
            graphv.erase(jj);
            graphv.erase(ii);
            club(g1,g2,g3);
            (*g3).computeScore();
            //std::cout << "g3_score= " << (*g3).getScore() << "\n";
            //std::cout << "graphv size= " << graphv.size() << "\n";
            i=0;
            flag1=0;
            for(kk=graphv.begin(); kk!=graphv.end(); kk++,i++) {
                //std::cout << (*ii)->getScore() << "  ";
                if((*g3).getScore() > (*ii)->getScore()) {
                    graphv.insert(graphv.begin() + i, g3);
                    flag1=1;
                    break;
                }
            }
            if(flag1==0) {
                graphv.insert(graphv.end(), g3);
                //std::cout << "insert end" << "\n";
            }
            //std::cout << "graphv size= " << graphv.size() << "\n";
            //std::cout<<"\n";
        }
        else {
            std::cout << "non-overlap count = " << counts2 << "\n";
            counts2++;
            ii=graphv.begin();
            complexes.push_back(*ii);
            graphv.erase(ii);
        }
    }
    ii=graphv.begin();
    complexes.push_back(*ii);
}

bool checkOverlap(Graph* g1, Graph *g2) {
    int flag=0, i, j, n1, n2, counts;
    double intscore;
#if 0
    std::vector<std::pair<int, int>> g1edge = (*g1).getEdgeList();
    std::vector<std::pair<int, int>> g2edge = (*g2).getEdgeList();
    //std::cout << "listsize2=" << g2edge.size()<< std::endl;
    std::vector<std::pair<int, int>>::iterator ii;
    std::vector<std::pair<int, int>>::iterator jj;
    /*for(ii=g1edge.begin(); ii!=g1edge.end(); ii++) {
        std::cout << (*g1).maps((*ii).first) << "  " << (*g1).maps((*ii).second) << "\n";
	}
	for(ii=g2edge.begin(); ii!=g2edge.end(); ii++) {
        std::cout << (*g2).maps((*ii).first) << "  " << (*g2).maps((*ii).second) << "\n";
	}*/
    for(ii=g1edge.begin(); flag==0,ii!=g1edge.end(); ii++) {
        for(jj=g2edge.begin(); jj!=g2edge.end(); jj++) {
            if ((*g1).maps((*ii).first)==(*g2).maps((*ii).first) &&
                (*g1).maps((*ii).second)==(*g2).maps((*ii).second)) {
                    flag=1;
            }
        }
    }
#endif
    n1=(*g1).numberNodes();
    n2=(*g2).numberNodes();
    counts=0;
    for (i=0; i < n1; i++) {
        for (j=0; j < n2; j++) {
            if((*g1).maps(i)==(*g2).maps(j)) counts++;
        }
    }
    if (counts>1) {
        int c1=0, c2=0;
        double is1=0, is2=0;
        n1=(*g1).numberNodes();
        n2=(*g2).numberNodes();
        /*for (i=0; i < n1; i++) {
            std::cout<< (*g1).maps(i) << " ";
        }
        std::cout<<"\n";
        for (i=0; i < n2; i++) {
            std::cout<< (*g2).maps(i) << " ";
        }
        std::cout<<"\n";*/
        for (i=0; i < n1; i++) {
            int p=0;
            for (j=0; j < n2; j++) {
                if((*g1).maps(i)==(*g2).maps(j)) {
                    p=1;
                    break;
                }
            }
            if (p==1) continue;
            c1++;
            for (j=0; j < n2; j++) {
                is1= is1+intadjustcd(g1,i,g2,j);
            }
        }
        //std::cout<< "is1= " << is1 << "\n";
        for (j=0; j < n2; j++) {
            int p=0;
            for (i=0; i < n1; i++) {
                if((*g1).maps(i)==(*g2).maps(j)) {
                    p=1;
                    break;
                }
            }
            if (p==1) continue;
            c2++;
            for (i=0; i < n1; i++) {
                is2= is2+intadjustcd(g1,i,g2,j);
            }
        }
        //std::cout<< "is2= " << is2 << "\n";
        //std::cout<< "c1= " << c1 << " c2= " << c2 << " n1= " << n1 << " n2= " << n2 << "\n";
        if (c1>0 && c2>0) {
            intscore=sqrt(is1*is2/(c1*n2*c2*n1));
        }
        else if(c1==0) {
            intscore=sqrt(is2/(c2*n1));
        }
        else if(c2==0) {
            intscore=sqrt(is1/(c1*n2));
        }
        //std::cout<< "intscore= " << intscore << "\n";
        //std::cout<< "overlap_threshold= " << overlap_threshold << "\n";
        if(intscore<overlap_threshold) {
            return true;
        }
    }
    return false;
}

void club(Graph* g1, Graph *g2, Graph *g3) {
    int n1, n2, i, j, k;
    int a,b;
    std::map<int, int> nodes;
    n1=(*g1).numberNodes();
    n2=(*g2).numberNodes();
    int map1[n1+n2];
    for (i=0; i < n1; i++) {
        //std::cout<< (*g1).maps(i) << "  ";
        map1[i]=(*g1).maps(i);
        nodes[(*g1).maps(i)]=i;
    }
    k=i;
    for (j=0; j < n2; j++) {
        for (i=0; i < n1; i++) {
            if((*g2).maps(j)==(*g1).maps(i)) {
                break;
            }
        }
        if (i<n1) continue;
        map1[k]=(*g2).maps(j);
        nodes[(*g2).maps(j)]=k;
        k++;
    }
    (*g3).createGraph(k);
#if 0
    std::vector<std::pair<int, int>> g1edge = (*g1).getEdgeList();
    //std::cout << "listsize=" << g1edge.size()<< std::endl;
    std::vector<std::pair<int, int>> g2edge = (*g2).getEdgeList();
    std::vector<std::pair<int, int>>::iterator ii;
    std::vector<std::pair<int, int>>::iterator jj;
    for(ii=g1edge.begin(); ii!=g1edge.end(); ii++) {
        a=(*g1).maps((*ii).first);
        b=(*g1).maps((*ii).second);
        //std::cout<< "a=" << a << " b=" << b << "\n";
        //std::cout<< "nodes(a)=" << nodes[a] << " nodes(b)=" << nodes[b] << "\n";
        (*g3).addEdge(nodes[a], nodes[b]);
    }
    for(jj=g2edge.begin(); jj!=g2edge.end(); jj++) {
        a=(*g2).maps((*jj).first);
        b=(*g2).maps((*jj).second);
        if ((*g3).hasEdge(nodes[a], nodes[b])) {
            //printf("Repeated connection: %d %d\n", a, b);
        } else {
            (*g3).addEdge(nodes[a], nodes[b]);
        }
    }
#endif
    for (i=0; i < n1-1; i++) {
        for (j=i+1; j < n1; j++) {
            if((*g1).hasEdge(i,j)) {
                a=(*g1).maps(i);
                b=(*g1).maps(j);
                //std::cout<< "a=" << a << " b=" << b << "\n";
                //std::cout<< "nodes(a)=" << nodes[a] << " nodes(b)=" << nodes[b] << "\n";
                (*g3).addEdge(nodes[a], nodes[b]);
            }
        }
    }

    for (i=0; i < n2-1; i++) {
        for (j=i+1; j < n2; j++) {
            if((*g2).hasEdge(i,j)) {
                a=(*g2).maps(i);
                b=(*g2).maps(j);
                //std::cout<< "a=" << a << " b=" << b << "\n";
                //std::cout<< "nodes(a)=" << nodes[a] << " nodes(b)=" << nodes[b] << "\n";
                if ((*g3).hasEdge(nodes[a], nodes[b])) {
                    //printf("Repeated connection: %d %d\n", a, b);
                } else {
                    (*g3).addEdge(nodes[a], nodes[b]);
                }
            }
        }
    }


    (*g3).setMap(map1);
    (*g3).sortNeighbours();
    (*g3).makeNeighboursArray();
    //(*g3).makeEdgeList();
    (*g1).deleteGraph();
    (*g2).deleteGraph();
    //std::vector<std::pair<int, int>> g3edge = (*g3).getEdgeList();
    //std::cout << "g3listsize=" << g3edge.size()<< std::endl;
}


double intadjustcd(Graph* g1, int a, Graph *g2, int b) {
    int i,j;
    int na=(*g1).numberNeighbours(a);
    int nb=(*g2).numberNeighbours(b);
    //std::cout << "na= " << na << " nb= " << nb << "\n";
    int *nalist = (*g1).getNeighboursArray(a);
    int *nblist = (*g2).getNeighboursArray(b);
    /*for (i=0; i< na; i++) {
        std::cout<< nalist[i] << "  ";
    }
    std::cout<< "\n";
    for (i=0; i< nb; i++) {
        std::cout<< nblist[i] << "  ";
    }
    std::cout<< "\n";*/
    int ncommon=0;
    for (i=0; i< na; i++) {
        for (j=0; j< nb; j++) {
            if ((*g1).maps(nalist[i])==(*g2).maps(nblist[j])) {
                ncommon++;
                break;
            }
        }
    }
    //std::cout << "ncommon= " << ncommon << "\n";
    double navga=0.0, navgb=0.0, lambdaa=0, lambdab=0;

    for (i=0; i<(*g1).numberNodes(); i++) {
        navga += (*g1).numberNeighbours(i);
    }
    navga=navga/(*g1).numberNodes();

    for (i=0; i<(*g2).numberNodes(); i++) {
        navgb += (*g2).numberNeighbours(i);
    }
    navgb=navgb/(*g2).numberNodes();

    if (navga-na > 0) lambdaa=navga-na;
    if (navgb-nb > 0) lambdab=navgb-nb;
    //std::cout << "lambdaa= " << lambdaa << " lambdab= " << lambdab << "\n";
    return (2*ncommon)/(na+nb+lambdaa+lambdab);
}

void buildHeap() {
    Graph* maxHeapv[graphv.size()];
    int i, heapSize=0;
    std::vector<Graph *>::iterator ii;
    for(ii=graphv.begin(); ii!=graphv.end(); ii++) {
        insertHeap(maxHeapv, heapSize, *ii);
	}
}
void insertHeap(Graph* maxHeapv[], int heapSize, Graph *gr) {
    maxHeapv[heapSize]=gr;
    int i=heapSize, p;
    while (i>0) {
        p=(i-1)/2;
        if(maxHeapv[p]->getScore() > maxHeapv[i]->getScore()) {
            Graph *t;
            t=maxHeapv[p];
            maxHeapv[p]=maxHeapv[i];
            maxHeapv[i]=t;
        }
    }
    heapSize++;
}

// Populate subgraph vector with subgraphs of 'size' read from file 'subgraphs_file'
void read_subgraphs() {
  char buf[MAX_BUF];
  FILE *f = fopen(subgraphs_file, "r");
  if (!f) {
    std::cout << "invalid file" << std::endl;
    exit(1);
  }
  while (fscanf(f, "%s", buf)==1) {
    motif_size=sqrt(strlen(buf));
    std::cout << "motif_size=" << motif_size << "\n";
    Subgraph *sg = new Subgraph(motif_size);
    sg->canstr=new char[strlen(buf)+1];
    int map2[motif_size];
    GraphIsomor::canonicalOrderWithMap(buf, sg->canstr, motif_size, map2);
    sgv.push_back(sg);
  }
  fclose(f);
}

// Compute random networks and result
void compute_results() {
    int i, j;
    std::map<std::string, int>::iterator ii;
    // Create map and init results
    std::map<std::string, int> realMap;
    int resultSize;
    if (strcmp(method, "modet") == 0) {
        realET->inducedFreq(motif_size, g, &realSG, threshold);
        realSG.populateMap(&realMap, motif_size);
    } else if (strcmp(method, "mdet") == 0) {
        realDET->inducedFreq(motif_size, g, &realDSG, threshold);
        realDSG.populateMap(&realMap, motif_size);
    }
    ResultType result[realMap.size()];
    for (ii=realMap.begin(), i=0; ii != realMap.end(); ii++) {
        if (ii->second > 0 && ii->second >= threshold) {
            result[i].s = strdup((ii->first).c_str());
            result[i].freq = ii->second;
            result[i].z_score = 0;
            result[i].avgf_rand = 0;
            result[i].devf_rand = 0;
            i++;
        }
    }
    resultSize = i;

    // Do we have random networks to compute?
    if (rand_number > 0) {
        std::map<std::string, int> randMap[rand_number];
        // Generate all random networks
        printf("Computing random networks: ");
        for (i=0; i < rand_number; i++) {
            std::cout << "random count = " << i+1 << std::endl;
            // Create new random network from previous one
            int nswap = 0;
            Graph::randomGraph(g, num_exchanges, num_tries, &nswap);
            g->sortNeighbours();
            g->makeNeighboursArray();
            // cout << "after:" << nswap << " number of edge swappings\n";
            SubgraphTree randSG;
            SubgraphTree randDSG;
            SubgraphTree randIsomorSG;
            if (strcmp(method, "modet") == 0) {
                randET = new ExpansionTree();
                randET->create(motif_size, &randIsomorSG);
                randET->census(motif_size, g, &randSG);
                randET->inducedFreq(motif_size, g, &randSG, threshold);
                randSG.populateMap(&randMap[i], motif_size);
                delete randET;
            } else if (strcmp(method, "mdet") == 0) {
                randDET = new DET();
                randDET->census(motif_size, g, &randDSG,
                                &randIsomorSG, threshold);
                randDET->inducedFreq(motif_size, g, &randDSG, threshold);
                randDSG.populateMap(&randMap[i], motif_size);
                delete randDET;
            }
        }
        // Compute significance
        for (i=0; i < resultSize; i++) {
            // Average frequency
            double avg = 0;
            for (j=0; j < rand_number; j++)
                avg += randMap[j][result[i].s];
            avg /= rand_number;
            // Standard deviation
            double dev = 0;
            for (j=0; j < rand_number; j++)
                dev += (randMap[j][result[i].s]-avg)*
                       (randMap[j][result[i].s]-avg)/
                       static_cast<double>(rand_number-1);
            dev = sqrt(dev);
            double zscore;
            if (dev != 0)
                zscore = (result[i].freq - avg)/dev;
            else
                zscore = 0;
            result[i].avgf_rand = avg;
            result[i].devf_rand = dev;
            result[i].z_score = fabs(zscore);
            // result[i].z_score = zscore;
        }
    }
    // Sort results
    qsort(result, resultSize, sizeof(ResultType), compare_results);
    FILE *fp;
    fp = fopen("result.txt", "w");
    for (i=0; i < resultSize; i++) {
        if (rand_number > 0) {
            fprintf(fp, "%s : %d : %f\n",
                    result[i].s, result[i].freq, result[i].z_score);
            std::cout << result[i].s << " : " << result[i].freq << " : ";
            std::cout << result[i].z_score << std::endl;
        } else {
            fprintf(fp, "%s : %d\n", result[i].s, result[i].freq);
            std::cout << result[i].s << " : " << result[i].freq << std::endl;
        }
    }
    fprintf(fp, "\nnumber of motifs = %d\n", resultSize);
    std::cout << "\nnumber of motifs = " << resultSize << std::endl;
    fclose(fp);
}

// Compare two different motif results (for sorting)
int compare_results(const void *a, const void *b) {
    ResultType *r1 = (ResultType *)a;
    ResultType *r2 = (ResultType *)b;

    if (r1->z_score > r2->z_score) return -1;
    if (r1->z_score < r2->z_score) return +1;
    if (r1->freq > r2->freq) return -1;
    if (r1->freq < r2->freq) return +1;
    return 0;
}
