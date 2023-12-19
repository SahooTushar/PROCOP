// Copyright © 2019 Sabyasachi Patra, IIIT Bhubaneswar, All rights reserved.
#ifndef GRAPH_H_
#define GRAPH_H_

#include <stdio.h>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>

class Graph {
 private:
     int numNodes;
     int numEdges;
     int *numNeighbours;
     bool **adjMatrix;
     int  **neighboursArray;
     int *map1;
     std::vector<int> *neighbours;
     std::vector<std::pair<int, int>> edges;
     void initGraph();
     double score;

 public:
     Graph();
     ~Graph();
     void deleteGraph();
     int numberNodes() {
         return numNodes;
     }
     int numberEdges() {
         return numEdges;
     }
     bool **adjacencyMatrix() {
         return adjMatrix;
     }
     void createGraph(int n);
     void addEdge(int a, int b);
     void rmEdge(int a, int b);
     bool hasEdge(int a, int b) {
         return adjMatrix[a][b];
     }
     int maps(int a) {
        return map1[a];
     }
     double getScore() {
        return score;
     }
     double adjustcd(int a, int b);
     bool isConnected(int a, int b) {
         return adjMatrix[a][b] || adjMatrix[b][a];
     }
     int numberNeighbours(int a) {
         return numNeighbours[a];
     }
     int *numberNeighbours() {
         return numNeighbours;
     }
     void sortNeighbours();
     void sortNeighboursArray();
     void makeNeighboursArray();
     void makeNeighboursVector();
     void makeEdgeList();
     void setMap(int *);
     void computeScore();
     std::vector<int> *getNeighbours(int a) {
         return &neighbours[a];
     }
     int **getNeighboursArray() {
         return neighboursArray;
     }
     int *getNeighboursArray(int a) {
         return neighboursArray[a];
     }
     std::vector<std::pair<int, int>> getEdgeList() {
         return edges;
     }
     static int compareInt(const void *a, const void *b);
     static void readGraphFile(Graph *g, const char *s);
     static void strToGraph(Graph *g, const char *s, int size);
     static void randomGraph(Graph *g, int num_ex, int num_tr, int *nswap);
     static int getRandom(int a, int b);
};
#endif  // GRAPH_H_
