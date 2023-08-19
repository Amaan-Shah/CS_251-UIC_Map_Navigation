// Amaan Shah, UIC, Spring 2021
// graph.h
//
// Graph class using adjacency list representation.
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_map>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
    // New implementation of adjacency list using nested maps
    unordered_map <VertexT, unordered_map<VertexT, WeightT>> AdjList;
    unordered_map <VertexT, set<VertexT>> NeighborsList;
    vector <VertexT> Vertices;
    int verts;
    int edges;

 public:
  // constructor:
  //
  graph() {
     verts = 0;
     edges = 0;
  }

  // NumVertices:
  // Total Complexity: O(1)
  //
  // Returns the # of vertices currently in the graph.
  int NumVertices() const {
    return verts;
  }

  // NumEdges:
  // Total Complexity: O(1)
  //
  // Returns the # of edges currently in the graph.
  int NumEdges() const {
    return edges;
  }

  // addVertex:
  // Total Complexity: O(logN)
  //
  // Adds the vertex v to the graph, if vertex is already
  // in the graph returns false, else returns true
  bool addVertex(VertexT v) {
    if (AdjList.count(v) == 1) {  // O(logN)
      return false;
    }
    unordered_map <VertexT, WeightT> List;
    set<VertexT>  S;
    NeighborsList.emplace(v, S);  // O(1)
    AdjList.emplace(v, List);  // O(1)
    Vertices.push_back(v);  // O(1)
    verts++;
    return true;
  }

  // addEdge:
  // Total Complexity: O(logN)
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist, false is returned.
  // If the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    if (AdjList.count(from) == 0 || AdjList.count(to) == 0) {  // O(logN)
        return false;
    }
    // insert / update the edge:
    if (AdjList.at(from).count(to) == 0) {  //
        AdjList.at(from).emplace(to, weight);  // 2 * O(1)
        NeighborsList.at(from).emplace(to);  // 2 * O(1)
        edges++;
    } else {
        AdjList.at(from).at(to) = weight;  // 2 * O(1)
    }
    return true;
  }

  // getWeight:
  // Total Complexity: // O(logN)
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    if (AdjList.count(from) == 0) {  // O(logN)
        return false;
    } else if (AdjList.count(to) == 0) {  // O(logN)
        return false;
    }
    if (AdjList.at(from).count(to) != 1) {  // O(1) + O(logN)
      return false;
    }
    weight = AdjList.at(from).at(to);  // 2 * O(1)
    return true;
  }

  // neighbors:
  // Total Complexity: O(logN)
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT>  S;
    if (AdjList.count(v) == 0) {  // O(logN)
        return S;
    }
    return NeighborsList.at(v);  // O(1)
  }

  // getVertices:
  // Total Complexity: O(1)
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  vector<VertexT> getVertices() const {
    return Vertices;
  }

  // dump:
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    A: (A,B,80) (A,C,100), …
  //    B: (B,A,100) (B,F,123), …
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    for (int i = 0; i < this->NumVertices(); ++i) {
      output << " " << i << ". " << this->Vertices[i] << endl;
    }

    output << endl;
    output << "**Edges:" << endl;
    for (auto &x : AdjList) {
        output << x.first << ":";
        for (auto &y : AdjList.at(x.first)) {
            output << " (" << x.first << "," << y.first << "," << y.second << ")";
        }
        output << endl;
    }
    output << "**************************************************" << endl;
  }
};
