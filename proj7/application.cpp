// Amaan Shah, UIC, Spring 2021
// application.cpp
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <stack>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>
#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

// Helper Function Prototypes:
void printBuildingStats(BuildingInfo Building);
void printCoordinates(Coordinates c);
bool isVisited(vector<long long> visited, long long currentV);
void printPath(long long startVertex, long long destVertex,
    map<long long, long long> &pred);
Coordinates findNearestNode(map<long long, Coordinates> &Nodes,
    vector<FootwayInfo> &Footways, BuildingInfo &Building);
vector<long long> dijkstraSearch(graph<long long, double> &G,
    long long startVertex, map<long long, double> &distances,
    map<long long, long long> &pred);
bool locateBuilding(vector<BuildingInfo> &Buildings,
    string name, BuildingInfo &node);

// Functor for c++ Priority Queue:
class prioritize {
 public:
  	bool operator()(pair<long long, double> const& v1,
  	pair<long long, double> const& v2) {
  	  if (v1.second == v2.second) {
  	    return v1.first > v2.first;
  	  }
      return v1.second > v2.second;
    }
};


int main() {
  const double INF = numeric_limits<double>::max();
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  graph <long long, double> G;
  // insert verticies into the graph
  for (auto &x : Nodes) {
      G.addVertex(x.first);
  }
  // add edges between nodes
  for (auto &x : Footways) {
      int size = x.Nodes.size();
      for (int i = 0; i < size - 1; i++) {
          long long src = x.Nodes[i];
          long long dest = x.Nodes[i + 1];
          double lat1 = Nodes.at(src).Lat;
          double lon1 = Nodes.at(src).Lon;
          double lat2 = Nodes.at(dest).Lat;
          double lon2 = Nodes.at(dest).Lon;
          // find distance between src and dest
          double dist = distBetween2Points(lat1, lon1, lat2, lon2);
          // add edge from src to dest
          G.addEdge(src, dest, dist);
          // add edge from dest to src
          G.addEdge(dest, src, dist);
      }
  }
  // output graph stats
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;
  //
  // Navigation from building to building
  //
  string startBuilding, destBuilding;

  while (true) {
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);
    if (startBuilding == "#") {
        break;
    }
    cout << "Enter destination (partial name or abbreviation)> ";
    getline(cin, destBuilding);
    // Lookup buildings:
    BuildingInfo start;
    BuildingInfo dest;
    bool foundStart = locateBuilding(Buildings, startBuilding, start);
    bool foundDest = locateBuilding(Buildings, destBuilding, dest);
    if (!foundStart) {
        cout << "Start building not found" << endl << endl;
        continue;
    } else if (!foundDest) {
        cout << "Destination building not found" << endl << endl;
        continue;
    }
    cout << "Starting point:" << endl << " ";
    printBuildingStats(start);
    cout << "Destination point:" << endl << " ";
    printBuildingStats(dest);
    cout << endl;
    // Find nearest start and dest nodes:
    Coordinates nearestStart = findNearestNode(Nodes, Footways, start);
    Coordinates nearestDest = findNearestNode(Nodes, Footways, dest);
    cout << "Nearest start node:" << endl << " ";
    printCoordinates(nearestStart);
    cout << "Nearest destination node:" << endl << " ";
    printCoordinates(nearestDest);
    cout << endl;
    // run Dijkstra's alg, output distance and path to destination:
    cout << "Navigating with Dijkstra..." << endl;
    long long startVertex = nearestStart.ID;
    long long destVertex = nearestDest.ID;
    map<long long, double> distances;
    map<long long, long long> pred;
    vector<long long> visited = dijkstraSearch(G, startVertex, distances, pred);
    
    if (distances[destVertex] == INF) {
        cout << "Sorry, destination unreachable" << endl << endl;
        continue;
    } else {
        cout << "Distance to dest: " << distances[destVertex];
        cout << " miles" << endl;
        cout << "Path: ";
        printPath(startVertex, destVertex, pred);
    }
  }

  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}

// Locates a building from the graph.
// If found returns true and info about the building.
// else returns false;
bool locateBuilding(vector<BuildingInfo> &Buildings,
string name, BuildingInfo &node) {
    for (auto &x : Buildings) {
        if (name == x.Abbrev) {
            node = x;
            return true;
        } else if (x.Fullname.find(name) != string::npos) {
            node = x;
            return true;
        }
    }
    return false;
}

// Used to print building stats
void printBuildingStats(BuildingInfo b) {
    cout << b.Fullname << endl;
    cout << " (" << b.Coords.Lat << ", " << b.Coords.Lon << ")" << endl;
}

// Finds the closest node to a given building
Coordinates findNearestNode(map<long long, Coordinates> &Nodes,
vector<FootwayInfo> &Footways, BuildingInfo &Building) {
    double bLat = Building.Coords.Lat;
    double bLon = Building.Coords.Lon;
    long long closestID = Footways[0].Nodes[0];
    double cLat = Nodes[closestID].Lat;
    double cLon = Nodes[closestID].Lon;
    double minDistance = distBetween2Points(bLat, bLon, cLat, cLon);
    
    for (auto &x : Footways) {
        for (auto &y : x.Nodes) {
            double newLat = Nodes[y].Lat;
            double newLon = Nodes[y].Lon;
            double newDist = distBetween2Points(bLat, bLon, newLat, newLon);
            if (minDistance > newDist) {
                closestID = y;
                cLat = newLat;
                cLon = newLon;
                minDistance = newDist;
            }
        }
    }
    Coordinates node(closestID, cLat, cLon);
    return node;
}

// Prints the stats of the nearest nodes
void printCoordinates(Coordinates c) {
    cout << c.ID << endl << " ";
    cout << "(" << c.Lat << ", " << c.Lon << ")" << endl;
}

// Checks if we have already visited a vertex
bool isVisited(vector<long long> visited, long long currentV) {
  int size = visited.size();
	for (int i = 0; i < size; i++) {
  	if (visited[i] == currentV)
      return true;
  }
  return false;
}

// Perfoms Dijkstra's Alg to find the optimal path
vector<long long> dijkstraSearch(graph<long long, double> &G,
long long startVertex, map<long long, double> &distances, 
map<long long, long long> &pred) {
    const double INF = numeric_limits<double>::max();
    vector<long long> visited;
    priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> unvisitedQueue;
    vector<long long> vertexVect = G.getVertices();
    int size = vertexVect.size();
    for (int i = 0; i < size; i++) {
        distances[vertexVect[i]] = INF;
        pred[vertexVect[i]] = 0;
        unvisitedQueue.push(make_pair(vertexVect[i], INF));
    }
    distances[startVertex] = 0;
    unvisitedQueue.push(make_pair(startVertex, 0));
    while (!unvisitedQueue.empty()) {
        long long currentV = unvisitedQueue.top().first;
        unvisitedQueue.pop();
        
        if (distances[currentV] == INF) {
            break;
        } else if (isVisited(visited, currentV)) {
            continue;
        } else {
            visited.push_back(currentV);
        }
        set<long long> neighbors = G.neighbors(currentV);
        for (auto &i : neighbors) {
            double edgeweight = 0.0;
            G.getWeight(currentV, i, edgeweight);
            double altDistance = distances[currentV] + edgeweight;
            if (altDistance < distances[i]) {
                distances[i] = altDistance;
                pred[i] = currentV;
                unvisitedQueue.push(make_pair(i,altDistance));
            }
        }
    }
    return visited;
}

// Uses the predecessor "array" to print path from
// source to destination
void printPath(long long startVertex, long long destVertex,
map<long long, long long> &pred) {
    stack<long long> path;
    path.push(destVertex);
    long long previousID = pred.at(destVertex);
    // adding verticies to stack in reverse order
    while (true) {
        if (previousID == 0) {
            break;
        }
        path.push(previousID);
        previousID = pred[previousID];
    }
    // print path stack verticies in order
    cout << path.top();
    path.pop();
    while(true) {
        if (path.empty()) {
            break;
        }
        cout << "->" << path.top();
        path.pop();
    }
    cout << endl << endl;
}














