#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <set>
#include <queue>
#include <stack>
#include <limits>

class DSU {
private:
    std::vector<int> parent;
    std::vector<int> rank;

public:
    DSU(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unite(int x, int y) {
        int px = find(x);
        int py = find(y);
        if (px == py) return;

        if (rank[px] < rank[py]) {
            parent[px] = py;
        }
        else if (rank[px] > rank[py]) {
            parent[py] = px;
        }
        else {
            parent[py] = px;
            rank[px]++;
        }
    }
};



class Graph {
public:

    enum Representation {
        ADJACENCY_MATRIX,
        ADJACENCY_LIST,
        EDGE_LIST 
    };


    enum GraphType {
        UNWEIGHTED,
        WEIGHTED
    };


    Graph() : representation_(ADJACENCY_MATRIX), graphType_(UNWEIGHTED), numVertices_(0), numEdges_(0) {}


    void readGraph(const std::string& fileName) {
        std::ifstream inputFile(fileName);
        if (!inputFile.is_open()) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            return;
        }

        char representationChar;
        inputFile >> representationChar;

        if (representationChar == 'C') {
            representation_ = ADJACENCY_MATRIX;
        }
        else if (representationChar == 'L') {
            representation_ = ADJACENCY_LIST;
        }
        else if (representationChar == 'E') {
            representation_ = EDGE_LIST;
        }
        else {
            std::cerr << "Invalid representation character in file." << std::endl;
            inputFile.close();
            return;
        }

        if (representation_ == ADJACENCY_MATRIX) {
            int n;
            inputFile >> n;
            numVertices_ = n;
            int directed, weighted;
            inputFile >> directed >> weighted;
            graphType_ = (weighted == 1) ? WEIGHTED : UNWEIGHTED;

            adjacencyMatrix_.resize(n, std::vector<int>(n, 0));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    inputFile >> adjacencyMatrix_[i][j];
                    if (adjacencyMatrix_[i][j] != 0) {
                        numEdges_++;
                    }
                }
            }
            if (directed == 0) {
                numEdges_ /= 2;
            }

        }

        else if (representation_ == ADJACENCY_LIST) {
            int n;
            inputFile >> n;
            numVertices_ = n;
            int directed, weighted;
            inputFile >> directed >> weighted;

            graphType_ = (weighted == 1) ? WEIGHTED : UNWEIGHTED;

            adjacencyList_.resize(n);
            for (int i = 0; i < n; ++i) {
                std::string line;
                std::getline(inputFile >> std::ws, line);
                std::stringstream ss(line);
                int neighbor, weight;
                if (graphType_ == WEIGHTED) {
                    while (ss >> neighbor >> weight) {
                        adjacencyList_[i].push_back({ neighbor - 1, weight });
                        numEdges_++;
                    }
                }
                else {
                    while (ss >> neighbor) {
                        adjacencyList_[i].push_back({ neighbor - 1, 1 });
                        numEdges_++;
                    }
                }
            }
            if (directed == 0) {
                numEdges_ /= 2;
            }
        }

        else if (representation_ == EDGE_LIST) {
            int n, m;
            inputFile >> n >> m;
            numVertices_ = n;
            numEdges_ = m;
            int directed, weighted;
            inputFile >> directed >> weighted;

            graphType_ = (weighted == 1) ? WEIGHTED : UNWEIGHTED;


            edgeList_.resize(m);
            for (int i = 0; i < m; ++i) {
                int u, v, weight = 1;
                if (graphType_ == WEIGHTED) {
                    inputFile >> u >> v >> weight;
                }
                else {
                    inputFile >> u >> v;
                }
                edgeList_[i] = { u - 1, v - 1, weight };
            }
        }
        inputFile.close();
    }



    int checkEuler(bool& circleExist) {
        Representation originalRep = representation_;
        transformToAdjList();
        std::vector<int> degrees(numVertices_, 0);
        int oddDegreeCount = 0;
        int startVertex = -1;
        for (int i = 0; i < numVertices_; ++i) {
            degrees[i] = adjacencyList_[i].size();
            if (degrees[i] % 2 == 1) {
                oddDegreeCount++;
                if (startVertex == -1) startVertex = i;
            }
            if (degrees[i] > 0 && startVertex == -1) {
                startVertex = i;
            }
        }
        DSU dsu(numVertices_);
        for (int u = 0; u < numVertices_; ++u) {
            for (const auto& edge : adjacencyList_[u]) {
                int v = edge.first;
                if (u < v) {
                    dsu.unite(u, v);
                }
            }
        }
        std::set<int> components;
        for (int i = 0; i < numVertices_; ++i) {
            if (degrees[i] > 0) {
                components.insert(dsu.find(i));
            }
        }
        if (originalRep == ADJACENCY_MATRIX) transformToAdjMatrix();
        else if (originalRep == EDGE_LIST) transformToListOfEdges();
        if (components.size() > 1) return 0;
        if (oddDegreeCount == 0) {
            circleExist = true;
            return startVertex + 1;
        }
        else if (oddDegreeCount == 2) {
            circleExist = false;
            return startVertex + 1;
        }
        return 0;
    }

    std::vector<int> getEuleranTourFleri() {
        std::vector<int> result;
        bool circleExist;
        int start = checkEuler(circleExist);
        if (start == 0) return result;

        Graph tempGraph = *this;
        tempGraph.transformToAdjList();
        std::stack<int> currPath;
        start--;
        currPath.push(start);

        while (!currPath.empty()) {
            int u = currPath.top();
            int nextV = -1;

            for (auto& edge : tempGraph.adjacencyList_[u]) {
                if (edge.second > 0) {
                    edge.second--;
                    for (auto& rev_edge : tempGraph.adjacencyList_[edge.first]) {
                        if (rev_edge.first == u && rev_edge.second > 0) {
                            rev_edge.second--;
                            break;
                        }
                    }

                    bool connected = tempGraph.isConnected();

                    edge.second++;
                    for (auto& rev_edge : tempGraph.adjacencyList_[edge.first]) {
                        if (rev_edge.first == u) {
                            rev_edge.second++;
                            break;
                        }
                    }

                    if (connected || nextV == -1) {
                        nextV = edge.first;
                        if (connected) break;
                    }
                }
            }

            if (nextV != -1) {
                for (auto& edge : tempGraph.adjacencyList_[u]) {
                    if (edge.first == nextV && edge.second > 0) {
                        edge.second--;
                        break;
                    }
                }
                for (auto& edge : tempGraph.adjacencyList_[nextV]) {
                    if (edge.first == u && edge.second > 0) {
                        edge.second--;
                        break;
                    }
                }

                currPath.push(nextV);
            }
            else {
                result.push_back(u + 1);
                currPath.pop();
            }
        }

        std::reverse(result.begin(), result.end());
        return result;
    }

    std::vector<int> getEuleranTourEffective() {
        std::vector<int> result;
        bool circleExist;
        int start = checkEuler(circleExist);
        if (start == 0) return result;

        std::vector<std::vector<int>> tempAdj;
        for (const auto& edges : adjacencyList_) {
            tempAdj.emplace_back();
            for (const auto& edge : edges) {
                tempAdj.back().push_back(edge.first);
            }
        }

        std::stack<int> currPath;
        start--;
        currPath.push(start);

        while (!currPath.empty()) {
            int u = currPath.top();

            if (!tempAdj[u].empty()) {
                int v = tempAdj[u].back();
                tempAdj[u].pop_back();

                for (auto it = tempAdj[v].begin(); it != tempAdj[v].end(); ++it) {
                    if (*it == u) {
                        tempAdj[v].erase(it);
                        break;
                    }
                }

                currPath.push(v);
            }
            else {
                result.push_back(u + 1);
                currPath.pop();
            }
        }

        std::reverse(result.begin(), result.end());
        return result;
    }
    

    void addEdge(int from, int to, int weight = 1) {

        --from;
        --to;

        if (from < 0 || from >= numVertices_ || to < 0 || to >= numVertices_) {
            std::cerr << "Invalid vertex index." << std::endl;
            return;
        }

        if (representation_ == ADJACENCY_MATRIX) {
            if (graphType_ == WEIGHTED) {
                adjacencyMatrix_[from][to] = weight;
                if (from != to) {
                    if (adjacencyMatrix_[to][from] != 0) {
                        numEdges_++;
                    }
                }
                else {
                    numEdges_++;
                }

            }
            else {
                if (adjacencyMatrix_[from][to] == 0) {
                    adjacencyMatrix_[from][to] = 1;
                    if (from != to) {
                        if (adjacencyMatrix_[to][from] != 0) {
                            numEdges_++;
                        }
                    }
                    else {
                        numEdges_++;
                    }
                }
            }

        }

        else if (representation_ == ADJACENCY_LIST) {
            bool edgeExists = false;
            for (auto& edge : adjacencyList_[from]) {
                if (edge.first == to) {
                    edgeExists = true;
                    break;
                }
            }
            if (!edgeExists) {
                adjacencyList_[from].push_back({ to, weight });
                numEdges_++;
            }
        }

        else if (representation_ == EDGE_LIST) {
            edgeList_.push_back({ from, to, weight });
            numEdges_++;
        }
    }


    void removeEdge(int from, int to) {
        --from;
        --to;

        if (from < 0 || from >= numVertices_ || to < 0 || to >= numVertices_) {
            std::cerr << "Invalid vertex index." << std::endl;
            return;
        }


        if (representation_ == ADJACENCY_MATRIX) {
            if (adjacencyMatrix_[from][to] != 0) {
                adjacencyMatrix_[from][to] = 0;
                numEdges_--;
            }

        }

        else if (representation_ == ADJACENCY_LIST) {
            for (auto it = adjacencyList_[from].begin(); it != adjacencyList_[from].end(); ++it) {
                if (it->first == to) {
                    adjacencyList_[from].erase(it);
                    numEdges_--;
                    break;
                }
            }
        }

        else if (representation_ == EDGE_LIST) {
            for (auto it = edgeList_.begin(); it != edgeList_.end(); ++it) {
                if (it->from == from && it->to == to) {
                    edgeList_.erase(it);
                    numEdges_--;
                    break;
                }
            }
        }
    }


    int changeEdge(int from, int to, int newWeight) {

        --from;
        --to;

        if (from < 0 || from >= numVertices_ || to < 0 || to >= numVertices_) {
            std::cerr << "Invalid vertex index." << std::endl;
            return -1;
        }

        int oldWeight = -1;

        if (representation_ == ADJACENCY_MATRIX) {
            oldWeight = adjacencyMatrix_[from][to];
            adjacencyMatrix_[from][to] = newWeight;
        }
        else if (representation_ == ADJACENCY_LIST) {
            for (auto& edge : adjacencyList_[from]) {
                if (edge.first == to) {
                    oldWeight = edge.second;
                    edge.second = newWeight;
                    break;
                }
            }
        }
        else if (representation_ == EDGE_LIST) {
            for (auto& edge : edgeList_) {
                if (edge.from == from && edge.to == to) {
                    oldWeight = edge.weight;
                    edge.weight = newWeight;
                    break;
                }
            }
        }
        return oldWeight;
    }

    void transformToAdjList() {
        if (representation_ == ADJACENCY_LIST) {
            return;
        }

        adjacencyList_.clear();
        adjacencyList_.resize(numVertices_);
        graphType_ = (graphType_ == WEIGHTED) ? WEIGHTED : UNWEIGHTED;

        if (representation_ == ADJACENCY_MATRIX) {
            for (int i = 0; i < numVertices_; ++i) {
                for (int j = 0; j < numVertices_; ++j) {
                    if (adjacencyMatrix_[i][j] != 0) {
                        adjacencyList_[i].push_back({ j, adjacencyMatrix_[i][j] });
                    }
                }
            }
        }
        else if (representation_ == EDGE_LIST) {
            for (const auto& edge : edgeList_) {
                adjacencyList_[edge.from].push_back({ edge.to, edge.weight });
            }
        }
        representation_ = ADJACENCY_LIST;
    }


    void transformToAdjMatrix() {
        if (representation_ == ADJACENCY_MATRIX) {
            return;
        }

        adjacencyMatrix_.clear();
        adjacencyMatrix_.resize(numVertices_, std::vector<int>(numVertices_, 0));

        if (representation_ == ADJACENCY_LIST) {
            for (int i = 0; i < numVertices_; ++i) {
                for (const auto& edge : adjacencyList_[i]) {
                    adjacencyMatrix_[i][edge.first] = edge.second;
                }
            }
        }
        else if (representation_ == EDGE_LIST) {
            for (const auto& edge : edgeList_) {
                adjacencyMatrix_[edge.from][edge.to] = edge.weight;
            }
        }
        representation_ = ADJACENCY_MATRIX;
    }


    void transformToListOfEdges() {
        if (representation_ == EDGE_LIST) {
            return;
        }

        edgeList_.clear();

        if (representation_ == ADJACENCY_MATRIX) {
            for (int i = 0; i < numVertices_; ++i) {
                for (int j = 0; j < numVertices_; ++j) {
                    if (adjacencyMatrix_[i][j] != 0) {
                        edgeList_.push_back({ i, j, adjacencyMatrix_[i][j] });
                    }
                }
            }
        }
        else if (representation_ == ADJACENCY_LIST) {
            for (int i = 0; i < numVertices_; ++i) {
                for (const auto& edge : adjacencyList_[i]) {
                    edgeList_.push_back({ i, edge.first, edge.second });
                }
            }
        }
        representation_ = EDGE_LIST;
    }


    void writeGraph(const std::string& fileName) {
        std::ofstream outputFile(fileName);
        if (!outputFile.is_open()) {
            std::cerr << "Error opening file for writing: " << fileName << std::endl;
            return;
        }

        if (representation_ == ADJACENCY_MATRIX) {
            outputFile << "C " << numVertices_ << std::endl;
            outputFile << "0 ";
            if (graphType_ == WEIGHTED)
                outputFile << "1" << std::endl;
            else
                outputFile << "0" << std::endl;

            for (int i = 0; i < numVertices_; ++i) {
                for (int j = 0; j < numVertices_; ++j) {
                    outputFile << adjacencyMatrix_[i][j] << " ";
                }
                outputFile << std::endl;
            }
        }
        else if (representation_ == ADJACENCY_LIST) {
            outputFile << "L " << numVertices_ << std::endl;
            outputFile << "0 ";
            if (graphType_ == WEIGHTED)
                outputFile << "1" << std::endl;
            else
                outputFile << "0" << std::endl;

            for (int i = 0; i < numVertices_; ++i) {
                for (const auto& edge : adjacencyList_[i]) {
                    outputFile << edge.first + 1 << " ";
                    if (graphType_ == WEIGHTED)
                        outputFile << edge.second << " ";
                }
                outputFile << std::endl;
            }
        }
        else if (representation_ == EDGE_LIST) {
            outputFile << "E " << numVertices_ << " " << numEdges_ << std::endl;
            outputFile << "0 ";
            if (graphType_ == WEIGHTED)
                outputFile << "1" << std::endl;
            else
                outputFile << "0" << std::endl;
            for (const auto& edge : edgeList_) {
                outputFile << edge.from + 1 << " " << edge.to + 1 << " " << edge.weight << std::endl;
            }
        }
        outputFile.close();
    }


private:
    Representation representation_;
    GraphType graphType_;
    int numVertices_;
    int numEdges_;

    std::vector<std::vector<int>> adjacencyMatrix_;
    std::vector<std::vector<std::pair<int, int>>> adjacencyList_;

    struct Edge {
        int from;
        int to;
        int weight;
    };
    std::vector<Edge> edgeList_;

    bool isConnected() {
        if (numVertices_ == 0) return true;

        int start = -1;
        for (int i = 0; i < numVertices_; ++i) {
            if (!adjacencyList_[i].empty()) {
                start = i;
                break;
            }
        }
        if (start == -1) return true;

        std::vector<bool> visited(numVertices_, false);
        std::queue<int> q;
        q.push(start);
        visited[start] = true;

        int count = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            count++;
            for (const auto& edge : adjacencyList_[u]) {
                int v = edge.first;
                if (edge.second > 0 && !visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }

        return count == std::count_if(adjacencyList_.begin(), adjacencyList_.end(),
                                      [](const std::vector<std::pair<int, int>>& edges) { return !edges.empty(); });
    }
};
int main(){

	freopen("output.txt", "w", stdout);

    Graph graph;

	graph.readGraph("input.txt");

	

	bool circleExist = false;

	int v = graph.checkEuler(circleExist);



	if(v == 0) {

		std::cout<<"0\n";

		exit(0);

	}

	else

	{

	  //  std::vector<int> path = graph.getEuleranTourEffective();

	   std::vector<int> path = graph.getEuleranTourFleri();

        std::cout<<path[0]<<"\n";

    	for (int i = 0; i < path.size(); ++i) {

    	    std::cout << path[i] << ((i + 1 < path.size()) ? " ":"\n");

    	}

	}

	

	return 0;

}