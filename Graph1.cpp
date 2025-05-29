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
    explicit DSU(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        while (x != parent[x]) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX == rootY) return;

        if (rank[rootX] < rank[rootY]) {
            parent[rootX] = rootY;
        }
        else {
            parent[rootY] = rootX;
            if (rank[rootX] == rank[rootY]) {
                ++rank[rootX];
            }
        }
    }
};

class Graph {
public:
    enum Representation { ADJACENCY_MATRIX, ADJACENCY_LIST, EDGE_LIST };
    enum GraphType { UNWEIGHTED, WEIGHTED };

    Graph()
        : representation_(ADJACENCY_MATRIX),
          graphType_(UNWEIGHTED),
          numVertices_(0),
          numEdges_(0) {}

    void readGraph(const std::string& fileName) {
        std::ifstream in(fileName);
        if (!in.is_open()) {
            std::cerr << "Could not open file: " << fileName << "\n";
            return;
        }

        char repChar{};
        in >> repChar;
        if (repChar == 'C') representation_ = ADJACENCY_MATRIX;
        else if (repChar == 'L') representation_ = ADJACENCY_LIST;
        else if (repChar == 'E') representation_ = EDGE_LIST;
        else {
            std::cerr << "Unknown graph representation symbol.\n";
            in.close();
            return;
        }

        if (representation_ == ADJACENCY_MATRIX) {
            loadAdjacencyMatrix(in);
        }
        else if (representation_ == ADJACENCY_LIST) {
            loadAdjacencyList(in);
        }
        else {
            loadEdgeList(in);
        }

        in.close();
    }

    int checkEuler(bool& circleExist) {
        auto originalRep = representation_;
        transformToAdjList();

        std::vector<int> degrees(numVertices_, 0);
        int oddCount = 0;
        int startVertex = -1;

        for (int i = 0; i < numVertices_; ++i) {
            degrees[i] = static_cast<int>(adjacencyList_[i].size());
            if (degrees[i] % 2 != 0) {
                oddCount++;
                if (startVertex == -1) startVertex = i;
            }
            else if (degrees[i] > 0 && startVertex == -1) {
                startVertex = i;
            }
        }

        DSU dsu(numVertices_);
        for (int u = 0; u < numVertices_; ++u) {
            for (auto const& e : adjacencyList_[u]) {
                int v = e.first;
                if (u < v) dsu.unite(u, v);
            }
        }

        std::set<int> compSet;
        for (int i = 0; i < numVertices_; ++i) {
            if (degrees[i] > 0) compSet.insert(dsu.find(i));
        }

        if (originalRep == ADJACENCY_MATRIX) transformToAdjMatrix();
        else if (originalRep == EDGE_LIST) transformToListOfEdges();

        if (compSet.size() > 1) return 0;

        if (oddCount == 0) {
            circleExist = true;
            return startVertex + 1;
        }
        if (oddCount == 2) {
            circleExist = false;
            return startVertex + 1;
        }
        return 0;
    }

    std::vector<int> getEuleranTourFleri() {
        std::vector<int> path;
        bool circleExist = false;
        int start = checkEuler(circleExist);
        if (start == 0) return path;

        Graph copyGraph = *this;
        copyGraph.transformToAdjList();

        std::stack<int> stack;
        --start;
        stack.push(start);

        while (!stack.empty()) {
            int u = stack.top();
            int nextVertex = -1;

            for (auto& e : copyGraph.adjacencyList_[u]) {
                if (e.second > 0) {
                    e.second--;
                    for (auto& revE : copyGraph.adjacencyList_[e.first]) {
                        if (revE.first == u && revE.second > 0) {
                            revE.second--;
                            break;
                        }
                    }

                    bool isConn = copyGraph.isConnected();

                    e.second++;
                    for (auto& revE : copyGraph.adjacencyList_[e.first]) {
                        if (revE.first == u) {
                            revE.second++;
                            break;
                        }
                    }

                    if (isConn || nextVertex == -1) {
                        nextVertex = e.first;
                        if (isConn) break;
                    }
                }
            }

            if (nextVertex != -1) {
                for (auto& e : copyGraph.adjacencyList_[u]) {
                    if (e.first == nextVertex && e.second > 0) {
                        e.second--;
                        break;
                    }
                }
                for (auto& e : copyGraph.adjacencyList_[nextVertex]) {
                    if (e.first == u && e.second > 0) {
                        e.second--;
                        break;
                    }
                }
                stack.push(nextVertex);
            }
            else {
                path.push_back(u + 1);
                stack.pop();
            }
        }

        std::reverse(path.begin(), path.end());
        return path;
    }

    std::vector<int> getEuleranTourEffective() {
        std::vector<int> path;
        bool circleExist = false;
        int start = checkEuler(circleExist);
        if (start == 0) return path;

        std::vector<std::vector<int>> tmpAdj(numVertices_);
        for (int i = 0; i < numVertices_; ++i) {
            for (const auto& e : adjacencyList_[i]) {
                tmpAdj[i].push_back(e.first);
            }
        }

        std::stack<int> stack;
        --start;
        stack.push(start);

        while (!stack.empty()) {
            int u = stack.top();
            if (!tmpAdj[u].empty()) {
                int v = tmpAdj[u].back();
                tmpAdj[u].pop_back();

                auto& vec = tmpAdj[v];
                auto it = std::find(vec.begin(), vec.end(), u);
                if (it != vec.end()) {
                    vec.erase(it);
                }

                stack.push(v);
            }
            else {
                path.push_back(u + 1);
                stack.pop();
            }
        }
        std::reverse(path.begin(), path.end());
        return path;
    }

    void addEdge(int from, int to, int weight = 1) {
        --from; --to;
        if (from < 0 || from >= numVertices_ || to < 0 || to >= numVertices_) {
            std::cerr << "Invalid vertex index.\n";
            return;
        }
        if (representation_ == ADJACENCY_MATRIX) {
            if (graphType_ == WEIGHTED) {
                adjacencyMatrix_[from][to] = weight;
                if (from != to) {
                    if (adjacencyMatrix_[to][from] != 0) numEdges_++;
                }
                else {
                    numEdges_++;
                }
            }
            else {
                if (adjacencyMatrix_[from][to] == 0) {
                    adjacencyMatrix_[from][to] = 1;
                    if (from != to) {
                        if (adjacencyMatrix_[to][from] != 0) numEdges_++;
                    }
                    else numEdges_++;
                }
            }
        }
        else if (representation_ == ADJACENCY_LIST) {
            bool found = false;
            for (const auto& e : adjacencyList_[from]) {
                if (e.first == to) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                adjacencyList_[from].emplace_back(to, weight);
                numEdges_++;
            }
        }
        else {
            edgeList_.push_back({ from, to, weight });
            ++numEdges_;
        }
    }

    void removeEdge(int from, int to) {
        --from; --to;
        if (from < 0 || from >= numVertices_ || to < 0 || to >= numVertices_) {
            std::cerr << "Invalid vertex index.\n";
            return;
        }
        if (representation_ == ADJACENCY_MATRIX) {
            if (adjacencyMatrix_[from][to] != 0) {
                adjacencyMatrix_[from][to] = 0;
                numEdges_--;
            }
        }
        else if (representation_ == ADJACENCY_LIST) {
            auto& adj = adjacencyList_[from];
            for (auto it = adj.begin(); it != adj.end(); ++it) {
                if (it->first == to) {
                    adj.erase(it);
                    numEdges_--;
                    break;
                }
            }
        }
        else {
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
        --from; --to;
        if (from < 0 || from >= numVertices_ || to < 0 || to >= numVertices_) {
            std::cerr << "Invalid vertex index.\n";
            return -1;
        }
        int oldWeight = -1;
        if (representation_ == ADJACENCY_MATRIX) {
            oldWeight = adjacencyMatrix_[from][to];
            adjacencyMatrix_[from][to] = newWeight;
        }
        else if (representation_ == ADJACENCY_LIST) {
            for (auto& e : adjacencyList_[from]) {
                if (e.first == to) {
                    oldWeight = e.second;
                    e.second = newWeight;
                    break;
                }
            }
        }
        else {
            for (auto& e : edgeList_) {
                if (e.from == from && e.to == to) {
                    oldWeight = e.weight;
                    e.weight = newWeight;
                    break;
                }
            }
        }
        return oldWeight;
    }

    void transformToAdjList() {
        if (representation_ == ADJACENCY_LIST) return;

        adjacencyList_.clear();
        adjacencyList_.resize(numVertices_);

        if (representation_ == ADJACENCY_MATRIX) {
            for (int i = 0; i < numVertices_; ++i) {
                for (int j = 0; j < numVertices_; ++j) {
                    if (adjacencyMatrix_[i][j] != 0) {
                        adjacencyList_[i].emplace_back(j, adjacencyMatrix_[i][j]);
                    }
                }
            }
        }
        else if (representation_ == EDGE_LIST) {
            for (auto const& e : edgeList_) {
                adjacencyList_[e.from].emplace_back(e.to, e.weight);
            }
        }
        representation_ = ADJACENCY_LIST;
    }

    void transformToAdjMatrix() {
        if (representation_ == ADJACENCY_MATRIX) return;

        adjacencyMatrix_.assign(numVertices_, std::vector<int>(numVertices_, 0));

        if (representation_ == ADJACENCY_LIST) {
            for (int i = 0; i < numVertices_; ++i) {
                for (auto const& e : adjacencyList_[i]) {
                    adjacencyMatrix_[i][e.first] = e.second;
                }
            }
        }
        else if (representation_ == EDGE_LIST) {
            for (auto const& e : edgeList_) {
                adjacencyMatrix_[e.from][e.to] = e.weight;
            }
        }
        representation_ = ADJACENCY_MATRIX;
    }

    void transformToListOfEdges() {
        if (representation_ == EDGE_LIST) return;

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
                for (auto const& e : adjacencyList_[i]) {
                    edgeList_.push_back({ i, e.first, e.second });
                }
            }
        }
        representation_ = EDGE_LIST;
    }

    void printGraph() const {
        if (representation_ == ADJACENCY_MATRIX) {
            std::cout << "Adjacency Matrix:\n";
            for (const auto& row : adjacencyMatrix_) {
                for (int val : row) {
                    std::cout << val << ' ';
                }
                std::cout << "\n";
            }
        }
        else if (representation_ == ADJACENCY_LIST) {
            std::cout << "Adjacency List:\n";
            for (int i = 0; i < numVertices_; ++i) {
                std::cout << i + 1 << ": ";
                for (auto const& e : adjacencyList_[i]) {
                    std::cout << "(" << e.first + 1 << ", " << e.second << ") ";
                }
                std::cout << "\n";
            }
        }
        else {
            std::cout << "Edge List:\n";
            for (auto const& e : edgeList_) {
                std::cout << e.from + 1 << " -> " << e.to + 1 << " (" << e.weight << ")\n";
            }
        }
    }

    bool isConnected() {
        std::vector<bool> visited(numVertices_, false);
        int startVertex = 0;
        for (; startVertex < numVertices_; ++startVertex) {
            if (!adjacencyList_[startVertex].empty()) break;
        }
        if (startVertex == numVertices_) return true;

        std::queue<int> q;
        q.push(startVertex);
        visited[startVertex] = true;
        int count = 1;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (auto const& e : adjacencyList_[u]) {
                if (!visited[e.first]) {
                    visited[e.first] = true;
                    q.push(e.first);
                    count++;
                }
            }
        }

        int edgesPresent = 0;
        for (int i = 0; i < numVertices_; ++i) {
            if (!adjacencyList_[i].empty()) edgesPresent++;
        }
        return count == edgesPresent;
    }

private:
    struct Edge {
        int from;
        int to;
        int weight;
    };

    Representation representation_;
    GraphType graphType_;

    int numVertices_;
    int numEdges_;

    std::vector<std::vector<int>> adjacencyMatrix_;
    std::vector<std::vector<std::pair<int, int>>> adjacencyList_;
    std::vector<Edge> edgeList_;

    void loadAdjacencyMatrix(std::ifstream& in) {
        in >> numVertices_;
        adjacencyMatrix_.assign(numVertices_, std::vector<int>(numVertices_, 0));
        numEdges_ = 0;
        for (int i = 0; i < numVertices_; ++i) {
            for (int j = 0; j < numVertices_; ++j) {
                in >> adjacencyMatrix_[i][j];
                if (adjacencyMatrix_[i][j] != 0) numEdges_++;
            }
        }
    }

    void loadAdjacencyList(std::ifstream& in) {
        in >> numVertices_;
        adjacencyList_.clear();
        adjacencyList_.resize(numVertices_);
        int edgesCount = 0;
        for (int i = 0; i < numVertices_; ++i) {
            int deg;
            in >> deg;
            while (deg--) {
                int to, weight = 1;
                in >> to;
                if (graphType_ == WEIGHTED) in >> weight;
                adjacencyList_[i].emplace_back(to - 1, weight);
                edgesCount++;
            }
        }
        numEdges_ = edgesCount;
    }

    void loadEdgeList(std::ifstream& in) {
        int vertices, edges;
        in >> vertices >> edges;
        numVertices_ = vertices;
        edgeList_.clear();
        for (int i = 0; i < edges; ++i) {
            int from, to, weight = 1;
            in >> from >> to;
            if (graphType_ == WEIGHTED) in >> weight;
            edgeList_.push_back({ from - 1, to - 1, weight });
        }
        numEdges_ = edges;
    }
};

int main() {
    freopen("output.txt", "w", stdout);

    Graph graph;
    graph.readGraph("input.txt");

    bool circleExist = false;
    int v = graph.checkEuler(circleExist);

    if (v == 0) {
        std::cout << "0\n";
        return 0;
    }
    else {
        // std::vector<int> path = graph.getEuleranTourEffective();
        std::vector<int> path = graph.getEuleranTourFleri();

        std::cout << path[0] << "\n";

        for (int i = 0; i < path.size(); ++i) {
            std::cout << path[i] << ((i + 1 < path.size()) ? " " : "\n");
        }
    }

    return 0;
}
