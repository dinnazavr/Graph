#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <queue>
#include <limits>

class DSU
{
private:
    std::vector<int> leader;
    std::vector<int> depth;

public:
    DSU(int N)
        : leader(N), depth(N, 0)
    {
        std::iota(leader.begin(), leader.end(), 0);
    }

    int find(int x)
    {
        if (x == leader[x])
            return x;
        return leader[x] = find(leader[x]);
    }

    void unite(int x, int y)
    {
        int leftRoot = find(x);
        int rightRoot = find(y);

        // If both elements have different roots, unite them
        if (leftRoot != rightRoot)
        {
            // Make the root with lesser depth point to the root with greater depth
            if (depth[leftRoot] < depth[rightRoot])
                leader[leftRoot] = rightRoot;
            else
            {
                leader[rightRoot] = leftRoot;
                // If both roots have the same depth, increment the depth of the new root
                if (depth[leftRoot] == depth[rightRoot])
                    ++depth[leftRoot];
            }
        }
    }
};

class Graph
{
public:
    enum Format
    {
        EDGE_MATRIX,
        CONNECTIVITY_LIST,
        VERTEX_LINKS
    };

    enum GraphKind
    {
        NO_WEIGHTS,
        WEIGHTED
    };

    Graph() : format_(EDGE_MATRIX), graphkind_(NO_WEIGHTS), vertexCount(0), edgeCount(0) {}

    void readGraph(const std::string &fileName)
    {
        std::ifstream dataSource{fileName};
        if (!dataSource.is_open())
        {
            std::cerr << "This file does not open" << fileName << '\n';
            return;
        }

        char symbol;
        dataSource >> symbol;

        // Determine the graph format based on the symbol
        switch (symbol)
        {
        case 'C':
            format_ = EDGE_MATRIX;
            break;
        case 'L':
            format_ = CONNECTIVITY_LIST;
            break;
        case 'E':
            format_ = VERTEX_LINKS;
            break;
        default:
            std::cerr << "There is no such symbol to represent a graph!";
            return;
        }

        dataSource >> vertexCount;
        int oriented, weighted;
        dataSource >> oriented >> weighted;
        graphkind_ = weighted ? WEIGHTED : NO_WEIGHTS;

        if (format_ == EDGE_MATRIX)
        {
            // Load the adjacency matrix
            edgeMatrix_.assign(vertexCount, std::vector<int>(vertexCount));
            for (int i = 0; i < vertexCount; ++i)
            {
                for (int j = 0; j < vertexCount; ++j)
                {
                    dataSource >> edgeMatrix_[i][j];
                    edgeCount += (edgeMatrix_[i][j] != 0);
                }
            }
            if (!oriented)
            {
                edgeCount /= 2; // Account for undirected edges
            }
        }
        else if (format_ == CONNECTIVITY_LIST)
        {
            // Load the adjacency list
            adjacentVerticesList_.assign(vertexCount, {});
            std::string line;
            std::getline(dataSource, line); // Consume the remainder of the current line
            for (int i = 0; i < vertexCount; ++i)
            {
                std::getline(dataSource, line);
                std::istringstream ss{line};
                int nhbr, weight = 1;
                while (ss >> nhbr)
                {
                    if (graphkind_ == WEIGHTED)
                    {
                        ss >> weight;
                    }
                    adjacentVerticesList_[i].emplace_back(nhbr - 1, weight);
                    ++edgeCount;
                }
            }
            if (!oriented)
            {
                edgeCount /= 2; // Account for undirected edges
            }
        }
        else if (format_ == VERTEX_LINKS)
        {
            // Load the edge list
            edgeList_.reserve(edgeCount);
            for (int i = 0; i < edgeCount; ++i)
            {
                int curr_peak, v, weight = 1;
                dataSource >> curr_peak >> v;
                if (graphkind_ == WEIGHTED)
                {
                    dataSource >> weight;
                }
                edgeList_.emplace_back(Edge{curr_peak - 1, v - 1, weight});
            }
        }
    }

    int checkEuler(bool &circleExist)
    {
        // Count the number of connected components in the graph
        const auto numComponents = countComponents();

        // If there are more than one connected components, the graph is not Eulerian
        if (numComponents > 1)
            return 0;

        // Count the number of vertices with odd degree
        int oddDegreeVertexCount = 0;
        int initialVertex = -1;
        for (size_t i = 0; i < adjacentVerticesList_.size(); ++i)
        {
            const auto valence = adjacentVerticesList_[i].size();
            if (valence % 2 == 1)
            {
                ++oddDegreeVertexCount;
                if (initialVertex == -1)
                    initialVertex = static_cast<int>(i);
            }
            if (valence > 0 && initialVertex == -1)
                initialVertex = static_cast<int>(i);
        }

        // A graph is Eulerian if it has no odd degree vertices
        circleExist = (oddDegreeVertexCount == 0);
        return initialVertex + 1;
    }

    int countComponents() const
    {
        // Mark all vertices as unvisited
        std::vector<bool> marked(vertexCount, false);

        // Count the number of connected components
        int components = 0;
        for (int i = 0; i < vertexCount; ++i)
        {
            if (!marked[i])
            {
                // Perform a depth-first search from the current vertex
                dfs(i, marked);
                ++components;
            }
        }

        return components;
    }

    void dfs(int vertex, std::vector<bool> &marked) const
    {
        // Mark the current vertex as visited
        marked[vertex] = true;

        // Recursively visit all the neighbors of the current vertex
        for (const auto &edge : adjacentVerticesList_[vertex])
        {
            const auto nhbr = edge.first;
            if (!marked[nhbr])
            {
                dfs(nhbr, marked);
            }
        }
    }

    std::vector<int> getEuleranTourFleri()
    {
        std::vector<int> response;
        bool circleExist;
        int init = checkEuler(circleExist);
        if (init == 0)
            return response;

        // Create a temporary copy of the graph
        Graph graphBuffer = *this;
        graphBuffer.transformToAdjList();

        // Use a stack to keep track of the active route
        std::stack<int> activeRoute;
        init--;
        activeRoute.push(init);

        // Loop until the stack is empty
        while (!activeRoute.empty())
        {
            int curr_peak = activeRoute.top();
            int nextNode = -1;

            // Find the next node in the Eulerian tour
            auto it = graphBuffer.adjacentVerticesList_[curr_peak].begin();
            while (it != graphBuffer.adjacentVerticesList_[curr_peak].end())
            {
                if (it->second > 0)
                {
                    nextNode = it->first;
                    it->second--;
                    auto &revEdges = graphBuffer.adjacentVerticesList_[nextNode];
                    auto revIt = std::find_if(revEdges.begin(), revEdges.end(), [curr_peak](const auto &edge)
                                              { return edge.first == curr_peak; });
                    revIt->second--;
                    break;
                }
                ++it;
            }

            // If there is a next node, add it to the active route
            if (nextNode == -1)
            {
                response.push_back(curr_peak + 1);
                activeRoute.pop();
            }
            else
            {
                activeRoute.push(nextNode);
            }
        }

        // Reverse the response vector to get the correct order
        std::reverse(response.begin(), response.end());
        return response;
    }

    std::vector<int> getEuleranTourEffective()
    {
        std::vector<int> response;
        bool circleExist;
        int init = checkEuler(circleExist);

        // If there is no Eulerian path or circuit, return an empty response
        if (init == 0)
            return response;

        std::vector<bool> attendanced(vertexCount, false); // Track visited vertices
        std::stack<int> activeRoute;
        init--;
        activeRoute.push(init);

        while (!activeRoute.empty())
        {
            int curr_peak = activeRoute.top();

            // Mark the vertex as visited if not already done
            if (!attendanced[curr_peak])
            {
                attendanced[curr_peak] = true;
                response.push_back(curr_peak + 1);
            }

            // If no adjacent vertices remain, backtrack
            if (adjacentVerticesList_[curr_peak].empty())
            {
                activeRoute.pop();
            }
            else
            {
                // Otherwise, continue to the next vertex
                int v = adjacentVerticesList_[curr_peak].back().first;
                adjacentVerticesList_[curr_peak].pop_back();
                activeRoute.push(v);
            }
        }

        return response;
    }

    void removeEdge(int from, int to)
    {
        if (from < 0 || from >= vertexCount || to < 0 || to >= vertexCount)
        {
            std::cerr << "Uncorrect vertex index" << std::endl;
            return;
        }

        // Remove the edge from the adjacency matrix representation
        if (format_ == EDGE_MATRIX)
        {
            if (edgeMatrix_[from][to] != 0)
            {
                edgeMatrix_[from][to] = 0;
                edgeCount--;
            }
        }

        // Remove the edge from the adjacency list representation
        else if (format_ == CONNECTIVITY_LIST)
        {
            auto &adj = adjacentVerticesList_[from];
            auto it = std::lower_bound(adj.begin(), adj.end(), to, [](const auto &l, int r)
                                       { return l.first < r; });
            if (it != adj.end() && it->first == to)
            {
                adj.erase(it);
                edgeCount--;
            }
        }

        // Remove the edge from the edge list representation
        else if (format_ == VERTEX_LINKS)
        {
            auto it = std::remove_if(edgeList_.begin(), edgeList_.end(), [from, to](const Edge &edge)
                                     { return edge.from == from && edge.to == to; });
            if (it != edgeList_.end())
            {
                edgeList_.erase(it, edgeList_.end());
                edgeCount--;
            }
        }
    }

    // Transforms the graph representation to an adjacency list.
    void transformToAdjList()
    {
        if (format_ == CONNECTIVITY_LIST)
        {
            return;
        }

        adjacentVerticesList_.clear();
        adjacentVerticesList_.resize(vertexCount);
        graphkind_ = (graphkind_ == WEIGHTED) ? WEIGHTED : NO_WEIGHTS;

        if (format_ == EDGE_MATRIX)
        {
            // Iterate over the adjacency matrix
            for (int i = 0; i < vertexCount; ++i)
            {
                for (int j = 0; j < vertexCount; ++j)
                {
                    if (edgeMatrix_[i][j] != 0)
                    {
                        // Add the edge to the adjacency list
                        adjacentVerticesList_[i].push_back({j, edgeMatrix_[i][j]});
                    }
                }
            }
        }
        else if (format_ == VERTEX_LINKS)
        {
            // Iterate over the edge list
            for (const auto &edge : edgeList_)
            {
                // Add the edge to the adjacency list
                adjacentVerticesList_[edge.from].push_back({edge.to, edge.weight});
            }
        }
        format_ = CONNECTIVITY_LIST;
    }

private:
    Format format_;
    GraphKind graphkind_;
    int vertexCount;
    int edgeCount;

    std::vector<std::vector<int>> edgeMatrix_;
    std::vector<std::vector<std::pair<int, int>>> adjacentVerticesList_;

    struct Edge
    {
        int from;
        int to;
        int weight;
    };
    std::vector<Edge> edgeList_;

    bool isConnected() const
    {
        // If there are no vertices, the graph is considered connected
        if (vertexCount == 0)
        {
            return true;
        }

        // Find an initial vertex with an edge
        int initialVertex = -1;
        for (int vertex = 0; vertex < vertexCount; ++vertex)
        {
            if (!adjacentVerticesList_[vertex].empty())
            {
                initialVertex = vertex;
                break;
            }
        }

        // If no initial vertex is found, the graph is considered connected
        if (initialVertex == -1)
        {
            return true;
        }

        // Vector to keep track of visited vertices
        std::vector<bool> marked(vertexCount, false);
        // Stack for implementing DFS
        std::stack<int> verticesToVisit;
        verticesToVisit.push(initialVertex);
        marked[initialVertex] = true;

        int markedCount = 0;
        while (!verticesToVisit.empty())
        {
            int currentVertex = verticesToVisit.top();
            verticesToVisit.pop();
            ++markedCount;

            // Visit all adjacent unmarked vertices
            for (const auto &edge : adjacentVerticesList_[currentVertex])
            {
                int neighbor = edge.first;
                if (edge.second > 0 && !marked[neighbor])
                {
                    marked[neighbor] = true;
                    verticesToVisit.push(neighbor);
                }
            }
        }

        // Check if all vertices with edges are marked
        return markedCount == std::count_if(adjacentVerticesList_.begin(), adjacentVerticesList_.end(),
                                            [](const auto &edges)
                                            { return !edges.empty(); });
    }
};
int main()
{
    freopen("output.txt", "w", stdout);
    Graph graph;
    graph.readGraph("input.txt");
    bool circleExist = false;
    int v = graph.checkEuler(circleExist);
    if (v == 0)
    {
        std::cout << "0\n";
        exit(0);
    }
    else
    {
        //  std::vector<int> path = graph.getEuleranTourEffective();
        std::vector<int> path = graph.getEuleranTourFleri();
        std::cout << path[0] << "\n";
        for (int i = 0; i < path.size(); ++i)
        {
            std::cout << path[i] << ((i + 1 < path.size()) ? " " : "\n");
        }
    }
    return 0;
}