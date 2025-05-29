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
    DSU(int N) : leader(N), depth(N, 0)
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
        if (leftRoot != rightRoot)
        {
            if (depth[leftRoot] < depth[rightRoot])
                leader[leftRoot] = rightRoot;
            else
            {
                leader[rightRoot] = leftRoot;
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
            std::cerr << "Не получается открыть файл" << fileName << '\n';
            return;
        }

        char symbol;
        dataSource >> symbol;

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
            std::cerr << "Неизвестный символ представления графа";
            return;
        }

        dataSource >> vertexCount;
        int oriented, weighted;
        dataSource >> oriented >> weighted;
        graphkind_ = weighted ? WEIGHTED : NO_WEIGHTS;

        if (format_ == EDGE_MATRIX)
        {
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
                edgeCount /= 2;
            }
        }
        else if (format_ == CONNECTIVITY_LIST)
        {
            adjacentVerticesList_.assign(vertexCount, {});
            std::string line;
            std::getline(dataSource, line);
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
                edgeCount /= 2;
            }
        }
        else if (format_ == VERTEX_LINKS)
        {
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
        const auto numComponents = countComponents();
        if (numComponents > 1)
            return 0;

        int oddDegreeVertexCount = 0;
        int initialVertex = -1;
        for (const auto &edges : adjacentVerticesList_)
        {
            const auto valence = edges.size();
            if (valence % 2 == 1)
            {
                oddDegreeVertexCount++;
                if (initialVertex == -1)
                    initialVertex = &edges - &adjacentVerticesList_[0];
            }
            if (valence > 0 && initialVertex == -1)
            {
                initialVertex = &edges - &adjacentVerticesList_[0];
            }
        }
        circleExist = (oddDegreeVertexCount == 0);
        return initialVertex + 1;
    }

    int countComponents() const
    {
        std::vector<bool> marked(vertexCount, false);
        int components = 0;
        for (int i = 0; i < vertexCount; ++i)
        {
            if (!marked[i])
            {
                dfs(i, marked);
                ++components;
            }
        }
        return components;
    }

    void dfs(int vertex, std::vector<bool> &marked) const
    {
        marked[vertex] = true;
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

        Graph graphBuffer = *this;
        graphBuffer.transformToAdjList();
        std::stack<int> activeRoute;
        init--;
        activeRoute.push(init);

        while (!activeRoute.empty())
        {
            int curr_peak = activeRoute.top();
            int nextNode = -1;
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

        std::reverse(response.begin(), response.end());
        return response;
    }

    std::vector<int> getEuleranTourEffective()
    {
        std::vector<int> response;
        bool circleExist;
        int init = checkEuler(circleExist);
        if (init == 0)
            return response;

        std::vector<std::vector<int>> temporaryAdjacency(adjacentVerticesList_.size());
        for (size_t i = 0; i < adjacentVerticesList_.size(); ++i)
        {
            for (const auto &edge : adjacentVerticesList_[i])
            {
                temporaryAdjacency[i].push_back(edge.first);
            }
        }

        std::stack<int> activeRoute;
        init--;
        activeRoute.push(init);

        while (!activeRoute.empty())
        {
            int curr_peak = activeRoute.top();

            if (!temporaryAdjacency[curr_peak].empty())
            {
                int v = temporaryAdjacency[curr_peak].back();
                temporaryAdjacency[curr_peak].pop_back();
                activeRoute.push(v);
            }
            else
            {
                response.push_back(curr_peak + 1);
                activeRoute.pop();
            }
        }

        std::reverse(response.begin(), response.end());
        return response;
    }

    // void addEdge(int from, int to, int weight = 1)
    // {

    //     --from;
    //     --to;

    //     if (from < 0 || from >= vertexCount || to < 0 || to >= vertexCount)
    //     {
    //         std::cerr << "Неверный индекс вершины" << std::endl;
    //         return;
    //     }

    //     switch (format_)
    //     {
    //     case EDGE_MATRIX:
    //     {
    //         if (graphkind_ == WEIGHTED)
    //         {
    //             edgeMatrix_[from][to] = weight;
    //         }
    //         else
    //         {
    //             edgeMatrix_[from][to] = 1;
    //         }
    //         if (from != to)
    //         {
    //             if (edgeMatrix_[to][from] != 0)
    //             {
    //                 edgeCount++;
    //             }
    //         }
    //         else
    //         {
    //             edgeCount++;
    //         }
    //         break;
    //     }
    //     case CONNECTIVITY_LIST:
    //     {
    //         auto &edges = adjacentVerticesList_[from];
    //         if (std::none_of(edges.begin(), edges.end(), [to](const auto &edge)
    //                          { return edge.first == to; }))
    //         {
    //             edges.emplace_back(to, weight);
    //             edgeCount++;
    //         }
    //         break;
    //     }
    //     case VERTEX_LINKS:
    //     {
    //         edgeList_.push_back({from, to, weight});
    //         edgeCount++;
    //         break;
    //     }
    //     }
    // }

    void removeEdge(int from, int to)
    {
        if (from < 0 || from >= vertexCount || to < 0 || to >= vertexCount)
        {
            std::cerr << "Неверный индекс вершины" << std::endl;
            return;
        }

        if (format_ == EDGE_MATRIX)
        {
            if (edgeMatrix_[from][to] != 0)
            {
                edgeMatrix_[from][to] = 0;
                edgeCount--;
            }
        }

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

    // int changeEdge(int from, int to, int newWeight)
    // {
    //     --from;
    //     --to;

    //     if (from < 0 || from >= vertexCount || to < 0 || to >= vertexCount)
    //     {
    //         std::cerr << "Неверный индекс вершины" << std::endl;
    //         return -1;
    //     }

    //     auto &adj = adjacentVerticesList_[from];
    //     auto it = std::lower_bound(adj.begin(), adj.end(), to, [](const auto &l, int r)
    //                                { return l.first < r; });

    //     if (it != adj.end() && it->first == to)
    //     {
    //         int oldWeight = it->second;
    //         it->second = newWeight;
    //         return oldWeight;
    //     }

    //     return -1;
    // }

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
            for (int i = 0; i < vertexCount; ++i)
            {
                for (int j = 0; j < vertexCount; ++j)
                {
                    if (edgeMatrix_[i][j] != 0)
                    {
                        adjacentVerticesList_[i].push_back({j, edgeMatrix_[i][j]});
                    }
                }
            }
        }
        else if (format_ == VERTEX_LINKS)
        {
            for (const auto &edge : edgeList_)
            {
                adjacentVerticesList_[edge.from].push_back({edge.to, edge.weight});
            }
        }
        format_ = CONNECTIVITY_LIST;
    }

    // void transformToAdjMatrix()
    // {
    //     if (format_ == EDGE_MATRIX)
    //         return;

    //     edgeMatrix_.resize(vertexCount, std::vector<int>(vertexCount, 0));

    //     if (format_ == CONNECTIVITY_LIST)
    //     {
    //         for (int vertex = 0; vertex < vertexCount; ++vertex)
    //         {
    //             auto &edges = adjacentVerticesList_[vertex];
    //             for (auto it = edges.begin(); it != edges.end(); ++it)
    //             {
    //                 int target = it->first;
    //                 int weight = it->second;
    //                 edgeMatrix_[vertex][target] = weight;
    //             }
    //         }
    //     }
    //     else if (format_ == VERTEX_LINKS)
    //     {
    //         for (const auto &edge : edgeList_)
    //         {
    //             edgeMatrix_[edge.from][edge.to] = edge.weight;
    //         }
    //     }
    //     format_ = EDGE_MATRIX;
    // }

    // void transformToListOfEdges()
    // {
    //     if (format_ == VERTEX_LINKS)
    //         return;

    //     edgeList_.clear();
    //     edgeList_.reserve(edgeCount);

    //     if (format_ == EDGE_MATRIX)
    //     {
    //         for (int from = 0; from < vertexCount; ++from)
    //         {
    //             for (int to = 0; to < vertexCount; ++to)
    //             {
    //                 int weight = edgeMatrix_[from][to];
    //                 if (weight != 0)
    //                 {
    //                     edgeList_.emplace_back(Edge{from, to, weight});
    //                 }
    //             }
    //         }
    //     }
    //     else if (format_ == CONNECTIVITY_LIST)
    //     {
    //         for (int from = 0; from < vertexCount; ++from)
    //         {
    //             for (const auto &edge : adjacentVerticesList_[from])
    //             {
    //                 edgeList_.emplace_back(Edge{from, edge.first, edge.second});
    //             }
    //         }
    //     }
    //     format_ = VERTEX_LINKS;
    // }

    // void writeGraph(const std::string &filename) const
    // {
    //     std::ofstream file{filename, std::ios::binary};
    //     if (!file.is_open())
    //     {
    //         throw std::runtime_error{"Ошибка открытия файла для следующей записи" + filename};
    //     }

    //     char representationCode = format_ == EDGE_MATRIX ? 'C' : format_ == CONNECTIVITY_LIST ? 'L'
    //                                                                                                             : 'E';
    //     file.write(&representationCode, sizeof(representationCode));
    //     auto writeInt = [&file](int value)
    //     {
    //         file.write(reinterpret_cast<const char *>(&value), sizeof(value));
    //     };
    //     writeInt(vertexCount);
    //     if (format_ == VERTEX_LINKS)
    //     {
    //         writeInt(edgeCount);
    //     }
    //     char graphTypeCode = graphkind_ == WEIGHTED ? '1' : '0';
    //     file.write(&graphTypeCode, sizeof(graphTypeCode));

    //     if (format_ == EDGE_MATRIX)
    //     {
    //         file.write(reinterpret_cast<const char *>(edgeMatrix_.data()),
    //                    vertexCount * vertexCount * sizeof(int));
    //     }
    //     else if (format_ == CONNECTIVITY_LIST)
    //     {
    //         for (const auto &edges : adjacentVerticesList_)
    //         {
    //             writeInt(edges.size());
    //             file.write(reinterpret_cast<const char *>(edges.data()), edges.size() * sizeof(std::pair<int, int>));
    //         }
    //     }
    //     else if (format_ == VERTEX_LINKS)
    //     {
    //         file.write(reinterpret_cast<const char *>(edgeList_.data()), edgeCount * sizeof(Edge));
    //     }
    // }

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
        if (vertexCount == 0)
        {
            return true;
        }

        int initialVertex = -1;
        for (int vertex = 0; vertex < vertexCount; ++vertex)
        {
            if (!adjacentVerticesList_[vertex].empty())
            {
                initialVertex = vertex;
                break;
            }
        }

        if (initialVertex == -1)
        {
            return true;
        }

        std::vector<bool> marked(vertexCount, false);
        std::stack<int> verticesToVisit;
        verticesToVisit.push(initialVertex);
        marked[initialVertex] = true;

        int markedkol = 0;
        while (!verticesToVisit.empty())
        {
            int currentVertex = verticesToVisit.top();
            verticesToVisit.pop();
            ++markedkol;

            for (const auto &edge : adjacentVerticesList_[currentVertex])
            {
                int nhbr = edge.first;
                if (edge.second > 0 && !marked[nhbr])
                {
                    marked[nhbr] = true;
                    verticesToVisit.push(nhbr);
                }
            }
        }

        return markedkol == std::count_if(adjacentVerticesList_.begin(), adjacentVerticesList_.end(),
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