#ifndef NETWORK_H_
#define NETWORK_H_

class Node;
class Arc;
class Instance;

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS,
		Node, Arc> TreeGraph;
typedef TreeGraph::edge_descriptor Edge;
typedef TreeGraph::vertex_descriptor Vertex;
typedef std::set<Edge> EdgeSet;
typedef TreeGraph::edge_iterator Edge_iterator;
typedef TreeGraph::out_edge_iterator Edge_out_iterator;
typedef TreeGraph::in_edge_iterator Edge_in_iterator;
typedef TreeGraph::vertex_iterator Vertex_iterator;

class Network {

	Instance &instance;
	TreeGraph graph;
	Vertex root;

//	std::vector<unsigned> customerSum;
//	std::vector<unsigned> numSucessors;
//	std::vector<unsigned> numNodesNextBranches;
//	std::vector<double> failureTimeSum;
//	std::vector<double> loadSum;
//	std::vector<double> minENS;
//	std::vector<double> minENSBranches;

	std::map<unsigned long, Vertex> mapNodes;
	std::unordered_map<Vertex, unsigned long> mapIdxNodes;
	std::map<Edge, unsigned long> mapArcs;
	std::unordered_map<Vertex, std::vector<Edge>> mapIncidents;
	std::unordered_map<Vertex, Edge> mapInEdges;
	std::unordered_map<Vertex, std::vector<Vertex>> mapOutVertices;
	std::vector<Vertex> nextSibling;

//	std::vector<Vertex> sucessors;
	std::vector<Vertex> predecessors;
	std::vector<Vertex> depth;
	std::vector<std::vector<Vertex>> verticesInTree;
	void addVertices(Vertex &current, std::vector<Vertex> &verticesInTree);
	void fillVerticesInTree();

	unsigned short calcClientsSum(Vertex &v);
	unsigned short calcNumSucessors(Vertex &v);
	double calcFailureTimeSum(Vertex &v);
	double calcDemandSum(Vertex &v);
	double calcENSMin(Vertex &v);
	void findNextBranchesNodes(Vertex &v);

	void clearAll();

public:
	Network(Instance &inst);
	~Network();

	inline Instance& getInstance() {
		return instance;
	}
	inline TreeGraph& getGraph() {
		return graph;
	}
	inline Vertex& getDepth(Vertex vert) {
		return depth[vert];
	}
	inline Vertex& getPredecessor(Vertex vert) {
		return predecessors[vert];
	}
	inline unsigned short getNumArcs() {
		return boost::num_edges(graph);
	}
	inline Vertex getSource(Edge edge) {
		return boost::source(edge, graph);
	}
	inline Vertex getTarget(Edge edge) {
		return boost::target(edge, graph);
	}
	inline Vertex& getRoot() {
		return root;
	}
	inline unsigned long getArcIndex(Edge edge) {
		return mapArcs[edge];
	}
	inline Vertex& getVertex(unsigned long idx) {
		return mapNodes[idx];
	}
	inline std::vector<Edge>& getIncidentArcs(Vertex v) {
		return mapIncidents[v];
	}
	inline double getArcLength(Edge edge) {
		return graph[edge].getLength();
	}
	inline Edge& getInEdge(Vertex v) {
		return mapInEdges[v];
	}
	inline std::vector<Vertex>& getOutVertices(Vertex v) {
		return mapOutVertices[v];
	}
	inline Vertex& getFirstSuc(Vertex v) {
		return mapOutVertices[v][0];
	}
	inline Node& getNode(Vertex v) {
		return graph[v];
	}
	inline Arc getArc(Edge edge) {
		return graph[edge];
	}
	inline Edge getEdge(Vertex source, Vertex target) {
		return boost::edge(source, target, graph).first;
	}
	inline unsigned long getNodeIndex(Vertex v) {
		return mapIdxNodes[v];
	}
	inline Vertex& getNextSibling(Vertex &v) {
		return nextSibling[getNodeIndex(v)];
	}
	inline unsigned getOutDegree(Vertex &v) {
		return boost::out_degree(v, graph);
	}            //return mapOutVertices[v].size();}
	inline unsigned short getNumNodes() {
		return boost::num_vertices(graph);
	}
	inline std::vector<Vertex>& getVerticesInTree(unsigned idx) {
		return verticesInTree[idx];
	}

	inline void setMapNodes(std::map<unsigned long, Vertex> map) {
		mapNodes = map;
	}
	inline void setIdxMapNodes(std::unordered_map<Vertex, unsigned long> map) {
		mapIdxNodes = map;
	}
	inline void setMapArcs(std::map<Edge, unsigned long> map) {
		mapArcs = map;
	}

	void initialize();
	void createStructures();
	void computeAcumulateValues();
	void print(Vertex &current);

};

#endif /* NETWORK_H_ */
