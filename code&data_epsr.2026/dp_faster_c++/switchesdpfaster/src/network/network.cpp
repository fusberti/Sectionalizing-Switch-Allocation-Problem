/*
 * Network.cpp
 *
 *  Created on: 24 de jun de 2020
 *      Author: jfvg
 */
#include <iostream>
//#include <vector>
//#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "../instance/node.h"
#include "../instance/arc.h"
#include "../instance/instance.h"
#include "network.h"

using namespace std;
//using namespace boost;

Network::Network(Instance &inst) :
		instance(inst) {
	initialize(); //initialize all the structures needed
	computeAcumulateValues(); //compute instance
}

void Network::initialize() {
	//construct all graph structures: sucessors, predecessors, depht vectors
	Instance inst = getInstance();
	NodeList nodes = inst.getNodes();
	ArcList arcs = inst.getArcs();

	for (Node node : nodes) { //inserting nodes
		Vertex v = boost::add_vertex(node, graph);
		mapNodes.insert(make_pair(node.getNumber(), v));
	}
	for (Arc arc : arcs) { //inserting arcs
		mapArcs.insert(
				make_pair(
						boost::add_edge(mapNodes[arc.getSource()],
								mapNodes[arc.getTarget()], arc, graph).first,
						arc.getNumber()));
	}
	unsigned idx = 0;
	for (Vertex vt : boost::make_iterator_range(boost::vertices(graph))) {
		mapIdxNodes.insert(make_pair(vt, idx++));
	}

	root = *vertices(graph).first;		//the root node
	createStructures();
}

void Network::createStructures() {
	depth.clear();
	predecessors.clear();
	nextSibling.clear();
	mapInEdges.clear();
	mapIncidents.clear();
	mapOutVertices.clear();
	verticesInTree.clear();
	//cout << "numvertices: " << boost::num_vertices(graph) << endl;
	depth.assign(boost::num_vertices(graph), TreeGraph::null_vertex());
	predecessors.assign(boost::num_vertices(graph), TreeGraph::null_vertex());

	boost::breadth_first_search(graph, root,
			boost::visitor(
					boost::make_bfs_visitor(
							boost::record_distances(depth.data(),
									boost::on_tree_edge { }))));

	for (auto d = depth.begin(); d != depth.end(); ++d)
		if (*d == Instance::NULLS)
			*d = 0;
		else
			(*d)++;

	breadth_first_search(graph, root,
			visitor(
					boost::make_bfs_visitor(
							boost::record_predecessors(predecessors.data(),
									boost::on_tree_edge { }))));

	nextSibling.resize(graph.m_vertices.size(), TreeGraph::null_vertex());//initializing the sibling structure

	Edge_out_iterator eo, eo_end;
	Edge_in_iterator ein, ein_end;
	for (Vertex vd : boost::make_iterator_range(boost::vertices(graph))) {
		//cout << "vd: " << vd << endl;
		std::vector<Edge> incidents;
		std::vector<Vertex> outVertices;
		Vertex last = TreeGraph::null_vertex();
		for (boost::tie(eo, eo_end) = boost::out_edges(vd, graph); eo != eo_end;
				++eo) {
			//cout << "\tedge: " << *eo << " target: " << getTarget(*eo) << " last: " << last << endl;

			incidents.push_back(*eo);
			outVertices.push_back(getTarget(*eo));
			if (last != TreeGraph::null_vertex())
				nextSibling[this->getNodeIndex(last)] = getTarget(*eo);
			last = getTarget(*eo);
		}
		for (boost::tie(ein, ein_end) = boost::in_edges(vd, graph);
				ein != ein_end; ++ein) {
			incidents.push_back(*ein);
			mapInEdges.insert(std::make_pair(vd, *ein));
		}
		mapIncidents.insert(std::make_pair(vd, incidents));
		mapOutVertices.insert(std::make_pair(vd, outVertices));

	}
	fillVerticesInTree();

}

void Network::addVertices(Vertex &current,
		std::vector<Vertex> &verticesInTree) {
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();
	while (itVertices != getOutVertices(current).end()) {
		Vertex next = *itVertices;
		verticesInTree.push_back(next);
		addVertices(next, verticesInTree);
		itVertices++;
	}
}

void Network::fillVerticesInTree() {
	verticesInTree.resize(getNumNodes());
	for (Vertex vd : boost::make_iterator_range(boost::vertices(graph))) {
		std::vector<Vertex> vertices;
		addVertices(vd, vertices);
		verticesInTree[getNodeIndex(vd)] = vertices;
	}
}

void Network::clearAll() {
	for (Vertex vd : boost::make_iterator_range(boost::vertices(graph))) {
		getNode(vd).setClientsSum(0);
		getNode(vd).setNumSuc(0);
		getNode(vd).setNumNodesNextBranches(0);
		getNode(vd).setFailureTimeSum(0);
		getNode(vd).setActiveLoadSum(0);
		getNode(vd).setMinENS(0);
		getNode(vd).setMinENSBranches(0);
	}
}

void Network::computeAcumulateValues() {
	//clearAll();
	calcNumSucessors(getRoot());
	calcClientsSum(getRoot());
	calcDemandSum(getRoot());
	calcFailureTimeSum(getRoot());
	calcENSMin(getRoot());
	findNextBranchesNodes(getRoot());
	//print(getRoot());
}

unsigned short Network::calcClientsSum(Vertex &current) {

	getNode(current).setClientsSum(getNode(current).getNumClients());
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();

	while (itVertices != getOutVertices(current).end()) {
		Vertex next = *itVertices;
		getNode(current).setClientsSum(
				getNode(current).getClientsSum() + calcClientsSum(next));
		itVertices++;
	}
	return getNode(current).getClientsSum();
}

double Network::calcDemandSum(Vertex &current) {
	getNode(current).setActiveLoadSum(getNode(current).getActiveLoad());
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();

	while (itVertices != getOutVertices(current).end()) {
		Vertex next = *itVertices;
		getNode(current).setActiveLoadSum(
				getNode(current).getActiveLoadSum() + calcDemandSum(next));
		itVertices++;
	}
	return getNode(current).getActiveLoadSum();
}

double Network::calcFailureTimeSum(Vertex &current) {
	getNode(current).setFailureTimeSum(getNode(current).getFailureTime());
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();

	while (itVertices != getOutVertices(current).end()) {
		Vertex next = *itVertices;
		getNode(current).setFailureTimeSum(
				getNode(current).getFailureTimeSum()
						+ calcFailureTimeSum(next));
		itVertices++;
	}
	return getNode(current).getFailureTimeSum();
}

unsigned short Network::calcNumSucessors(Vertex &current) {
	getNode(current).setNumSuc(0);
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();
	while (itVertices != getOutVertices(current).end()) {
		Vertex next = *itVertices;
		getNode(current).setNumSuc(
				1 + getNode(current).getNumSuc() + calcNumSucessors(next));
		itVertices++;
	}
	return getNode(current).getNumSuc();
}

double Network::calcENSMin(Vertex &current) {
	getNode(current).setMinENS(
			getNode(current).getFailureTime()
					* getNode(current).getActiveLoadSum());
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();
	while (itVertices != getOutVertices(current).end()) {
		Vertex next = *itVertices;
		getNode(current).setMinENS(
				getNode(current).getMinENS() + calcENSMin(next));
		itVertices++;
	}
	return getNode(current).getMinENS();
}

void Network::print(Vertex &current) {
	std::vector<Vertex>::const_iterator itVertices =
			getOutVertices(current).begin();
	cout << "node: " << current << "(" << getNode(current).getNumber()
			<< ") num nodes branches: "
			<< getNode(current).getNumNodesNextBranches() << " ENS Branches: "
			<< getNode(current).getMinENSBranches() << endl;
	cout << "node: " << current << "(" << getNode(current).getNumber()
			<< ") first: ";
	if (getOutDegree(current) != 0)
		cout << getFirstSuc(current);
	else
		cout << "null";
	cout << " next sibling: " << getNextSibling(current) << endl;
	while (itVertices != getOutVertices(current).end()
			&& *itVertices != TreeGraph::null_vertex()) {
		Vertex next = *itVertices;
		print(next);
		itVertices++;
	}
}

void Network::findNextBranchesNodes(Vertex &current) {
//	cout << "current: " << current << " (" << getNode(current).getNumber()
//			<< ")" << endl;

	getNode(current).setNumNodesNextBranches(1);
	getNode(current).setMinENSBranches(
			getNode(current).getFailureTime()
					* getNode(current).getActiveLoadSum());
	if (getOutDegree(current) != 0) {
		Vertex &next = getFirstSuc(current);
//		cout << "first: " << next << " (" << getNode(next).getNumber() << ")"
//				<< endl;

		findNextBranchesNodes(next);
		getNode(current).setNumNodesNextBranches(
				getNode(current).getNumNodesNextBranches()
						+ getNode(next).getNumNodesNextBranches());
		getNode(current).setMinENSBranches(
				getNode(current).getMinENSBranches()
						+ getNode(next).getMinENSBranches());
	}

	if (getNextSibling(current) != TreeGraph::null_vertex()) {
		Vertex &next = getNextSibling(current);
//		cout << "next sibling: " << next << " (" << getNode(next).getNumber() << ")"
//				<< endl;
		findNextBranchesNodes(next);
		getNode(current).setNumNodesNextBranches(
				getNode(current).getNumNodesNextBranches()
						+ getNode(next).getNumNodesNextBranches());
		getNode(current).setMinENSBranches(
				getNode(current).getMinENSBranches()
						+ getNode(next).getMinENSBranches());
	}

}

Network::~Network() {
// TODO Auto-generated destructor stub
}

