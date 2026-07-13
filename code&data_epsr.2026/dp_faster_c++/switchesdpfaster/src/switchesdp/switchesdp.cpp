#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <limits>
#include <chrono>
#include <iomanip>
#include <tsl/sparse_map.h>


#include "../exceptions/restexception.h"
#include "../instance/node.h"
#include "../instance/arc.h"
#include "../inputoutput/read.h"
#include "../inputoutput/readinstance.h"
#include "../instance/instance.h"
#include "../network/network.h"
#include "switchesdp.h"


#define HYBRID_SEARCH_ON



using namespace std;

SwitchesDP::SwitchesDP(Network &graph) :
		tree(graph) {
}

SwitchesDP::~SwitchesDP() {
	// TODO Auto-generated destructor stub
}

unsigned long SwitchesDP::hash(unsigned i, unsigned j, unsigned k) {
	return (unsigned long) i * 2500000000 + j * 50000 + k;
}

double SwitchesDP::computeENSTilde(Vertex &i, Vertex &j, unsigned k) {
	double sol;
	// first check memoization table
	unsigned long key = hash(tree.getNodeIndex(i), tree.getNodeIndex(j), k);
	tsl::sparse_map<unsigned long, float>::const_iterator mapIt =
			memoENStilde.find(key);
	recursionTilde++;
	if (mapIt == memoENStilde.end()) {
		// no more branches
		if (tree.getNextSibling(j) == TreeGraph::null_vertex()) {
			sol = computeENS(i, j, k);
		} else {

			#ifdef HYBRID_SEARCH_ON
			//Complexity O(log k) for the ternary search of the best projection
			int l = 0;
			int r = k;
			while (r > l) {

				//Find the mid1 and mid2
				int m1 = l + (r - l) / 3;
				int m2 = r - (r - l) / 3;

				double f1 = computeENS(i, j, m1) + computeENSTilde(i, tree.getNextSibling(j), k - m1); //evaluates the function at m1
				double f2 = computeENS(i, j, m2) +computeENSTilde(i, tree.getNextSibling(j), k - m2); //evaluates the function at m2

				if (f1 > f2)
					l = m1 + 1;
				else
					r = m2 - 1;

			}
			sol = computeENS(i, j, r) + computeENSTilde(i, tree.getNextSibling(j), k - r);

			double best = sol;
			int firstIndex = r;
			int bestIndex = r;

			// explore to the right of local minimum
			for (int right = firstIndex+1; right <= k && (computeENS(i, j, k) + computeENSTilde(i, tree.getNextSibling(j), k-right) < best); right++) {
				double now = computeENS(i, j, right) + computeENSTilde(i, tree.getNextSibling(j), k - right);
				if (best > now) {
					best = now;
					bestIndex = right;
				}
			}

			// explore to the left or local minimum
			for (int left = firstIndex-1; left >= 0 && (computeENS(i, j, left) + computeENSTilde(i, tree.getNextSibling(j), k) < best); left--) {
				double now = computeENS(i, j, left) + computeENSTilde(i, tree.getNextSibling(j), k - left);
				if (best > now) {
					best = now;
					bestIndex = left;
				}
			}

			sol = computeENS(i, j, bestIndex) + computeENSTilde(i, tree.getNextSibling(j), k - bestIndex);
			#else
			// Complexity O(k) for the linear search of the best projection
			sol = std::numeric_limits<double>::max();
			for (unsigned q = 0; q <= k; q++) {
				sol = std::min(sol,
						computeENS(i, j, q)
								+ computeENSTilde(i, tree.getNextSibling(j),
										k - q));
			}
			#endif

		}
		// insert solution value into the memoization table
		memoENStilde[key] = sol;
	} else
		sol = mapIt->second;		//return solution from table
	return sol;

}

double SwitchesDP::computeENS(Vertex &i, Vertex &j, unsigned k) {

	double sol;
	// first check memoization table
	unsigned long key = hash(tree.getNodeIndex(i), tree.getNodeIndex(j), k);
	tsl::sparse_map<unsigned long, float>::const_iterator mapIt =
			memoENS.find(key);
//	std::unordered_map<unsigned long, float>::const_iterator mapIt =
//				memoENS.find(key);
	recursion++;
	if (mapIt == memoENS.end()) {
		unsigned degree = tree.getOutDegree(j);
		// if current node is the root or it is a substation, then its
		// preceding arc already contains a switch (circuit breaker).
//		if (tree.getNode(j).getType() == Node::SUBSTATION
//				|| tree.getNode(j).getType() == Node::ROOT) {
//
//			if (degree == 1) {
//
//				sol = computeENS(j, tree.getFirstSuc(j), k)
//						+ tree.getNode(j).getFailureTime()
//								* tree.getNode(j).getActiveLoadSum();
//
//			} else {
//
//				sol = computeENSTilde(j, tree.getFirstSuc(j), k)
//						+ tree.getNode(j).getFailureTime()
//								* tree.getNode(j).getActiveLoadSum();
//
//			}
//		} else //{
		//{
		// if number of switches is greater than the number of arcs in subtree
		if (k > tree.getNode(j).getNumSuc()) {

			sol = tree.getNode(j).getMinENS();

			// if there are no more available switches
		} else if (k == 0) {

			sol = tree.getNode(j).getFailureTimeSum()
					* tree.getNode(i).getActiveLoadSum();

			// j is not branched (and k > 0)
		} else if (degree == 1) {

			sol = std::min(
					computeENS(j, tree.getFirstSuc(j), k - 1)
							+ tree.getNode(j).getFailureTime()
									* tree.getNode(j).getActiveLoadSum(),
					computeENS(i, tree.getFirstSuc(j), k)
							+ tree.getNode(j).getFailureTime()
									* tree.getNode(i).getActiveLoadSum());

			// j is branched (and k > 0)
		} else {
			sol = std::min(
					computeENSTilde(j, tree.getFirstSuc(j), k - 1)
							+ tree.getNode(j).getFailureTime()
									* tree.getNode(j).getActiveLoadSum(),
					computeENSTilde(i, tree.getFirstSuc(j), k)
							+ tree.getNode(j).getFailureTime()
									* tree.getNode(i).getActiveLoadSum());
		}
		//}
		// insert solution value into the memoization table
		memoENS[key] = sol;
	} else
		sol = mapIt->second;		//return solution from table
	return sol;
}

void SwitchesDP::recoverSwitchesTilde(Vertex &i, Vertex &j, unsigned k) {
	double sol, solmin;
	unsigned qmin;
	if (tree.getNextSibling(j) == TreeGraph::null_vertex()) {
		recoverSwitches(i, j, k);
	}
	else {
		// Complexity O(k) for the linear search of the best projection
		solmin = std::numeric_limits<double>::max();
		qmin=0;
		for (unsigned q = 0; q <= k; q++) {
			sol = memoENS[hash(i, j, q)]	+ memoENStilde[hash(i, tree.getNextSibling(j), k - q)];
			if (solmin > sol){
				solmin = sol;
				qmin = q;
			}
		}
		recoverSwitches(i, j, qmin);
		recoverSwitchesTilde(i, tree.getNextSibling(j), k-qmin);
	}
}

void SwitchesDP::recoverSwitches(Vertex &i, Vertex &j, unsigned k) {
	double sol1, sol2;
	std::pair<double, std::vector<Edge>> compSol, compSol1, compSol2;

	unsigned degree = tree.getOutDegree(j);
	// if current node is the root or it is a substation, then its
	// preceding arc already contains a switch (circuit breaker).
//		if (tree.getNode(j).getType() == Node::SUBSTATION
//				|| tree.getNode(j).getType() == Node::ROOT) {
//
//			if (degree == 1) {
//
//				recoverSwitches(j, tree.getFirstSuc(j), k);
//
//			} else {
//
//				recoverSwitchesTilde(j, tree.getFirstSuc(j), k);
//
//			}
//		} else //{
	//{
	// if number of switches is greater than the number of arcs in subtree
	if (k > tree.getNode(j).getNumSuc()) {

		for (auto itv = tree.getVerticesInTree(j).begin();
				itv != tree.getVerticesInTree(j).end(); itv++)
			switches.push_back(tree.getInEdge(*itv));

		switches.push_back(tree.getInEdge(j));

		// if there are no more available switches
	} else if (k == 0) {

		//do nothing

		// j is not branched (and k > 0)
	} else if (degree == 1) {

		sol1 = memoENS[hash(j, tree.getFirstSuc(j), k - 1)];
		sol1 += tree.getNode(j).getFailureTime()
				* tree.getNode(j).getActiveLoadSum();

		sol2 = memoENS[hash(i, tree.getFirstSuc(j), k)];
		sol2 += tree.getNode(j).getFailureTime()
				* tree.getNode(i).getActiveLoadSum();

		if (sol1 < sol2) {
			switches.push_back(tree.getInEdge(j));
			recoverSwitches(j, tree.getFirstSuc(j), k - 1);
		} else {
			recoverSwitches(i, tree.getFirstSuc(j), k);
		}
		// j is branched (and k > 0)
	} else {
		sol1 = memoENStilde[hash(j, tree.getFirstSuc(j), k - 1)];
		sol1 += tree.getNode(j).getFailureTime()
				* tree.getNode(j).getActiveLoadSum();

		sol2 = memoENStilde[hash(i, tree.getFirstSuc(j), k)];
		sol2 += tree.getNode(j).getFailureTime()
				* tree.getNode(i).getActiveLoadSum();

		if (sol1 < sol2) {
			switches.push_back(tree.getInEdge(j));
			recoverSwitchesTilde(j, tree.getFirstSuc(j), k - 1);
		} else {
			recoverSwitchesTilde(i, tree.getFirstSuc(j), k);
		}
	}

}

void SwitchesDP::run(unsigned kmax) {
	double ENDmin = computeENS(tree.getRoot(), tree.getRoot(), kmax);
	std::cout << "k = " << kmax << " \t END = " << ENDmin << std::endl;
	recoverSwitches(tree.getRoot(), tree.getRoot(), kmax);
	std::cout << "switches in arcs: " << switches.size() << std::endl;
//	std::copy(switches.begin(), switches.end(),
//					std::ostream_iterator<Edge>(std::cout, ", "));
	for(unsigned k = 0; k < switches.size(); k++)
		//std::cout << tree.getArc(switches[k]).getNumber() << endl;
		std::cout <<  tree.getArc(switches[k]).getSource() << "\t" <<  tree.getArc(switches[k]).getTarget() << endl;

}

void SwitchesDP::run(unsigned kmax, unsigned limit_time) {
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
			stop - start);
	for (unsigned k = 0; k < kmax && duration.count() < limit_time * 1000;
			k++) {
		auto partialstart = std::chrono::high_resolution_clock::now();
		double minENS = computeENS(tree.getRoot(), tree.getRoot(), k);
		stop = std::chrono::high_resolution_clock::now();
		auto partialduration = std::chrono::duration_cast<
				std::chrono::milliseconds>(stop - partialstart);
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(
				stop - start);
		
		std::cout << fixed << std::setprecision(2) << "k = " << k << " \t ENS = " << minENS; //<< endl;
		recoverSwitches(tree.getRoot(), tree.getRoot(), k);
		std::cout << " switches in: ";
		std::cout << " : ";
		std::copy(switches.begin(), switches.end(),
				std::ostream_iterator<Edge>(std::cout, ", "));
		std::cout << std::setprecision(3) << " time: " << partialduration.count() / 1000.0 << " secs" << endl;
//				<< "memoENS.size(): " << memoENS.size()  << " memoENStilde.size(): " << memoENStilde.size() << endl;
//		recoverSwitches(tree.getRoot(), tree.getRoot(), k);
//		std::cout << "switches in arcs: " << switches.size() << std::endl;
//		for(unsigned k = 0; k < switches.size(); k++)
//			std::cout <<  tree.getArc(switches[k]).getSource() << "\t" <<  tree.getArc(switches[k]).getTarget() << endl;

		clearSwitches();

	}
	std::cout << "Total time: = " << duration.count() / 1000.0 << " secs"
			<< endl;

}

