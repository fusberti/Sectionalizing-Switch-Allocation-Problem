/*
 * DpSwitches.h
 *
 *  Created on: 15 de mai. de 2024
 *      Author: jfvg
 */

#ifndef DPSWITCHES_H
#define DPSWITCHES_H

class Network;

class SwitchesDP {

	Network &tree;

	tsl::sparse_map<unsigned long, float> memoENS;
	tsl::sparse_map<unsigned long, float> memoENStilde;

	std::vector<Edge> switches;

	unsigned long hash(unsigned i, unsigned j, unsigned k) ;
	double computeENS(Vertex &i, Vertex &j, unsigned k);
	double computeENSTilde(Vertex &i, Vertex &j, unsigned k);
	void recoverSwitches(Vertex &i, Vertex &j, unsigned k);
	void recoverSwitchesTilde(Vertex &i, Vertex &j, unsigned k);
	void clearSwitches(){switches.clear();};

	unsigned recursion = 0;
	unsigned recursionTilde = 0;

public:
	SwitchesDP(Network &graph);
	virtual ~SwitchesDP();
	void run(unsigned kmax);
	void run(unsigned kmax, unsigned limit_time);
};

#endif /*DPSWITCHES_H*/
