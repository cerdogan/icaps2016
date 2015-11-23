/** 
 * @author Can Erdogan
 * @file sampling.cpp
 * @date 2015-11-22
 */

#include <nlopt.h>
#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>

using namespace Eigen;
using namespace std;

#define sq(x) ((x) * (x))

/* ******************************************************************************************** */
int H_TYPE = 1;
static const int NUM_GEARS = 3;
vector <double> radii {1.4, 0.5, 1.2, 1.5, 1.2, 0.3, 1.2, 1.0, 0.8, 0.9, 1.4}; //, 1.6};
static const int NUM_STICKS = 3;
vector <double> lengths {8.0, 5.0}; 
Vector3d obj0 (0.5, 1.0, 1);
Vector3d objn (18.0, 2.0, 2);
vector <Vector3d> obs { Vector3d(5.0, 1.0, 1.0) };
Vector2d minP (1.0, -1.0), maxP (objn(0), 3.0);

/* ******************************************************************************************** */
struct Action_ {
	int use;
	bool check;
	Action_ () { check = false; use = -1; }
	string s () const {
		char buf [256];
		sprintf(buf, "use: %d, check: %d", use, (check ? 1 : 0));
		return buf;
	}
};

/* ******************************************************************************************** */
struct State {
	bool check;
	vector <bool> used, used_sticks;
	vector <int> order;
	double H;
	int d;
	string s () {
		char buf [1024];
		sprintf(buf, "\n  used: { ");
		for(int i = 0; i < used.size(); i++) sprintf(buf, "%s%d ", buf, (used[i] ? 1 : 0));
		sprintf(buf, "%s\b }\n  order: { ", buf);
		for(int i = 0; i < order.size(); i++) sprintf(buf, "%s%d ", buf, order[i]);
		sprintf(buf, "%s\b }\n  check: %d\n", buf, (check ? 1 : 0));
		return buf;
	}
};

/* ******************************************************************************************** */
struct mycomparison {
	bool reverse;
	mycomparison(const bool& revparam=1) {reverse=revparam; }
	bool operator() (State* lhs, State* rhs) const {
		if (reverse) return (lhs->H>rhs->H);
		else return (lhs->H<rhs->H);
	}
};

typedef std::priority_queue<State*,std::vector<State*>,mycomparison> pq;

/* ******************************************************************************************** */
State initialState () {
	State s;
	s.check = 0;
	s.d = 0;
	for(int i = 0; i < NUM_GEARS; i++) s.used.push_back(false);
	for(int i = 0; i < NUM_STICKS; i++) s.used_sticks.push_back(false);
	s.H = 0;
	return s;
}

/* ******************************************************************************************** */
void generateActions (const State& s, vector <Action_>& acts) {
	acts.clear();
	for(int i = 0; i < NUM_GEARS; i++) {
		if(s.used[i]) continue;
		Action_ a;
		a.use = i;
		acts.push_back(a);
		a.check = true;
		acts.push_back(a);
		int use = a.use;
		for(int j = 0; j < NUM_STICKS; j++) {
			if(s.used_sticks[j]) continue;
			a.use = -10 * j - use;
			a.check = false;
			acts.push_back(a);
			a.check = true;
			acts.push_back(a);
		}
	}
}

/* ********************************************************************************************* */
// Given two points p0, p2 and distances d01 and d12, find p1 such 
//that |p0p1|^2 = d01 and |p1p2|^2 = d12. Distances in L2 norm.
Vector2d solveL2Distance (const Vector2d& p0, const Vector2d& p2, double rA, double rB) {

	// Setup
	double xA = p0(0), yA = p0(1); 
	double xB = p2(0), yB = p2(1); 
	double d = sqrt((xB-xA)*(xB-xA) + (yB-yA)*(yB-yA));
	double K = 0.25 * sqrt(((rA+rB)*(rA+rB)-d*d)*(d*d-(rA-rB)*(rA-rB)));

	// Make a random choice between two solutions
	double rand_ (((double) rand()) / RAND_MAX);
	Vector2d sol;
	if(rand_ > 0.5) {
		double x1 = (0.5)*(xB+xA) + (0.5)*(xB-xA)*(rA*rA-rB*rB)/(d*d) - 2*(yB-yA)*K/(d*d);
		double y1 = (0.5)*(yB+yA) + (0.5)*(yB-yA)*(rA*rA-rB*rB)/(d*d) - (-2*(xB-xA)*K/(d*d));
		sol = Vector2d(x1,y1);
	}
	else {
		double x2 = (0.5)*(xB+xA) + (0.5)*(xB-xA)*(rA*rA-rB*rB)/(d*d) + 2*(yB-yA)*K/(d*d);
		double y2 = (0.5)*(yB+yA) + (0.5)*(yB-yA)*(rA*rA-rB*rB)/(d*d) + (-2*(xB-xA)*K/(d*d));
		sol = Vector2d(x2,y2);
	}

	// Check the math
	assert((abs((p0-sol).norm() - rA) < 1e-3) && (abs((p2-sol).norm() - rB) < 1e-3));
	return sol;
}

/* ******************************************************************************************** */
bool feasible (State& s, VectorXd& X, VectorXd& Y, int* totalIters = NULL) {

	if(s.order.empty()) return true;

	X = VectorXd (s.order.size());
	Y = VectorXd (s.order.size());
	for(int iter = 0; iter < 5*1e6; iter++) {

		// Check if the contact constraints can be achieved by (1) checking for x-axis distance between
		// consecutive gears, (2) computing the y-axis values randomly.
		bool fail = false;
		double lastX = obj0(0), lastY = obj0(1), lastR = obj0(2);
		int numGearsToPropagate = s.order.size() - (s.check ? 1 : 0);
		for(int i = 0; i < numGearsToPropagate; i++) {

			// Sample random x location
			double x = (((double) rand()) / RAND_MAX) * (maxP(0) - minP(0)) + minP(0);

			// Connect with a spur gear 
			double Rs, r;
			if(s.order[i] >= 0) {
				r = radii[s.order[i]];
				Rs = (lastR + r);
			}
			else {

				// Get info
				int stick_id = (-s.order[i]) / 10;
				int gear_id = (-s.order[i]) % 10;
				double length = lengths[stick_id];
				r = radii[gear_id];																				//< be careful with r
				Rs = sqrt(sq(length) + sq(lastR - r));
			}
		
			// Check distance
			if(abs(lastX - x) > Rs) {
				fail = true;
				break;
			}

			// Compute y location (y = sqrt((x1-x2)^2) 
			int random = 1 - 2 * round(((double) rand()) / RAND_MAX);
			assert((random == -1) || (random == 1));
			double y = random * sqrt(Rs*Rs - sq(lastX-x)) + lastY;

			// Check gear for collisions
			for(int o_i = 0; o_i < obs.size(); o_i++) {
				double distSQ = sq(x-obs[o_i](0)) + sq(y-obs[o_i](1));
				if(distSQ < sq(r+obs[o_i](2))) {
					fail = true;
					break;
				}
				if(fail) break;
			}

			// Check stick for collisions
			bool stickDir = true;
			if(s.order.back() < 0) {

				// Choose stick direction randomly
				stickDir = (round(((double) rand()) / RAND_MAX) == 0);
				int stick_id = (-s.order.back()) / 10;
				int gear_id = (-s.order.back()) % 10;
				double length = lengths[stick_id];
			
				// Get stick endpoints' coordinates
				Vector2d p (x,y), p0 (lastX, lastY);
				double alpha = atan2(p(1) - p0(1), p(0) - p0(0));
				double alpha2 = alpha + (stickDir ? 1 : -1) * asin(length/((p0-p).norm()));
				Vector2d p2 = p0 + Vector2d(cos(alpha2), sin(alpha2)) * lastR;
				Vector2d p2b = p + Vector2d(cos(alpha2), sin(alpha2)) * r;

				// Check for collisions
				for(int o_i = 0; o_i < obs.size(); o_i++) {
					Vector2d po = obs[o_i].topLeftCorner<2,1>();
					Vector2d v21 = (p2b - p2) / (p2b - p2).norm();
					Vector2d vPerp (-v21(1), v21(0));
					double distPara = (po - p2).dot(v21);
					double distPerp = (po - p2).dot(vPerp);
					double oR = obs[o_i](2) + 1.0;
					if((distPara > -oR) && (distPara < (length + oR))) {
						if(abs(distPerp) < oR) {
							fail = true;
							break;
						}
					}
				}

				if(fail) continue;
			}

			// Move
			X(i) = (stickDir ? 1 : -1) * x; Y(i) = y;
			lastY = y, lastX = x, lastR = r;
		}
	
		// Stop if failed
		if(fail) continue;

		// Check if should connect to the last part by placing the last gear in the order
		if(s.check) {

			// Connect with a spur gear 
			double r, r1, r2;
			if(s.order.back() >= 0) {
				r = radii[s.order.back()];
				r1 = lastR + r;
				r2 = r + objn(2);
			}
			else {

				// Get info
				int stick_id = (-s.order.back()) / 10;
				int gear_id = (-s.order.back()) % 10;
				double length = lengths[stick_id];
				r = radii[gear_id];
				r1 = sqrt(sq(length) + sq(lastR - r));
				r2 = r + objn(2);
			}
	
			// Check if the two from the last and the last gears can be connected
			double distSQ = sq(lastX - objn(0)) + sq(lastY - objn(1));
			if(distSQ > sq(r1+r2)) {
				continue;
			}

			// Compute the location of the one from the last gear
			Vector2d p = solveL2Distance(Vector2d(lastX,lastY),objn.topLeftCorner<2,1>(), r1, r2);

			// Check for collisions
			for(int o_i = 0; o_i < obs.size(); o_i++) {
				double distSQ = sq(p(0)-obs[o_i](0)) + sq(p(1)-obs[o_i](1));
				if(distSQ < sq(r+obs[o_i](2))) {
					fail = true;
					break;
				}
				if(fail) continue;
			}

			// Check stick for collisions
			bool stickDir = true;
			if(s.order.back() < 0) {

				// Choose stick direction randomly
				stickDir = (round(((double) rand()) / RAND_MAX) == 0);
				int stick_id = (-s.order.back()) / 10;
				int gear_id = (-s.order.back()) % 10;
				double length = lengths[stick_id];
			
				// Get stick endpoints' coordinates
				Vector2d p0 (lastX, lastY);
				double alpha = atan2(p(1) - p0(1), p(0) - p0(0));
				double alpha2 = alpha + (stickDir ? 1 : -1) * asin(length/((p0-p).norm()));
				Vector2d p2 = p0 + Vector2d(cos(alpha2), sin(alpha2)) * lastR;
				Vector2d p2b = p + Vector2d(cos(alpha2), sin(alpha2)) * r;

				// Check for collisions
				for(int o_i = 0; o_i < obs.size(); o_i++) {
					Vector2d po = obs[o_i].topLeftCorner<2,1>();
					Vector2d v21 = (p2b - p2) / (p2b - p2).norm();
					Vector2d vPerp (-v21(1), v21(0));
					double distPara = (po - p2).dot(v21);
					double distPerp = (po - p2).dot(vPerp);
					double oR = obs[o_i](2) + 1.0;
					if((distPara > -oR) && (distPara < (length + oR))) {
						if(abs(distPerp) < oR) {
							fail = true;
							break;
						}
					}
				}

				if(fail) continue;
			}

			X(X.rows()-1) = (stickDir ? 1 : -1) * p(0);
			Y(X.rows()-1) = p(1);
		}
		
		if(totalIters != NULL) *totalIters = iter;
		return true;
	}

	return false;
}

/* ******************************************************************************************** */
void heuristic (State& s) {

	int totalIters = 0;
	VectorXd X, Y;
	int trials = 10, success = 0;
	for(int i = 0; i < trials; i++) {
		int iters;
		if(feasible(s,X,Y,&iters)) {
			totalIters += iters;
			success++;
		}
	}
	s.H = (((double) success) / totalIters);
}

/* ******************************************************************************************** */
State* applyAction (const State& s, Action_& a) {

	// printf("\nApplying '%s'\n", a.s().c_str());
	assert(!s.check && "Should not reach here");

	// Copy state
	State* s2 = new State;
	s2->order = s.order;
	s2->check = s.check;
	s2->used = s.used;
	s2->used_sticks = s.used_sticks;
	
	// Make the changes
	if(a.use >= 0) {
		assert(!s2->used[a.use]);
		s2->used[a.use] = true;
		s2->order.push_back(a.use);
	}
	else {
		int gear_id = (-a.use) % 10;
		int stick_id = (-a.use) / 10;
		assert(!s2->used[gear_id]);
		s2->used[gear_id] = true;
		s2->used_sticks[stick_id] = true;
		s2->order.push_back(a.use);
	}
	if(a.check) s2->check = true;
		
	// Compute heuristic
	switch (H_TYPE) {
		case 0: heuristic(*s2); break;
		case 1: s2->H = s.H+1; break;
		case 2: s2->H = rand(); break;
	};
	s2->d = s.d+1;

	return s2;
}

/* ******************************************************************************************** */
void printResult (State& s, VectorXd& X, VectorXd& Y) {

	ofstream res ("res");
	res << obj0 << "\n0" << endl;
	for(int i = 0; i < X.rows(); i++) {
		if(s.order[i] < 0) {
			int stick_id = (-s.order[i]) / 10;
			int gear_id = (-s.order[i]) % 10;
			res << -abs(X(i)) << "\n" << Y(i) << "\n" << radii[gear_id] << "\n";
			if(X(i) < 0)
				res << -lengths[stick_id] << endl;
			else 
				res << lengths[stick_id] << endl;
		}
		else 
			res << X(i) << "\n" << Y(i) << "\n" << radii[s.order[i]] << "\n0" << endl;
	}
	res << objn << "\n0" << endl;
	for(int i = 0; i < obs.size(); i++) 
		res << obs[i](0) << "\n" << obs[i](1) << "\n" << -obs[i](2) << "\n0" << endl;
	res.close();
}

/* ******************************************************************************************** */
bool dfs (State* s0) {

	pq queue;
	queue.push(s0);
	int c_ = 0;
	while(!queue.empty()) {

		// Get the current state
		State* s = queue.top();
		// printf("s%d_%d: %s\n", c_, s->d, s->s().c_str());
		queue.pop();
		c_++;

		// Check if constraints are satisfied
		Eigen::VectorXd X, Y;
		if(!feasible(*s, X, Y)) continue;

		// Check if reached goal
		if(s->check) {
			printResult(*s, X,Y);
			printf("%d\n", c_);	
			return true;	
		}

		// Generate children
		vector <Action_> A;
		generateActions(*s, A);
		for(int i = 0; i < A.size(); i++) {
			State* sn = applyAction(*s, A[i]);
			queue.push(sn);
		}
	}
	return false;
}

/* ******************************************************************************************** */
int main (int argc, char* argv[]) {

	// Setup
	srand(time(NULL));
	int bla = rand();
	printf("random seed: %d\n", bla);
	srand(bla);
	// if(argc > 1) srand(atoi(argv[1]));
	if(argc > 1) H_TYPE = atoi(argv[1]);
	if(argc > 2) objn(0) = maxP(0) = atof(argv[2]);
	printf("H_TYPE: %d\n", H_TYPE);
	printf("x dist: %lf\n", objn(0));
		
	State s = initialState();
	// printf("S0: %s\n", s.s().c_str());
	//getchar();
//
//	vector <Action_> A0;
//	generateActions(s, A0);
//	for(int i = 0; i < A0.size(); i++) 
//		printf("%d: %s\n", i, ( A0[i]).s().c_str());
//	//getchar();
//	State* s2 = applyAction(s, A0[0]);
//	printf("S2: %s\n", s2->s().c_str());
//	//getchar();
//
//	generateActions(*s2, A0);
//	for(int i = 0; i < A0.size(); i++) 
//		printf("%d: %s\n", i, ( A0[i]).s().c_str());
//	//getchar();
//	State* s3 = applyAction(*s2, A0[5]);
//	printf("S3: %s\n", s3->s().c_str());
//	//getchar();
////
////	generateActions(*s3, A0);
////	for(int i = 0; i < A0.size(); i++) 
////		printf("%d: %s\n", i, ( A0[i]).s().c_str());
////	// getchar();
////	State* s4 = applyAction(*s3, A0[0]);
////	printf("S4: %s\n", s4->s().c_str());
////	//getchar();
////
////	generateActions(*s4, A0);
////	for(int i = 0; i < A0.size(); i++) 
////		printf("%d: %s\n", i, ( A0[i]).s().c_str());
////	//getchar();
////	State* s5 = applyAction(*s4, A0[1]);
////	printf("S5: %s\n", s5->s().c_str());
//////  getchar();
////
//	VectorXd X, Y;
//	bool res2= feasible(*s3, X,Y);
//	printf("res: %d\n", res2);
//	printResult(*s3, X,Y);

	bool res = dfs(&s);
	if(!res) printf("0\n");
}

