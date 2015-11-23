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

/* ******************************************************************************************** */
static const int NUM_GEARS = 4;
static const int NUM_OBJECTS = NUM_GEARS;
vector <double> radii {1.4, 1.5, 1.2, 1.5}; //, 1.6};
vector <double> lengths {1.2}; //, 1.6};
Vector3d obj0 (0.5, 1.0, 1);
Vector3d objn (14.0, 2.0, 2);

Vector3d obs1 (20.0, 0.0, -0.2);
Vector3d obs2 (20.0, 2.0, -0.2);
Vector2d minP (1.0, -1.0), maxP (objn(0), 3.0);

/* ******************************************************************************************** */
typedef double(*Function)(VectorXi, VectorXd);
typedef pair <VectorXi, vector <Function> > Constraint;

/* ******************************************************************************************** */
/// Functions 
double contact_offset = 1e-1;
double length_offset = 1;
inline double contact_low (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1));
	return (radii[i(0)] + radii[i(1)] - contact_offset) - (p1-p2).norm();
}
inline double contact_high (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1));
	return  (p1-p2).norm() - (radii[i(0)] + radii[i(1)] + contact_offset);
}
inline double contact_low_o0 (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (radii[i(0)] + obj0(2) - contact_offset) - (obj0.topLeftCorner<2,1>()-p).norm();
}
inline double contact_high_o0 (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (obj0.topLeftCorner<2,1>()-p).norm() - (radii[i(0)] + obj0(2) + contact_offset);
}
inline double contact_low_on (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	//cout << "p: " << p.transpose() << " vs. objn: " <<  objn.topLeftCorner<2,1>().transpose() << endl;
	//printf("dist: %lf\n", (radii[i(0)] + objn(2) - contact_offset));
	return (radii[i(0)] + objn(2) - contact_offset) - (objn.topLeftCorner<2,1>()-p).norm();
}
inline double contact_high_on (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (objn.topLeftCorner<2,1>()-p).norm() - (radii[i(0)] + objn(2) + contact_offset);
}
inline double collision1 (VectorXi i, VectorXd x){
	return -1;
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (radii[i(0)] + -obs1(2)) - (obs1.topLeftCorner<2,1>()-p).norm();
}
inline double collision2 (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (radii[i(0)] + -obs2(2)) - (obs2.topLeftCorner<2,1>()-p).norm();
}
inline double length_low (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1));
	return (lengths[(i(0)-2*NUM_GEARS)/4] - length_offset) - (p1-p2).norm();
}
inline double length_high (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1));
	return  (p1-p2).norm() - (lengths[(i(0)-2*NUM_GEARS)/4] + length_offset);
}
double perp_offset = 1.5;
inline double contact_perp_low (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1)), p3 (x(2*i(2)), x(2*i(2)+1));;
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d v31 = (p3 - p1) / (p3-p1).norm();
	return  - perp_offset - v21.dot(v31);
}
inline double contact_perp_high (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1)), p3 (x(2*i(2)), x(2*i(2)+1));;
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d v31 = (p3 - p1) / (p3-p1).norm();
	return v21.dot(v31) - perp_offset;
}
inline double contact_perp_o0_low (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1)), p3 = obj0.topLeftCorner<2,1>();
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d v31 = (p3 - p1) / (p3-p1).norm();
	return  - perp_offset - v21.dot(v31);
}
inline double contact_perp_o0_high (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1)), p3 = obj0.topLeftCorner<2,1>();
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d v31 = (p3 - p1) / (p3-p1).norm();
	return v21.dot(v31) - perp_offset;
}
inline double contact_perp_on_low (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1)), p3 = objn.topLeftCorner<2,1>();
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d v31 = (p3 - p1) / (p3-p1).norm();
	printf("v21 . v31: %lf\n", v21.dot(v31));
	return  - perp_offset - v21.dot(v31);
}
inline double contact_perp_on_high (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1)), p3 = objn.topLeftCorner<2,1>();
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d v31 = (p3 - p1) / (p3-p1).norm();
	return v21.dot(v31) - perp_offset;
}
inline double collision_stick1 (VectorXi i, VectorXd x){
	Vector2d p1 (x(2*i(0)), x(2*i(0)+1)), p2 (x(2*i(1)), x(2*i(1)+1));
	double length = lengths[(i(0)-2*NUM_GEARS)/4];
	Vector2d v21 = (p2 - p1) / (p2-p1).norm();
	Vector2d po = obs1.topLeftCorner<2,1>();
	double paraDist = v21.dot((po-p1));	
	if(!((paraDist > 0) && (paraDist < length))) return 0;
	Vector2d vPerp (-v21(1), v21(0));
	double perpDist = abs(vPerp.dot((po-p1)));
	if(perpDist < obs1(2)) return (perpDist - obs1(2));
	return 0.0;
}
string fs (Function f) {
	if(f == contact_low) return "c_low";
	else if(f == contact_high) return "c_high";
	else if(f == contact_low_o0) return "c_low_o0";
	else if(f == contact_high_o0) return "c_high_o0";
	else if(f == contact_low_on) return "c_low_on";
	else if(f == contact_high_on) return "c_high_on";
	else if(f == length_low) return "len_low";
	else if(f == length_high) return "len_high";
	else if(f == collision1) return "coll1";
	else if(f == collision2) return "coll2";
	else if(f == contact_perp_low) return "p_low";
	else if(f == contact_perp_high) return "p_high";
	else if(f == contact_perp_o0_low) return "p_low_o0";
	else if(f == contact_perp_o0_high) return "p_high_o0";
	else if(f == contact_perp_on_low) return "p_low_on";
	else if(f == contact_perp_on_high) return "p_high_on";
	else if(f == collision_stick1) return "coll_s1";
	else assert(false && "unknown f name");
}

/* ******************************************************************************************** */
struct Action_ {
	int type;
	VectorXi objIds;
	vector <Constraint> Cs; 
	Action_() { type = -1;}
	string s() { 
		char buf [1024]; 
		stringstream ss;
		ss << objIds.transpose();
		sprintf(buf, "(%s) (t: %d, #%lu):\n", ss.str().c_str(), type, Cs.size()); 
		for(int i = 0; i < Cs.size(); i++) {
			stringstream ss;
			ss << Cs[i].first.transpose();	
			sprintf(buf, "%s    %s:\n", buf, ss.str().c_str());
			for(int j = 0; j < Cs[i].second.size(); j++) {
				sprintf(buf, "%s        %s\n", buf, fs(Cs[i].second[j]).c_str());
			}
		}
		return buf;
	}
};

struct AddGear : Action_ {
	AddGear (int i, int j) {
		type = 0;
		objIds = Vector2i(i,j);
		if(i == 0 || j == 0) 
			Cs.push_back(make_pair(Vector2i(max(i,j)-1,-1), vector <Function> {&contact_low_o0, &contact_high_o0}));
		else if(i == NUM_GEARS+1 || j == NUM_GEARS+1) 
			Cs.push_back(make_pair(Vector2i(min(i-1,j-1),-1), vector <Function> {&contact_low_on, &contact_high_on}));
		else 
			Cs.push_back(make_pair(Vector2i(i-1,j-1), vector <Function> {&contact_low, &contact_high}));
		Cs.push_back(make_pair(Vector2i(i-1,-1), vector <Function> {&collision1}));
		Cs.push_back(make_pair(Vector2i(i-1,-1), vector <Function> {&collision2}));
	}
};

struct ConnectGear : Action_ {
	ConnectGear (int i, int j) {
		type = 1;
		objIds = Vector2i(i,j);
		if(i == 0 || j == 0) 
			Cs.push_back(make_pair(Vector2i(max(i,j)-1,-1), vector <Function> {&contact_low_o0, &contact_high_o0}));
		else if(i == NUM_GEARS+1 || j == NUM_GEARS+1) 
			Cs.push_back(make_pair(Vector2i(min(i,j)-1,-1), vector <Function> {&contact_low_on, &contact_high_on}));
		else 
			Cs.push_back(make_pair(Vector2i(i-1,j-1), vector <Function> {&contact_low, &contact_high}));
	}
};

struct AddWorm : Action_ {
	AddWorm (int i, int j, int k) {
		type = 2;
		objIds = Vector3i(i,j,k);		// assume i < j 
		int w_id = k - (NUM_GEARS+2);  	// 0 if first
		int x_id = (NUM_GEARS*2 + w_id*4)/2, x_id2 = x_id+1;
		if(i == 0) {
			Cs.push_back(make_pair(Vector2i(x_id,-1), vector <Function> {&contact_low_o0, &contact_high_o0}));
			Cs.push_back(make_pair(Vector2i(x_id,x_id2), vector <Function> {&contact_perp_o0_low, &contact_perp_o0_high}));
			Cs.push_back(make_pair(Vector2i(x_id2,j-1), vector <Function> {&contact_low, &contact_high}));
			Cs.push_back(make_pair(Vector3i(x_id,x_id2,j-1), vector <Function> {&contact_perp_low, &contact_perp_high}));
		}
		else if(j == NUM_GEARS+1) {
			Cs.push_back(make_pair(Vector2i(x_id,i-1), vector <Function> {&contact_low, &contact_high}));
//			Cs.push_back(make_pair(Vector3i(x_id,x_id2,i-1), vector <Function> {&contact_perp_low, &contact_perp_high}));
			Cs.push_back(make_pair(Vector2i(x_id2,-1), vector <Function> {&contact_low_on, &contact_high_on}));
			Cs.push_back(make_pair(Vector2i(x_id,x_id2), vector <Function> {&contact_perp_on_low, &contact_perp_on_high}));
		}
		else {
			Cs.push_back(make_pair(Vector2i(x_id,i-1), vector <Function> {&contact_low, &contact_high}));
			Cs.push_back(make_pair(Vector3i(x_id,x_id2,i-1), vector <Function> {&contact_perp_low, &contact_perp_high}));
			Cs.push_back(make_pair(Vector2i(x_id2,j-1), vector <Function> {&contact_low, &contact_high}));
			Cs.push_back(make_pair(Vector3i(x_id,x_id2,j-1), vector <Function> {&contact_perp_low, &contact_perp_high}));
		}
		Cs.push_back(make_pair(Vector2i(x_id,x_id2), vector <Function> {&length_low, &length_high, &collision_stick1}));
//		if(j == 0) 
//			Cs.push_back(make_pair(Vector2i(max(i,j)-1,-1), vector <Function> {&contact_low_o0, &contact_high_o0}));
//		else if(i == NUM_OBJECTS+1 || j == NUM_OBJECTS+1) 
//			Cs.push_back(make_pair(Vector2i(min(i,j)-1,-1), vector <Function> {&contact_low_on, &contact_high_on}));
//		else 
//			Cs.push_back(make_pair(Vector2i(i-1,j-1), vector <Function> {&contact_low, &contact_high}));
	}
};

struct State {
	vector <int> numContacts;
	vector <Constraint> Cs; 
	double H;
	int d;
	string s () {
		char buf [1024];
		sprintf(buf, "contacts: {");
		for(int i = 0; i < numContacts.size(); i++) sprintf(buf, "%s%d ", buf, numContacts[i]);
		sprintf(buf, "%s\b}\n constraints: \n", buf);
		for(int i = 0; i < Cs.size(); i++) {
			stringstream ss;
			ss << Cs[i].first.transpose();	
			sprintf(buf, "%s    %s:\n", buf, ss.str().c_str());
			for(int j = 0; j < Cs[i].second.size(); j++) {
				sprintf(buf, "%s        %s\n", buf, fs(Cs[i].second[j]).c_str());
			}
		}
		return buf;
	}
};

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
	s.numContacts;
	s.d = 0;
	for(int i = 0; i < 2+NUM_OBJECTS; i++) s.numContacts.push_back(0);
	s.numContacts[0] += 1;
	s.numContacts[NUM_GEARS+1] += 1;
	s.H = 0;
//	s.numContacts[1] = 1;
	return s;
}

/* ******************************************************************************************** */
void heuristic (State& s) {
}

/* ******************************************************************************************** */
State* applyAction (const State& s, const Action_& a) {

	// Copy state
	State* s2 = new State;
	s2->numContacts = s.numContacts;
	s2->Cs = s.Cs;
	
	// Add the new constraints
	for(int i = 0; i < a.Cs.size(); i++) s2->Cs.push_back(a.Cs[i]);

	// Add gear type
	if(a.type == 0) {
		s2->numContacts[a.objIds(0)] = 1;
		s2->numContacts[a.objIds(1)] = 2;
	}
	else if(a.type == 1) {
		s2->numContacts[a.objIds(0)] = 2;
		s2->numContacts[a.objIds(1)] = 2;
	}
	else if(a.type == 2) {
		s2->numContacts[a.objIds(0)] = 2;
		s2->numContacts[a.objIds(1)] = 2;
		s2->numContacts[a.objIds(2)] = 2;
	}

	// Compute heuristic
	heuristic(*s2);
	s2->H = s.H-1;
	s2->d = s.d+1;

	return s2;
}

/* ******************************************************************************************** */
void generateActions (const State& s, vector <Action_>& acts) {

	// See if two gears can be connected with a worm gear
	acts.clear();
	for(int i = 0; i < NUM_GEARS+2; i++) {
		if(s.numContacts[i] != 1) continue;		
		for(int j = i+1; j < NUM_GEARS+2; j++) {
			if(s.numContacts[j] != 1) continue;		
			for(int k = NUM_GEARS+2; k < s.numContacts.size(); k++) {
				if(s.numContacts[k] != 0) continue;
				if((i == 0) && (j == NUM_GEARS+1)) continue;
				AddWorm a (i,j,k);
				acts.push_back(a);
			}
		}
	}

	// See what can be done with each object
	for(int i = 0; i < NUM_GEARS+2; i++) {
		
		// If already connected to two gears, nothing
		if(s.numContacts[i] >= 2) continue;		

		// If exists, connect to one of the gears
		if(s.numContacts[i] == 1) {
			for(int j = i+1; j < NUM_GEARS+2; j++) {
				if(s.numContacts[j] != 1) continue;
				if((i == 0) && (j == NUM_GEARS+1)) continue;
				ConnectGear a (i,j);
				acts.push_back(a);
			}
		}

		// Add this gear to connect to an existing gear
		if(s.numContacts[i] == 0) {
			for(int j = 0; j < NUM_GEARS+2; j++) {
				if(i == j) continue;
				if(s.numContacts[j] != 1) continue;
				AddGear a (i,j);
				acts.push_back(a);
			}
		}
	}
}

/* ******************************************************************************************** */
Vector2d random (const Vector4d& lims) {
	double rx = (((double) rand()) / RAND_MAX) * (lims(1) - lims(0)) + lims(0);
	double ry = (((double) rand()) / RAND_MAX) * (lims(3) - lims(2)) + lims(2);
	return Vector2d(rx,ry);
}

/* ******************************************************************************************** *
void computeFeasibleRatio (constraint C[], const Vector4d& lims) {

	int p_ = 0, n_ = 0;
	for(int i = 0; i < 1e5; i++) {
		Vector2d p = random(lims);
		bool fail = false;
		for(int j = 0; j < 4; j++) {
			if(C[j](p(0), p(1)) > 0) {
				fail = true;
				break;
			}
		}
		if(!fail) p_++;
		else n_++;
	}
	printf("Ratio: %lf\n", ((double) p_) / (p_ + n_));
}

/* ******************************************************************************************** */
bool feasible (vector <Constraint>& Cs, VectorXd& X) {

	X = VectorXd (NUM_GEARS * 2 + (NUM_OBJECTS - NUM_GEARS) * 4);
	for(int i = 0; i < 1e9; i++) {

		// Generate random value
		for(int j = 0; j < (NUM_GEARS + 2 * (NUM_OBJECTS - NUM_GEARS)); j++) {
			X(2*j) = (((double) rand()) / RAND_MAX) * (maxP(0) - minP(0)) + minP(0);
			X(2*j+1) = (((double) rand()) / RAND_MAX) * (maxP(1) - minP(1)) + minP(1);
		}

		// Check against each constraint
		bool fail = false;
		for(int j = 0; j < Cs.size(); j++) {
			VectorXi& ids = Cs[j].first;
			for(int k = 0; k < Cs[j].second.size(); k++) {
				double val = Cs[j].second[k](ids, X);
				if(val > 0) { // printf("(%d,%d): %lf\n", j, k, val); 
fail = true; break;}
			}
			if(fail) break;
		}

		// Return positive if all the constraints passed
		if(!fail) return true;
	}

	return false;
}

/* ******************************************************************************************** */
bool dfs (State* s0) {

	pq queue;
	queue.push(s0);
	int c_ = 0;
	while(!queue.empty()) {

		// Get the current state
		State* s = queue.top();
		printf("s%d_%d: %s\n", c_, s->d, s->s().c_str());
		queue.pop();
		c_++;

		// Check if constraints are satisfied
		Eigen::VectorXd X;
		if(!feasible(s->Cs, X)) continue;

		// Check if reached goal
		bool all2s = true;
		for(int i = 0; i < s->numContacts.size(); i++) {
			if(s->numContacts[i] == 1) { all2s = false; break; }
		}
		if(all2s) {
			ofstream res ("res");
			res << obj0 << endl;
			for(int i = 0; i < NUM_OBJECTS; i++)
				res << X.block<2,1>(2*i,0) << "\n" << radii[i] <<  endl;
			res << objn << endl;
			res << obs1 << "\n" << obs2 << endl;
			res.close();
			printf("Reached goal!\n");	
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
	State s = initialState();
	printf("S0: %s\n", s.s().c_str());
	//getchar();

 vector <Action_> A0;
 generateActions(s, A0);
 for(int i = 0; i < A0.size(); i++) 
 	printf("%s\n", ( A0[i]).s().c_str());
 //getchar();
 State* s2 = applyAction(s, A0[0]);
 printf("S2: %s\n", s2->s().c_str());
 //getchar();

 generateActions(*s2, A0);
 for(int i = 0; i < A0.size(); i++) 
 	printf("%s\n", ( A0[i]).s().c_str());
 //getchar();
 State* s3 = applyAction(*s2, A0[1]);
 printf("S3: %s\n", s3->s().c_str());
 //getchar();

 generateActions(*s3, A0);
 for(int i = 0; i < A0.size(); i++) 
 	printf("%s\n", ( A0[i]).s().c_str());
 // getchar();
 State* s4 = applyAction(*s3, A0[1]);
 printf("S4: %s\n", s4->s().c_str());
 //getchar();

generateActions(*s4, A0);
 for(int i = 0; i < A0.size(); i++) 
 	printf("%s\n", ( A0[i]).s().c_str());
 //getchar();
 State* s5 = applyAction(*s4, A0[1]);
 printf("S5: %s\n", s5->s().c_str());
//  getchar();

generateActions(*s5, A0);
 for(int i = 0; i < A0.size(); i++) 
 	printf("%s\n", ( A0[i]).s().c_str());
 //getchar();
 State* s6 = applyAction(*s5, A0[0]);
 printf("S6: %s\n", s6->s().c_str());
  getchar();





 VectorXd X;
 bool res2= feasible(s6->Cs, X);
 printf("res: %d\n", res2);
 cout << X << endl;

 ofstream res ("res");
 res << obj0 << "\n-12" << endl;
 for(int i = 0; i < NUM_GEARS; i++)
 	res << X.block<2,1>(2*i,0) << "\n" << radii[i] << "\n-12" << endl;
 for(int i = 2*NUM_GEARS; i < X.rows(); i++)
 	res << X(i) << endl;
 res << objn << "\n-12" << endl;
 res << obs1 << "\n-12" << "\n" << obs2 << "\n-12" << endl;
 res.close();

//
//	return 1;
//	State s2 = applyAction(s, A0[2]);
//	printf("S2: %s\n", s2.s().c_str());
	
//	vector <Action_> A2;
//	generateActions(s2, A2);
//	for(int i = 0; i < A2.size(); i++) 
//		printf("%s\n", ((AddGear&) A2[i]).s().c_str());
	
//	bool res = dfs(&s);
//	if(!res) printf("Failure\n");
//
//	// Generate the minimal simplex
//	vector <Vector2d> H;
//	vector <Vector2i> k;
//	generateInitialHull(C, lims, H, k);
}

