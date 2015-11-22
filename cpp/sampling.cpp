/** 
 * @author Can Erdogan
 * @file sampling.cpp
 * @date 2015-11-22
 */

#include <nlopt.h>
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>

using namespace Eigen;
using namespace std;

/* ******************************************************************************************** */
static const int NUM_OBJECTS = 2;
vector <double> radii;
Vector3d obj0 (0.5, 1.0, 1);
Vector3d objn (7.0, 2.0, 2);

Vector3d obs1 (5.0, 2.0, 5.0);
Vector2d minP (1.0, -3.0), maxP (7.0, 3.0);

/* ******************************************************************************************** */
typedef double(*Function)(VectorXi, VectorXd);
typedef pair <VectorXi, vector <Function> > Constraint;
struct Action_ {
	VectorXi objIds;
	vector <Constraint> Cs; 
	string s() {return "";}
};
/* ******************************************************************************************** */
/// Functions 
double contact_offset = 1e-1;
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
	return (radii[i(0)] + objn(2) - contact_offset) - (objn.topLeftCorner<2,1>()-p).norm();
}
inline double contact_high_on (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (objn.topLeftCorner<2,1>()-p).norm() - (radii[i(0)] + objn(2) + contact_offset);
}
inline double collision1 (VectorXi i, VectorXd x){
	Vector2d p (x(2*i(0)), x(2*i(0)+1));
	return (radii[i(0)] + obs1(2)) - (obs1.topLeftCorner<2,1>()-p).norm();
}
string fs (Function f) {
	if(f == contact_low) return "c_low";
	else if(f == contact_high) return "c_high";
	else if(f == contact_low_o0) return "c_low_o0";
	else if(f == contact_high_o0) return "c_high_o0";
	else if(f == contact_low_on) return "c_low_on";
	else if(f == contact_high_on) return "c_high_on";
	else if(f == collision1) return "c_coll1";
	else assert(false && "unknown f name");
}

/* ******************************************************************************************** */
struct AddGear : Action_ {
	string s() { char buf [256]; sprintf(buf, "(%d,%d) (#%lu):\n", objIds(0), objIds(1), Cs.size()); 
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

struct State {
	vector <int> numContacts;
	vector <Constraint> Cs; 
};

/* ******************************************************************************************** */
State initialState () {
	State s;
	s.numContacts;
	for(int i = 0; i < 2+NUM_OBJECTS; i++) s.numContacts.push_back(0);
	s.numContacts[0] += 2;
	s.numContacts[NUM_OBJECTS+1] += 1;
	s.numContacts[1] = 1;
	return s;
}

/* ******************************************************************************************** */
void generateActions (const State& s, vector <Action_>& acts) {

	for(int i = 0; i < s.numContacts.size(); i++) {
		if(s.numContacts[i] != 0) continue;
		for(int j = 0; j < s.numContacts.size(); j++) {
			if(i == j) continue;
			if(s.numContacts[j] != 1) continue;
			AddGear a;
			a.objIds = Vector2i(i,j);
			if(i == 0 || j == 0) 
				a.Cs.push_back(make_pair(Vector2i(max(i,j),-1), vector <Function> {&contact_low_o0, &contact_high_o0}));
			else if(i == NUM_OBJECTS+1 || j == NUM_OBJECTS+1) 
				a.Cs.push_back(make_pair(Vector2i(min(i,j),-1), vector <Function> {&contact_low_on, &contact_high_on}));
			else 
				a.Cs.push_back(make_pair(Vector2i(i,j), vector <Function> {&contact_low, &contact_high}));
			a.Cs.push_back(make_pair(Vector2i(i,-1), vector <Function> {&collision1}));
			acts.push_back(a);
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
int main (int argc, char* argv[]) {

	// Setup
	srand(time(NULL));
	State s = initialState();
	vector <Action_> A0;
	generateActions(s, A0);
	for(int i = 0; i < A0.size(); i++) 
		printf("%s\n", ((AddGear&) A0[i]).s().c_str());
	
	
//
//	// Generate the minimal simplex
//	vector <Vector2d> H;
//	vector <Vector2i> k;
//	generateInitialHull(C, lims, H, k);
}

