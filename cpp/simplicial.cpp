/** 
 * @author Can Erdogan
 * @file simplicial.cpp
 * @date 2015-11-21
 */

#include <nlopt.h>
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>

using namespace Eigen;
using namespace std;

typedef double(*constraint)(double, double);

/* ******************************************************************************************** */
// Constants 
int DIMS = 2;

/* ******************************************************************************************** */
/// Constraints
double optimOffset = 0.0;
inline double eq1 (double x1, double x2) {
  return (-(-pow(0.2*(-x1-5),4.0) + pow((0.4 * -x1),3) + 1.5*pow((-x1),2) +10) - x2) + optimOffset;
}
inline double eq2 (double x1, double x2) {
	return -(-pow((0.25*(x1)+3),2)+20 - x2) + optimOffset;
}
inline double eq3 (double x1, double x2) {
	return -(-pow((0.2*x1+0.5),5) + 2*pow((0.2*x1+0.5),4) + 3*pow((0.2*x1+0.5),2) - x2) + optimOffset;
}
inline double eq4 (double x1, double x2) {
	return (-pow((0.25*x1+0.5),3)-20 - x2) + optimOffset;
}

/* ******************************************************************************************** */
Vector2d random (const Vector4d& lims) {
	double rx = (((double) rand()) / RAND_MAX) * (lims(1) - lims(0)) + lims(0);
	double ry = (((double) rand()) / RAND_MAX) * (lims(3) - lims(2)) + lims(2);
	return Vector2d(rx,ry);
}

/* ******************************************************************************************** */
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
Vector2d findFeasibleBySampling (constraint C[], const Vector4d& lims) {

	for(int i = 0; i < 1e5; i++) {
		Vector2d p = random(lims);
		bool fail = false;
		for(int j = 0; j < 4; j++) {
			if(C[j](p(0), p(1)) > 0) {
				fail = true;
				break;
			}
		}
		if(!fail) return p;
	}

	assert(false && "Sampling should have found one solution");
	return Vector2d(0,0);
}

/* ******************************************************************************************** */
bool addPointForHull0(const Vector2d& p0, vector <Vector2d> H, constraint C[], 
		const Vector4d& lims, Vector2d& newP) {

	return false;
}

/* ******************************************************************************************** */
void generateInitialHull (constraint C[], const Vector4d& lims, vector <Vector2d> H, 
		vector <Vector2i> k) {

	// Find a feasible sample
	optimOffset = 4.0;
	Vector2d p0 = findFeasibleBySampling(C, lims); 
	optimOffset = 0.0;

  // Choose a random line, find intersections with other inequalities,
  // find the closest one to the solution which we assume is in the
  // feasible space, collect three such points.
	while(H.size() <= DIMS) {

    // Try to find a new point to add to the boundary
		Vector2d newP;
		if(!addPointForHull0(p0, H, C, lims, newP)) continue;

   	// Add the point
		H.push_back(newP);
	}
}

/* ******************************************************************************************** */
int main (int argc, char* argv[]) {

	// Setup
	srand(time(NULL));
	Vector4d lims (-30, 30, -30, 30);
	constraint C [] = {eq1, eq2, eq3, eq4};

	// Generate the minimal simplex
	vector <Vector2d> H;
	vector <Vector2i> k;
	generateInitialHull(C, lims, H, k);
}

/* ******************************************************************************************** *
double nlopt_cost (unsigned n, const double *x, double *grad, void *data) {
	Vector2d* d = (Vector2d*) data;
	if(grad) {
		grad[0] = 2 * (x[0] - (*d)(0)); 
		grad[1] = 2 * (x[1] - (*d)(1)); 
	}
	return pow((x[0] - (*d)(0)),2) + pow((x[1] - (*d)(1)),2);
}

/* ******************************************************************************************** *
double nlopt_bounds (unsigned n, const double *x, double *grad, void *data) {
	int* eqIdx = (int*) data;
	double x1 = x[0];
	switch (*eqIdx) {
		case 1: {
			if(grad) {
				grad[0] = x1*-3.0+pow(x1*(1.0/5.0)+1.0,3.0)*(4.0/5.0)+(x1*x1)*(2.4E1/1.25E2);
				grad[1] = -1.0;
			}
			return eq1(x[0], x[1]);
		} break;
		case 2: {
			if(grad) {
				grad[0] = x1*(1.0/8.0)+3.0/2.0;
				grad[1] = 1.0;
			}
			return eq1(x[0], x[1]);
		} break;
		case 3: {
			if(grad) {
				grad[0] = x1*(-6.0/2.5E1)-pow(x1*(1.0/5.0)+1.0/2.0,3.0)*(8.0/5.0)+pow(x1*(1.0/5.0)+1.0/2.0,4.0)-3.0/5.0;
				grad[1] = 1.0;
			}
			return eq1(x[0], x[1]);
		} break;
		case 4: {
			if(grad) {
				grad[0] = pow(x1*(1.0/4.0)+1.0/2.0,2.0)*(-3.0/4.0);
				grad[1] = -1.0;
			}
			return eq1(x[0], x[1]);
		} break;
	};
}

/* ******************************************************************************************** *
Vector2d nlopt_solve (constraint C[] , const Vector4d& lims, Vector2d& x0) {

	// The optimization engine with the cost function
	nlopt_opt opt = nlopt_create(NLOPT_GN_ISRES, 2);
	double lb[] = {lims(0), lims(2)};
	double ub[] = {lims(1), lims(3)};
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, nlopt_cost, &x0);

	// The two constraints
	int data[] = {1,2,3,4}; 
	nlopt_add_inequality_constraint(opt, nlopt_bounds, &data[0], 1e-4);
	nlopt_add_inequality_constraint(opt, nlopt_bounds, &data[1], 1e-4);
	nlopt_add_inequality_constraint(opt, nlopt_bounds, &data[2], 1e-4);
	nlopt_add_inequality_constraint(opt, nlopt_bounds, &data[3], 1e-4);

	// The stopping criteria
	nlopt_set_xtol_rel(opt, 1e-4);

	// Optimize!
	double minf; 
	double x[2] = {x0(0), x0(1) }; 
	int res = nlopt_optimize(opt, x, &minf);
	assert(res > 0 && "Weird optimization failure");
	nlopt_destroy(opt);

	// Check if the value is in the feasible region 
	bool feasible = true;
	for(int i = 0; i < 4; i++) {
		double val = C[i](x[0], x[1]);
		if(val > 0) {
			printf("val: %lf (for %d)\n", val, i);
			feasible = false;
			break;
		}
	}
	assert(feasible && "Optimization value is not in feasible region.");
	return Vector2d(x[0], x[1]);
}
/* ******************************************************************************************** */


