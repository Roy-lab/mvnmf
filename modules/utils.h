#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>
#ifndef _utils_
#define _utils_
using namespace std;

namespace utils
{
	double get_frobenius_norm(gsl_matrix*);
	int get_neighborhood_graph(gsl_matrix*, int);
	int get_laplacian(gsl_matrix*, gsl_matrix*);
	int get_inverse(gsl_matrix*);
};
#endif
