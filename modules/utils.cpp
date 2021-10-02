#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <list>
#include <vector>
#include <math.h>
#include "utils.h"

int utils::get_laplacian(gsl_matrix* W, gsl_matrix* L) {
	int n = W->size1;
	gsl_vector* D = gsl_vector_alloc(n);
	for (int i = 0; i < n; i++) {
		double *val = &(D->data[i]);	
		gsl_vector_view row = gsl_matrix_row(W, i);
		*val = gsl_blas_dasum(&row.vector);
	}
	gsl_vector_view diag = gsl_matrix_diagonal(L);
	gsl_vector_memcpy(&diag.vector, D);
	gsl_matrix_sub(L, W);
	gsl_vector_free(D);
	return 0;
}

int utils::get_neighborhood_graph(gsl_matrix* W, int neighborhood_radius) {
	gsl_vector_view diag = gsl_matrix_diagonal(W);
	gsl_vector_set_all(&diag.vector, 1);
	for (int i = 1; i <= neighborhood_radius; i++) {
		gsl_vector_view lower = gsl_matrix_subdiagonal(W, i);
		gsl_vector_set_all(&lower.vector, 1);
		gsl_vector_view upper = gsl_matrix_superdiagonal(W, i);
		gsl_vector_set_all(&upper.vector, 1);
	}
	return 0;
}

double utils::get_frobenius_norm(gsl_matrix* X) {
	int n = X->size1;
	int m = X->size2;
	double sum = 0;
	gsl_vector* row_copy = gsl_vector_alloc(m);
	for (int i = 0; i < n; i++) {
		gsl_vector_view row = gsl_matrix_row(X, i);
		gsl_vector_memcpy(row_copy, &row.vector); 
		gsl_vector_mul(row_copy, &row.vector);
		double rowsum = gsl_blas_dasum(row_copy);
		sum += rowsum;
	}
	gsl_vector_free(row_copy);
	return sum;
}

int utils::get_inverse(gsl_matrix* A) {
	gsl_linalg_cholesky_decomp1(A);
	gsl_linalg_cholesky_invert(A);
	return 0;
}
