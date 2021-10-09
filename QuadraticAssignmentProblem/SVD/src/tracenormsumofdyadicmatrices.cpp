
#include <cassert>


#include "QuadraticAssignmentProblem/SVD/tracenormsumofdyadicmatrices.h"
#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"
#include "QuadraticAssignmentProblem/SVD/colpivhouseholderqr.h"

// only works for square

double QAPSolver::internal::TraceNormSumOfDyadicMatrices(const Matrix& left_matrix, const Matrix& right_matrix) {
	// check that matrices  are compatible sizes
	int size = left_matrix.size1();
	int num_components = left_matrix.size2();

	assert((size == right_matrix.size2()) && (num_components == right_matrix.size1()));

	SquareMatrix product_matrix(size);

	if (size > num_components) {
		Matrix left_matrix_qr = left_matrix;
		std::vector<int> left_qr_column_order = ColPivHouseHolderQR(left_matrix_qr);
		assert(left_qr_column_order.size() == num_components);
		Matrix right_matrix_qr = right_matrix.Transpose();
		std::vector<int> right_qr_row_order = ColPivHouseHolderQR(right_matrix_qr);

		// compute product_matrix
		for (int i = 0; i < num_components; ++i) {
			int column_in_left_r = left_qr_column_order[i]; // 
			int row_in_right_r = right_qr_row_order[i];
			for (int j = 0; j <= column_in_left_r; ++j) {
				for (int k = 0; k <= row_in_right_r; ++k) {
					product_matrix(j, k) += left_matrix_qr(j, column_in_left_r) * right_matrix_qr(k, row_in_right_r);
				}
			}
		}
	}
	else {
		//SquareMatrix product_matrix(size);
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				double entry = 0;
				for (int k = 0; k < num_components; k++) {
					entry += left_matrix(i, k) * right_matrix(k, j);
				}
				product_matrix(i, j) = entry;
			}
		}
	}

	return product_matrix.TraceNorm();

}