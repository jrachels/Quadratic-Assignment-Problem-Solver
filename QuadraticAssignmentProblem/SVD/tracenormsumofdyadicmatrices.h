#ifndef QAP_SVD_TRACENORMSUMOFDYADICMATRIX_H_
#define QAP_SVD_TRACENORMSUMOFDYADICMATRIX_H_

#include "QuadraticAssignmentProblem/Matrices/matrix.h"

namespace QAPSolver {
	namespace internal {
		double TraceNormSumOfDyadicMatrices(const Matrix& left_matrix, const Matrix& right_matrix);
	}
}

#endif