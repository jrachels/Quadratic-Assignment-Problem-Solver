// TO DO: clean variable names. add underscores.
// get rid of the private constructor? just save the variable size?

#ifndef QAP_SVD_SQUARESVD_H_
#define QAP_SVD_SQUARESVD_H_

#include <vector>
#include <algorithm>

#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"
#include "QuadraticAssignmentProblem/SVD/Jacobi.h"
#include "QuadraticAssignmentProblem/SVD/realsvd2x2.h"

namespace QAPSolver {
	namespace internal {
		// an SVD algorithm specialized for small, square
		// matrices based on the Jacobi SVD algorithm.
		class SquareSVD
		{
		public:
			explicit SquareSVD(const SquareMatrix& m_workMatrix) : SquareSVD(m_workMatrix, m_workMatrix.size()) {}

			double TraceNorm();

		private:
			// SVD Algorithm based on Eigen's Jacobi SVD
			// I believe the copy here is necessary
			SquareSVD(SquareMatrix m_workMatrix, int size);

			std::vector<double> m_singularValues_;
		};

	}
}


#endif

