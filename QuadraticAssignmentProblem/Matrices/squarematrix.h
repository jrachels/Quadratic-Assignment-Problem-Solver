#ifndef QAP_MATRICES_SQUAREMATRIX_H_
#define QAP_MATRICES_SQUAREMATRIX_H_

#include "QuadraticAssignmentProblem/Matrices/matrix.h"

namespace QAPSolver {
	namespace internal {
		class JacobiRotation;
		//class SVD;
	}
	class SquareMatrix : public Matrix
	{
	public:
		//SquareMatrix(Matrix)

		explicit SquareMatrix(int n) : Matrix(n, n) {};

		SquareMatrix() = default;

		static SquareMatrix IdentityMatrix(int n) {
			SquareMatrix matrix(n);
			for (int i = 0; i < n; ++i) {
				matrix(i, i) = 1;
			}
		}

		double MaxAbsOfDiagonal() const;

		double TraceNorm() const;

		int size() const {
			return this->size1();
		}

		inline void applyOnTheLeft(int p, int q, const internal::JacobiRotation& j);

		inline void applyOnTheRight(int p, int q, const internal::JacobiRotation& j);
	};
}

#endif