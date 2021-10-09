#include "QuadraticAssignmentProblem/Matrices/squarematrix.h"
#include "QuadraticAssignmentProblem/SVD/squaresvd.h"

double QAPSolver::SquareMatrix::TraceNorm() const {
	return internal::SquareSVD(*this).TraceNorm();
}

double QAPSolver::SquareMatrix::MaxAbsOfDiagonal() const {
	double max = 0;
	int size = this->size1();
	for (int i = 0; i < size; i++) {
		double abs_of_entry = std::abs((*this)(i, i));
		if (abs_of_entry > max) {
			max = abs_of_entry;
		}
	}
	return max;
}