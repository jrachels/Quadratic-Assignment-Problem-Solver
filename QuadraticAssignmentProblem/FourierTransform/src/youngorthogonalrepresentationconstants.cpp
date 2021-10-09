#include <cmath>
#include <cassert>

#include "QuadraticAssignmentProblem/FourierTransform/src/youngorthogonalrepresentationconstants.h"

QAPSolver::YoungOrthogonalRepresentationConstants& QAPSolver::YoungOrthogonalRepresentationConstants::GrowYORConstants(int i) {
	int current_size = constants_.size();
	if (current_size < i) {
		constants_.resize(i);
		if (current_size == 0) {
			std::get<0>(constants_[0]) = std::numeric_limits<double>::infinity();
			std::get<1>(constants_[0]) = std::numeric_limits<double>::infinity();
			current_size += 1;
		}
	}
	while (current_size < i) {
		double current_size_d = static_cast<double>(current_size);
		std::get<0>(constants_[current_size]) = 1.0 / current_size_d;
		std::get<1>(constants_[current_size]) = std::sqrt(current_size_d * current_size_d - 1.0) / current_size_d;
		current_size += 1;
	}
	assert(current_size >= i);
	return *this;
}

double QAPSolver::YoungOrthogonalRepresentationConstants::current_size() {
	return constants_.size();
}
