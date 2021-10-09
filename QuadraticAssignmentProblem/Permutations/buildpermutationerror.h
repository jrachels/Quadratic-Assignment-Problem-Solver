#ifndef QAP_PERMUTATIONS_BUILDPERMUTATIONERROR_H_
#define QAP_PERMUTATIONS_BUILDPERMUTATIONERROR_H_

#include <stdexcept>
#include <string>

namespace QAPSolver {
	// error checking is done when permutation is build for speed
	// in the QAPSolver algorithm. 
	// If errors are expected, PermutationBuilder must keep track
	// of the permutation size and error checks should be done in AddRotation.
	class BuildPermutationError : public std::runtime_error {
	public:
		const int start_index_;
		const int end_index_;
		const int permutation_length_;

		BuildPermutationError(int i, int j, int n) : start_index_(i), end_index_(j), permutation_length_(n),
			runtime_error("Error building permutation. Permutation has length " + std::to_string(n) +
				". Start index is " + std::to_string(i) + ". End index is " + std::to_string(j)) {};

	};
}

#endif