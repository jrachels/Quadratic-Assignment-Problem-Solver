// TO DO: This was made for a slightly different kind of rotation than what QAPSolver uses. Might could optimize.

#ifndef QAP_PERMUTATIONS_PERMUTATION_H_
#define QAP_PERMUTATIONS_PERMUTATION_H_

#include <vector>

#include "QuadraticAssignmentProblem/Permutations/permutationbuilder.h"
#include "QuadraticAssignmentProblem/Permutations/buildpermutationerror.h"

namespace QAPSolver {
	//DON'T COMPUTE PERMUTATION. STORE PERMUTATION AS A LIST OF CYCLES
	class Permutation
	{
	public:
		Permutation(const PermutationBuilder& left_rotations, int n);

		[[nodiscard]] std::vector<int> perm_to_vec() const {
			return permutation_;
		}

		int operator[](int i) const {
			return permutation_[i];
		}

	private:
		std::vector<int> permutation_;
	};
}

#endif