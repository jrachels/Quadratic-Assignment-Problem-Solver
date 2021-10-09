#include "QuadraticAssignmentProblem/Permutations/permutation.h"
#include "QuadraticAssignmentProblem/Permutations/permutationbuilder.h"

namespace QAPSolver {
	Permutation::Permutation(const PermutationBuilder& left_rotations, int n) : permutation_(n) {
		permutation_.shrink_to_fit();
		for (int i = 0; i < n; i++) {
			permutation_[i] = i + 1;
		}
		int num_rotations = left_rotations.Size();

		for (int i = num_rotations-1; i > -1; i--) {
			// save locally?
			const int start_index = left_rotations[i][0] - 1;
			const int end_index = left_rotations[i][1] - 1;
			// check bounds
			// check order
			// i <= j <= n
			if ((start_index <= end_index) && (end_index <= (n - 1))) {
				// do rotations of perm
				// delegate to private function?
				int temp = permutation_[start_index];
				for (int i = start_index; i < end_index; ++i) {
					permutation_[i] = permutation_[i + 1];
				}
				permutation_[end_index] = temp;
			}
			else {
				throw BuildPermutationError(start_index, end_index, n);
			}
		}
	}
}