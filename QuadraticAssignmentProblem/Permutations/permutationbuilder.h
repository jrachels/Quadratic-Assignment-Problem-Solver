#ifndef QAP_PERMUTATIONS_PERMUTATIONBUILDER_H_
#define QAP_PERMUTATIONS_PERMUTATIONBUILDER_H_

#include <vector>
#include <array>

namespace QAPSolver {
	// currently only computes permutation as a product of rotations from the left,
	// because that is all that is needed in the QAPSolver program.
	class PermutationBuilder
	{
	public:
		PermutationBuilder() = default;

		// Doesn't check for bounds errors. Left to Permutation.
		void AddRotation(int i, int j) {
			left_rotations_.push_back({ i, j });
		}

		std::array<int, 2> operator[](int i) const {
			return left_rotations_[i];
		}

		int Size() const {
			return left_rotations_.size();
		}

	private:
		// each entry (i, j) is a rotation of elements i through j, inclusive
		// permutation is product of rotations.

		// note: For QAPSolver, storing the second integer is unnecessary as we already know what it is
		std::vector<std::array<int, 2>> left_rotations_;
	};
}

#endif